// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file    SimpleGeomGenerate.cxx
/// \author  Maja Kabus
/// \brief   Generating root simple geometry files for event display
/// \remarks Based on simple_geom_generate macro by Jeremi Niedziela on 09/02/2016.

#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationBase/GeometryManager.h"
#include "DetectorsBase/GeometryManager.h"
#include "CCDB/Manager.h"

#include <TApplication.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEveGeoNode.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TRegexp.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct Options {
  std::string detName;   // -m "TPC"
  std::string outputDir; // -d "O2/EventVisualisation/resources/geometry/run3"
  int runNumber;         // -r 3
};

void addNodes(TGeoNode* node, TEveGeoNode* parent, Int_t depth, Int_t depthmax, TObjArray* list)
{
  if (--depth <= 0)
    return;

  TObjArray* nlist = node->GetVolume()->GetNodes(); // all nodes in current level
  if (!nlist)
    return;

  TObjString* nname = (TObjString*)list->At(depthmax - depth); // name of required node in current level

  for (int i = 0; i < nlist->GetEntries(); i++) {
    // loop over nodes in current level and find the one with matching name
    TGeoNode* node2 = (TGeoNode*)nlist->At(i);

    if (strcmp(node2->GetName(), nname->GetString().Data()) == 0) {
      TEveGeoNode* son = dynamic_cast<TEveGeoNode*>(parent->FindChild(nname->GetName()));
      if (!son) {
        son = new TEveGeoNode(node2);
        parent->AddElement(son);
      }
      addNodes(node2, son, depth, depthmax, list);
    }
  }
}

void generateSimpleGeometry(const char* detectorName = "", const int runNumber = 0, const std::string outputDir = ".")
{
  if (strcmp(detectorName, "") == 0) {
    std::cout << "Give name of the detector as a first argument!" << std::endl;
    return;
  }

  // read all files with names matching "geom_list_XYZ.txt"
  std::vector<std::string> detectorsList;
  const char* o2basePath = gSystem->Getenv("ALIBUILD_WORK_DIR");
  const char* o2Path = Form("%s/../O2/EventVisualisation/Geometries/resources", o2basePath);
  TString inPath = Form("%s/geom_list_", o2Path);
  TSystemDirectory dir(o2Path, o2Path);
  TList* files = dir.GetListOfFiles();

  if (files) {
    TRegexp e("geom_list_[A-Z][A-Z][A-Z].txt");
    TRegexp e2("[A-Z][A-Z][A-Z]");

    TSystemFile* file;
    TString fname;
    TIter next(files);

    while ((file = (TSystemFile*)next())) {
      fname = file->GetName();
      if (fname.Contains(e)) {
        TString detName = fname(e2);
        detName.Resize(3);
        detectorsList.push_back(detName.Data());
      }
    }
  }

  // load geometry library
  gSystem->Load("libGeom");

  // create visualisation manager
  if (!TEveManager::Create()) {
    std::cout << "FATAL -- Could not create TEveManager!!" << std::endl;
    exit(0);
  }

  // load config file
  TEnv settings;
  o2::event_visualisation::ConfigurationManager::getInstance().getConfig(settings);

  // set OCDB path from config and set run number for which we want to generate geometry
  o2::ccdb::Manager::Instance()->setDefaultStorage(settings.GetValue("OCDB.default.path", "local:///home/maja/CERN/O2CDB"));
  o2::ccdb::Manager::Instance()->setRun(runNumber);

  // load geometry from OCDB
  o2::base::GeometryManager::loadGeometry();
  // load geometry from file
  //o2::base::GeometryManager::loadGeometry("O2geometry.root", "FAIRGeom");

  gGeoManager = o2::base::GeometryManager::getGeometry();
  gGeoManager->DefaultColors();

  // find main node for our detector
  TGeoNode* tnode = gGeoManager->GetTopNode();
  tnode->SetVisibility(kFALSE);

  std::string baseFileName(outputDir);
  baseFileName += "/simple_geom_";

  for (int i = 0; i < detectorsList.size(); i++) {
    TString path;
    TObjArray* list;
    Int_t depth;
    std::string line;

    const char* currentDetector;
    if (strcmp(detectorName, "ALL") == 0) {
      // if we update all detectors
      currentDetector = detectorsList[i].c_str();
    } else {
      // if we update specific detector, exit loop after one pass
      currentDetector = detectorName;
      i = detectorsList.size();
    }

    TEveGeoTopNode* eve_tnode = new TEveGeoTopNode(gGeoManager, tnode);
    eve_tnode->SetVisLevel(0);

    gEve->AddGlobalElement(eve_tnode);

    const char* detPath = Form("%s%s.txt", inPath.Data(), currentDetector);
    std::ifstream in(detPath, std::ios::in);
    std::cout << "Adding shapes from file:" << detPath << std::endl;

    int lineIter = 0;

    while (true) {
      std::getline(in, line, '\n');
      if (in.eof())
        break;

      if (line[0] == '#')
        continue;

      path = TString(line);

      if (!path.Contains("cave"))
        continue;

      list = path.Tokenize("/");
      depth = list->GetEntries();
      addNodes(tnode, eve_tnode, depth, depth, list);
      lineIter++;
    }
    in.close();

    if (lineIter == 0) {
      std::cout << "File for " << currentDetector << " is empty. Skipping..." << std::endl;
    } else {
      std::ostringstream str;
      str << baseFileName << currentDetector << ".root";
      eve_tnode->SaveExtract(str.str().c_str(), currentDetector, kTRUE);
    }
  }
}

void usage(char* name)
{
  std::cout << name << "[-d output directory] [-m detector name] [-r run number] [-h]" << std::endl;
}

Options* processCommandLine(int argc, char* argv[])
{
  static Options options;
  int opt;

  // put ':' at the beginning of the string so that program can distinguish between '?' and ':'
  while ((opt = getopt(argc, argv, ":d:m:r:h")) != -1) {
    switch (opt) {
      case 'd':
        options.outputDir = optarg;
        break;
      case 'm':
        options.detName = optarg;
        break;
      case 'r':
        options.runNumber = atoi(optarg);
        break;
      case 'h':
        usage(argv[0]);
        return nullptr;
      case ':':
        std::cout << "Option needs a value"
                  << "\n";
        usage(argv[0]);
        return nullptr;
      case '?':
        std::cout << "Unknown option: " << optarg << "\n";
        usage(argv[0]);
        return nullptr;
    }
  }

  // optind is for the extra arguments which are not parsed
  if (optind < argc) {
    std::cout << "Extra arguments:" << argv[optind] << "\n";
    usage(argv[0]);
    return nullptr;
  }

  return &options;
}

int main(int argc, char** argv)
{
  std::cout << "Launching simple geometry generator..." << std::endl;
  Options* options = processCommandLine(argc, argv);
  if (options == nullptr)
    exit(-1);

  // create ROOT application environment
  TApplication* app = new TApplication("o2-simple-geom-generate", &argc, argv);
  app->Connect("TEveBrowser", "CloseWindow()", "TApplication", app, "Terminate()");

  generateSimpleGeometry(options->detName.c_str(), options->runNumber, options->outputDir);

  TEveManager::Terminate();

  std::cout << "Simple geometries generated" << std::endl;

  return EXIT_SUCCESS;
}
