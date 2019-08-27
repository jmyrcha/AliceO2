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
/// \file    FillO2CDBGeometry.cxx
/// \author  Maja Kabus
/// \brief   Creating O2CDB entry for O2geometry.root

#include "CCDB/Manager.h"
#include "CCDB/ConditionId.h"
#include "CCDB/Condition.h"
#include "CCDB/ConditionMetaData.h"
#include "CCDB/IdRunRange.h"

#include "DetectorsBase/GeometryManager.h"

#include <TGeoManager.h>
#include <TSystem.h>

#include <iostream>
#include <sstream>

using namespace o2::ccdb;

void processCommandLine(int argc, char* argv[], char** dir, int* run)
{
  int opt;
  // put ':' at the beginning of the string so that program can distinguish between '?' and ':'
  while ((opt = getopt(argc, argv, ":d:r:")) != -1) {
    switch (opt) {
      case 'd':
        *dir = optarg;
        break;
      case 'r':
        *run = atoi(optarg);
        break;
      case ':':
        std::cout << "Option needs a value"
                  << "\n";
        break;
      case '?':
        std::cout << "Unknown option: " << optarg << "\n";
        break;
    }
  }

  // optind is for the extra arguments which are not parsed
  if (optind < argc) {
    std::cout << "Extra arguments:" << argv[optind] << "\n";
  }
}

// Based on AliRoot GRP/UpdateCDBIdealGeom.C
// and on O2 CCDB/example/fill_local_ocdb.C
int main(int argc, char** argv)
{
  char* dir = "./";
  int run = 0;
  processCommandLine(argc, argv, &dir, &run);
  std::ostringstream str;
  str << dir << "O2geometry.root";
  std::string path = str.str();

  Manager* cdb = Manager::Instance();
  cdb->setDefaultStorage("local:///home/maja/CERN/O2CDB");

  cdb->setRun(run);
  ConditionId* id = new ConditionId("GRP/Geometry/Data", 0, IdRunRange::Infinity(), 0);
  ConditionMetaData* md = new ConditionMetaData();

  if (gSystem->AccessPathName(path.data())) {
    std::cout << "No geometry file found!" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Loading file: " << path << std::endl;
  if (TGeoManager::IsLocked())
    TGeoManager::UnlockGeometry();

  o2::base::GeometryManager::loadGeometry(path);
  gGeoManager->DefaultColors();

  cdb->putObject(gGeoManager, *id, md);

  // This is to allow macros launched after this one in the same session to find the
  // newly produced geometry.
  //cdb->queryStorages();

  return EXIT_SUCCESS;
}
