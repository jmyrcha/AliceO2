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
/// \brief   Creating O2CDB entry for O2geometry.root
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#include "CCDB/CcdbApi.h"
#include "EventVisualisationBase/ConfigurationManager.h"
#include "DetectorsBase/GeometryManager.h"

#include "FairLogger.h"

#include <TGeoManager.h>
#include <TSystem.h>

#include <sstream>

using namespace o2::ccdb;
using namespace o2::event_visualisation;

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
        LOG(FATAL) << "Option needs a value";
      case '?':
        LOG(FATAL) << "Unknown option: " << optarg;
    }
  }

  // optind is for the extra arguments which are not parsed
  if (optind < argc) {
    LOG(WARNING) << "Extra arguments:" << argv[optind];
  }
}

// Based on AliRoot GRP/UpdateCDBIdealGeom.C
int main(int argc, char** argv)
{
  char* dir = "./";
  int run = 0;
  processCommandLine(argc, argv, &dir, &run);
  std::ostringstream str;
  str << dir << "O2geometry.root";
  std::string path = str.str();

  CcdbApi api;
  std::map<std::string, std::string> metadata; // can be empty
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);
  const std::string ocdbStorage = settings.GetValue("OCDB.default.path", "local://$ALICE_ROOT/OCDB"); // default path to OCDB
  api.init(ocdbStorage);

  if (gSystem->AccessPathName(path.data())) {
    LOG(FATAL) << "No geometry file found!";
  }

  LOG(INFO) << "Loading file: " << path;
  if (TGeoManager::IsLocked())
    TGeoManager::UnlockGeometry();

  o2::base::GeometryManager::loadGeometry(path);
  gGeoManager->DefaultColors();

  api.storeAsTFileAny(gGeoManager, "GRP/Geometry/Data", metadata);

  return EXIT_SUCCESS;
}
