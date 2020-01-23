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
/// \brief   Creating O2CDB entry for EMCAL ideal alignment matrices
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#include "CCDB/CcdbApi.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsBase/GeometryManager.h"
#include "EventVisualisationBase/ConfigurationManager.h"

#include <TEnv.h>
#include <TString.h>
#include <TClonesArray.h>

#include <string>
#include <map>

using namespace o2::ccdb;
using namespace o2::event_visualisation;

// Based on AliEMCALSetAlignment.C, AliPHOSSetAlignment.C
int main(int argc, char** argv)
{
  CcdbApi api;
  std::map<std::string, std::string> metadata; // can be empty
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);
  const std::string ocdbStorage = settings.GetValue("OCDB.default.path", "local://$ALICE_ROOT/OCDB"); // default path to OCDB
  api.init(ocdbStorage);

  TClonesArray* array = new TClonesArray("o2::detectors::AlignParam", 20);
  TClonesArray& alobj = *array;

  o2::detectors::AlignParam a;

  // Null shifts and rotations
  double dx = 0., dy = 0., dz = 0., dpsi = 0., dtheta = 0., dphi = 0.;
  // Dummy volume identity
  int uid = o2::base::GeometryManager::getSensID(o2::detectors::DetID::EMC, 0);

  // EMCAL
  TString basePath = "EMCAL/FullSupermodule";
  for (int iModule = 0; iModule < 10; iModule++) {
    TString newPath = basePath;
    newPath += iModule + 1;

    new (alobj[iModule]) o2::detectors::AlignParam(newPath.Data(),
                                                   uid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  }

  //1/3 SMs
  new (alobj[10]) o2::detectors::AlignParam("EMCAL/OneThrdSupermodule1",
                                            uid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new (alobj[11]) o2::detectors::AlignParam("EMCAL/OneThrdSupermodule2",
                                            uid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);

  // DCal
  basePath = "EMCAL/DCALSupermodule";
  for (int iModule = 0; iModule < 6; iModule++) {
    TString newPath = basePath;
    newPath += iModule + 1;
    new (alobj[12 + iModule]) o2::detectors::AlignParam(newPath.Data(),
                                                        uid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  }

  //1/3 SMs
  new (alobj[18]) o2::detectors::AlignParam("EMCAL/OneThrdSupermodule3",
                                            uid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new (alobj[19]) o2::detectors::AlignParam("EMCAL/OneThrdSupermodule4",
                                            uid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);

  // Store arbitrary user object in strongly typed manner
  api.storeAsTFileAny(array, "EMCAL/Align/Data", metadata);

  array->Clear();

  // PHOS
  array = new TClonesArray("o2::detectors::AlignParam", 20);
  alobj = *array;

  // Dummy volume identity
  uid = o2::base::GeometryManager::getSensID(o2::detectors::DetID::PHS, 0);

  basePath = "PHOS/Module";
  const int nModules = 5;

  for (int iModule = 1; iModule <= nModules; iModule++) {
    TString newPath = basePath;
    newPath += iModule;
    new (alobj[iModule - 1]) o2::detectors::AlignParam(newPath.Data(),
                                                       uid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  }

  api.storeAsTFileAny(array, "PHOS/Align/Data", metadata);

  array->Delete();

  return EXIT_SUCCESS;
}
