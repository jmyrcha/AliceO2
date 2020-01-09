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
/// \file    ConfigurationManager.cxx
/// \author  Jeremi Niedziela
/// \author  Julian Myrcha
///

#include "EventVisualisationBase/ConfigurationManager.h"

#include "FairLogger.h"

#include <TSystem.h>

namespace o2
{
namespace event_visualisation
{

ConfigurationManager& ConfigurationManager::getInstance()
{
  static ConfigurationManager instance;
  return instance;
}

void ConfigurationManager::getConfig(TEnv& settings) const
{
  TString fileName;
  if (settings.ReadFile(fileName = ".o2eve_config", kEnvUser) < 0) {
    //    LOG(WARNING) << "Could not find .o2eve_config in working directory! Trying .o2eve_config in home directory";
    if (settings.ReadFile(fileName = Form("%s/.o2eve_config", gSystem->Getenv("HOME")), kEnvUser) < 0) {
      //      LOG(WARNING) << "Could not find .o2eve_config in home directory! Trying o2eve_config in home directory";
      if (settings.ReadFile(fileName = Form("%s/o2eve_config", gSystem->Getenv("HOME")), kEnvUser) < 0) {
        //        LOG(WARNING) << "Could not find o2eve_config in home directory! Trying o2eve_config in O2/EventVisualisation";
        if (settings.ReadFile(fileName = Form("%s/EventVisualisation/o2eve_config",
                                              gSystem->Getenv("ALICEO2_INSTALL_PATH")),
                              kEnvUser) < 0) {
          LOG(FATAL) << "Could not find .o2eve_config or o2eve_config file!.";
          exit(0);
        }
      }
    }
  }
  // LOG(INFO) << Form("Using %s config settings", fileName.Data());
}

} // namespace event_visualisation
} // namespace o2
