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
  // TODO: we need a way to point to the O2 installation directory
  if (settings.ReadFile(Form("%s/.o2eve_config", gSystem->Getenv("HOME")), kEnvUser) < 0) {
    if (settings.ReadFile(Form("%s/o2eve_config", gSystem->Getenv("HOME")), kEnvUser) < 0) {
      LOG(WARNING) << "Could not find eve_config in home directory! Trying default one in O2/EventVisualisation/Base/";
      if (settings.ReadFile(Form("%s/EventVisualisation/o2eve_config", gSystem->Getenv("ALICEO2_INSTALL_PATH")), kEnvUser) < 0) {
        LOG(FATAL) << "Could not find eve_config file!";
      }
    }
  }
}

} // namespace event_visualisation
} // namespace o2
