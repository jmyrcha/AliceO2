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

#include <TSystem.h>

#include <iostream>

using namespace std;

namespace o2
{
namespace event_visualisation
{

ConfigurationManager& ConfigurationManager::getInstance()
{
  static ConfigurationManager instance;
  return instance;
}

void ConfigurationManager::getConfig(TEnv& settings, bool run2) const
{
  const std::string configName = run2 ? "eve_config" : "o2eve_config";
  const std::string methodName = "ConfigurationManager::getConfig";
  // TODO:
  // we need a way to point to the O2 installation directory
  //
  if (settings.ReadFile(Form(".%s", configName.c_str()), kEnvUser) < 0) {
    if (settings.ReadFile(Form("%s/.%s", gSystem->Getenv("HOME"), configName.c_str()), kEnvUser) < 0) {
      if (settings.ReadFile(Form("%s/%s", gSystem->Getenv("HOME"), configName.c_str()), kEnvUser) < 0) {
        std::cout << "WARNING: Could not find eve_config in home directory!\n";
        if (run2) {
          std::cout << "WARNING: Trying default one in $ALICE_ROOT/EVE/EveBase/\n";
          if (settings.ReadFile(Form("%s/EVE/EveBase/eve_config",
                                     gSystem->Getenv("ALICE_ROOT")),
                                kEnvUser) < 0) {
            std::cout << "ERROR: Could not find eve_config file!\n";
            exit(0);
          }
        } else {
          std::cout << "WARNING: Trying default one in $ALICEO2_INSTALL_PATH/EventVisualisation\n";
          if (settings.ReadFile(Form("%s/EventVisualisation/o2eve_config",
                                     gSystem->Getenv("ALICEO2_INSTALL_PATH")),
                                kEnvUser) < 0) {
            std::cout << "ERROR: Could not find o2eve_config file!\n";
            exit(0);
          }
        }
      } else {
        std::cout << "INFO: Using config at: " << gSystem->Getenv("HOME") << "/" << configName.c_str() << '\n';
      }
    } else {
      std::cout << "INFO: Using config at: " << gSystem->Getenv("HOME") << "/." << configName.c_str() << '\n';
    }
  } else {
    std::cout << "INFO: Using config at: " << configName.c_str() << '\n';
  }
}

} // namespace event_visualisation
} // namespace o2
