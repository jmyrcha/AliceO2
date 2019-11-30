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
/// \file    Initializer.cxx
/// \author  Jeremi Niedziela
/// \author  julian.myrcha@cern.ch
/// \author  p.nowakowski@cern.ch
///

#include "EventVisualisationView/Initializer.h"

#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationView/EventManager.h"
#include "EventVisualisationBase/GeometryManager.h"
#include "EventVisualisationView/MultiView.h"
#include "EventVisualisationBase/VisualisationConstants.h"
#include "EventVisualisationView/EventManagerFrame.h"
#include "EventVisualisationBase/DataSourceOffline.h"
#include "EventVisualisationDetectors/DataReaderAOD.h"
#include "EventVisualisationDetectors/DataReaderVSD.h"
#include "EventVisualisationDetectors/DataReaderITS.h"
#include "EventVisualisationDetectors/DataReaderTPC.h"
#include "EventVisualisationDetectors/DataInterpreterAOD.h"
#include "EventVisualisationDetectors/DataInterpreterVSD.h"
#include "EventVisualisationDetectors/DataInterpreterITS.h"
#include "EventVisualisationDetectors/DataInterpreterTPC.h"

#include <TGTab.h>
#include <TEnv.h>
#include <TEveBrowser.h>
#include <TEveManager.h>
#include <TRegexp.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TEveWindowManager.h>
#include <iostream>
#include <TFile.h>
#include <TGeoManager.h>

using namespace std;

namespace o2
{
namespace event_visualisation
{

void Initializer::setup(const Options options, EventManager::EDataSource defaultDataSource)
{
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);

  const bool fullscreen = settings.GetValue("fullscreen.mode", false);                           // hide left and bottom tabs
  const string ocdbStorage = settings.GetValue("OCDB.default.path", "local://$ALICE_ROOT/OCDB"); // default path to OCDB
  cout << "Initializer -- OCDB path:" << ocdbStorage << endl;

  auto& eventManager = EventManager::getInstance();
  eventManager.setDataSourceType(defaultDataSource);
  eventManager.setCdbPath(ocdbStorage);

  if (options.vsd)
    eventManager.registerDetector(new DataReaderVSD(), new DataInterpreterVSD(), EVisualisationGroup::VSD);
  if (options.its)
    eventManager.registerDetector(new DataReaderITS(), new DataInterpreterITS(), EVisualisationGroup::ITS);
  if (options.tpc)
    eventManager.registerDetector(new DataReaderTPC(), new DataInterpreterTPC(), EVisualisationGroup::TPC);
  if (options.aod)
    eventManager.registerDetector(new DataReaderAOD(), new DataInterpreterAOD(), EVisualisationGroup::AOD);

  eventManager.Open();

  // Setup windows size, fullscreen and focus
  TEveBrowser* browser = gEve->GetBrowser();
  browser->GetTabRight()->SetTab(1);
  browser->MoveResize(0, 0, gClient->GetDisplayWidth(), gClient->GetDisplayHeight() - 32);

  browser->StartEmbedding(TRootBrowser::kBottom);
  EventManagerFrame* frame = new EventManagerFrame(eventManager);
  browser->StopEmbedding("EventCtrl");

  if (fullscreen) {
    ((TGWindow*)gEve->GetBrowser()->GetTabLeft()->GetParent())->Resize(1, 0);
    ((TGWindow*)gEve->GetBrowser()->GetTabBottom()->GetParent())->Resize(0, 1);
  }
  gEve->GetBrowser()->Layout();

  gEve->AddEvent(&EventManager::getInstance());

  frame->refresh(true);
}

} // namespace event_visualisation
} // namespace o2
