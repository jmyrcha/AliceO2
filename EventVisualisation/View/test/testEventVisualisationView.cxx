// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#define BOOST_TEST_MODULE Test EventVisualisation View
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "EventVisualisationView/EventManager.h"
#include "EventVisualisationView/MultiView.h"
#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationBase/DataSource.h"
#include "EventVisualisationDetectors/DataReaderAOD.h"
#include "EventVisualisationDetectors/DataInterpreterAOD.h"
#include "EventVisualisationDetectors/DataReaderTPC.h"
#include "EventVisualisationDetectors/DataInterpreterTPC.h"
#include "EventVisualisationDetectors/DataReaderITS.h"
#include "EventVisualisationDetectors/DataInterpreterITS.h"
#include "EventVisualisationDetectors/DataReaderVSD.h"
#include "EventVisualisationDetectors/DataInterpreterVSD.h"

#include "FairLogger.h"

#include <TApplication.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEveCalo.h>
#include <TEveCaloData.h>

namespace o2
{
namespace event_visualisation
{
/// What is done before and after all tests
struct Fixture {
  Fixture()
  {
    TApplication* app = new TApplication("o2eve", nullptr, nullptr);
    app->Connect("TEveBrowser", "CloseWindow()", "TApplication", app, "Terminate()");
    LOG(INFO) << "Initializing TEveManager";
    if (!TEveManager::Create(kTRUE, "FI")) {
      LOG(FATAL) << "Could not create TEveManager!";
      exit(1);
    }
    auto& eventManager = EventManager::getInstance();
    eventManager.setDataSourceType(EventManager::EDataSource::SourceOffline);
    gEve->AddEvent(&eventManager);
    // FIXME: VSD breaks, where is any correct code??
    DataReader* readers[] = { new DataReaderAOD(), new DataReaderITS(), new DataReaderTPC, new DataReaderVSD };
    DataInterpreter* interpreters[] = { new DataInterpreterAOD(), new DataInterpreterITS(), new DataInterpreterTPC, new DataInterpreterVSD };
    EVisualisationGroup visualisationGroups[] = { EVisualisationGroup::AOD, EVisualisationGroup::ITS, EVisualisationGroup::TPC, EVisualisationGroup::VSD };
    for (int i = 0; i < 3; i++) {
      eventManager.registerDetector(readers[i], interpreters[i], visualisationGroups[i]);
    }
    eventManager.Open();
  }
  ~Fixture()
  {
    TEveManager::Terminate();
    gApplication->Terminate(0);
  }
};
BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_AUTO_TEST_CASE(OpenTest)
{
  auto& eventManager = EventManager::getInstance();
  BOOST_CHECK_EQUAL(eventManager.getDataSource()->GetEventCount(), 150);
}

BOOST_AUTO_TEST_CASE(GoToEventAODTest)
{
  auto& eventManager = EventManager::getInstance();

  // Event with all kind of AOD data
  eventManager.GotoEvent(140);

  TEveElementList* tracks = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Tracks"));
  BOOST_CHECK_NE(tracks, nullptr);

  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);

  TEveElementList* AODTracks;
  if (settings.GetValue("tracks.byPt.show", false)) {
    AODTracks = (TEveElementList*)(tracks->FindChild("AOD tracks by Pt"));
  } else if (settings.GetValue("tracks.byType.show", false)) {
    AODTracks = (TEveElementList*)(tracks->FindChild("AOD tracks by type"));
  }
  BOOST_CHECK_NE(AODTracks, nullptr);

  TEveElementList* calo = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Calo"));
  BOOST_CHECK_NE(calo, nullptr);

  TEveElementList* histList = (TEveElementList*)(calo->FindChild("3D Histogram"));
  BOOST_CHECK_NE(histList, nullptr);
  TEveCalo3D* calo3d = (TEveCalo3D*)histList->FirstChild();
  BOOST_CHECK_NE(calo3d, nullptr);
  TEveCaloDataHist* data = (TEveCaloDataHist*)calo3d->GetData();
  BOOST_CHECK_NE(data, nullptr);

  TEveElementList* EMCALList = (TEveElementList*)(calo->FindChild("EMCAL"));
  BOOST_CHECK_NE(EMCALList, nullptr);
  TEveElementList* PHOSList = (TEveElementList*)(calo->FindChild("PHOS"));
  BOOST_CHECK_NE(PHOSList, nullptr);

  TEveElementList* muonList = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Muon"));
  BOOST_CHECK_NE(muonList, nullptr);
  TEveTrackList* muonTracks = (TEveTrackList*)(muonList->FindChild("Ghost"));
  BOOST_CHECK_NE(muonTracks, nullptr);
}

BOOST_AUTO_TEST_CASE(GoToEventITSTest)
{
  auto& eventManager = EventManager::getInstance();

  // Event with all kind of ITS data
  eventManager.GotoEvent(1);

  TEveElementList* tracks = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Tracks"));
  BOOST_CHECK_NE(tracks, nullptr);
  TEveElementList* ITSTracks = (TEveElementList*)(tracks->FindChild("ITS tracks by type"));
  BOOST_CHECK_NE(ITSTracks, nullptr);

  TEveElementList* clusters = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Clusters"));
  BOOST_CHECK_NE(clusters, nullptr);
  TEveElementList* ITSClusters = (TEveElementList*)(clusters->FindChild("ITS"));
  BOOST_CHECK_NE(ITSClusters, nullptr);
}

BOOST_AUTO_TEST_CASE(GoToEventTPCTest)
{
  auto& eventManager = EventManager::getInstance();

  // Event with all kind of TPC data
  eventManager.GotoEvent(7);

  TEveElementList* tracks = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Tracks"));
  BOOST_CHECK_NE(tracks, nullptr);
  TEveElementList* TPCTracks = (TEveElementList*)(tracks->FindChild("TPC tracks by type"));
  BOOST_CHECK_NE(TPCTracks, nullptr);

  TEveElementList* clusters = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Clusters"));
  BOOST_CHECK_NE(clusters, nullptr);
  TEveElementList* TPCClusters = (TEveElementList*)(clusters->FindChild("TPC"));
  BOOST_CHECK_NE(TPCClusters, nullptr);
}

} // namespace event_visualisation
} // namespace o2
