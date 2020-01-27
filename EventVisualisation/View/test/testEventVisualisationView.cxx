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
/// \file    testEventVisualisationView.cxx
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#define BOOST_TEST_MODULE Test EventVisualisation View
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "EventVisualisationView/EventManager.h"
#include "EventVisualisationView/MultiView.h"
#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationBase/DataSource.h"
#include "EventVisualisationBase/VisualisationConstants.h"
#include "EventVisualisationDetectors/DataReaderAOD.h"
#include "EventVisualisationDetectors/DataInterpreterAOD.h"
#include "EventVisualisationDetectors/DataReaderTPC.h"
#include "EventVisualisationDetectors/DataInterpreterTPC.h"
#include "EventVisualisationDetectors/DataReaderITS.h"
#include "EventVisualisationDetectors/DataInterpreterITS.h"

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
    if (!TEveManager::Create(kTRUE, "FI")) {
      LOG(FATAL) << "Could not create TEveManager!";
    }
    auto& eventManager = EventManager::getInstance();
    eventManager.setDataSourceType(EDataSource::SourceOffline);
    gEve->AddEvent(&eventManager);

    DataReader* readers[] = {new DataReaderAOD(), new DataReaderITS(), new DataReaderTPC};
    DataInterpreter* interpreters[] = {new DataInterpreterAOD(), new DataInterpreterITS(), new DataInterpreterTPC};
    EVisualisationGroup visualisationGroups[] = {EVisualisationGroup::AOD, EVisualisationGroup::ITS, EVisualisationGroup::TPC};
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

BOOST_AUTO_TEST_CASE(Should_SetHighestEventCount_When_MultipleDataFilesRead_Test)
{
  LOG(STATE) << "Checking that event count is set to biggest event from all input files.";

  // Arrange, act - also in tests fixture constructor
  auto& eventManager = EventManager::getInstance();

  // Assert
  BOOST_CHECK_EQUAL(eventManager.getDataSource()->getEventCount(), 150);
}

BOOST_AUTO_TEST_CASE(Should_DisplayTracks_When_CorrectEvent_Test)
{
  LOG(STATE) << "Checking that tracks are displayed in eve.";

  // Arrange
  auto& eventManager = EventManager::getInstance();

  // Act
  // Event with tracks
  eventManager.GotoEvent(140);
  TEveElementList* tracks = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Tracks"));

  // Assert
  BOOST_CHECK_NE(tracks, nullptr);
}

BOOST_AUTO_TEST_CASE(Should_DisplayAODTracksAccordingToSettings_Test)
{
  LOG(STATE) << "Checking that AOD tracks are displayed according to settings.";

  // Arrange
  auto& eventManager = EventManager::getInstance();

  // Act
  // Event with AOD data only
  eventManager.GotoEvent(140);
  TEveElementList* tracks = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Tracks"));

  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);

  TEveElementList* AODTracks;
  if (settings.GetValue("tracks.byPt.show", false)) {
    AODTracks = (TEveElementList*)(tracks->FindChild("AOD tracks by Pt"));
  } else if (settings.GetValue("tracks.byType.show", false)) {
    AODTracks = (TEveElementList*)(tracks->FindChild("AOD tracks by type"));
  }

  // Assert
  BOOST_CHECK_NE(AODTracks, nullptr);
}

BOOST_AUTO_TEST_CASE(Should_DisplayITSTracksAccordingToSettings_Test)
{
  LOG(STATE) << "Checking that ITS tracks are displayed according to settings.";

  // Arrange
  auto& eventManager = EventManager::getInstance();

  // Act
  // Event with ITS tracks only
  eventManager.GotoEvent(17);
  TEveElementList* tracks = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Tracks"));

  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);

  TEveElementList* ITSTracks;
  if (settings.GetValue("tracks.byPt.show", false)) {
    ITSTracks = (TEveElementList*)(tracks->FindChild("ITS tracks by Pt"));
  } else {
    ITSTracks = (TEveElementList*)(tracks->FindChild("ITS"));
  }

  // Assert
  BOOST_CHECK_NE(ITSTracks, nullptr);
}

BOOST_AUTO_TEST_CASE(Should_DisplayTPCTracksAccordingToSettings_Test)
{
  LOG(STATE) << "Checking that TPC tracks are displayed according to settings.";

  // Arrange
  auto& eventManager = EventManager::getInstance();

  // Act
  // Event with TPC tracks
  eventManager.GotoEvent(7);
  TEveElementList* tracks = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Tracks"));

  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);

  TEveElementList* TPCTracks;
  if (settings.GetValue("tracks.byPt.show", false)) {
    TPCTracks = (TEveElementList*)(tracks->FindChild("TPC tracks by Pt"));
  } else {
    TPCTracks = (TEveElementList*)(tracks->FindChild("TPC"));
  }

  // Assert
  BOOST_CHECK_NE(TPCTracks, nullptr);
}

BOOST_AUTO_TEST_CASE(Should_DisplayClusters_When_CorrectEvent_Test)
{
  LOG(STATE) << "Checking that clusters are displayed in eve.";

  // Arrange
  auto& eventManager = EventManager::getInstance();

  // Act
  // Event with clusters
  eventManager.GotoEvent(1);
  TEveElementList* clusters = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Clusters"));

  // Assert
  BOOST_CHECK_NE(clusters, nullptr);
}

BOOST_AUTO_TEST_CASE(Should_DisplayITSClusters_When_CorrectEvent_Test)
{
  LOG(STATE) << "Checking that ITS clusters are displayed in eve.";

  // Arrange
  auto& eventManager = EventManager::getInstance();

  // Act
  // Event with ITS clusters
  eventManager.GotoEvent(1);
  TEveElementList* clusters = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Clusters"));
  TEveElementList* ITSClusters = (TEveElementList*)(clusters->FindChild("ITS"));

  // Assert
  BOOST_CHECK_NE(ITSClusters, nullptr);
}

BOOST_AUTO_TEST_CASE(Should_DisplayTPCClusters_When_CorrectEvent_Test)
{
  LOG(STATE) << "Checking that TPC clusters are displayed in eve.";

  // Arrange
  auto& eventManager = EventManager::getInstance();

  // Act
  // Event with TPC clusters
  eventManager.GotoEvent(1);
  TEveElementList* clusters = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Clusters"));
  TEveElementList* TPCClusters = (TEveElementList*)(clusters->FindChild("TPC"));

  // Assert
  BOOST_CHECK_NE(TPCClusters, nullptr);
}

BOOST_AUTO_TEST_CASE(Should_DisplayCaloCells_When_CorrectAODEvent_Test)
{
  LOG(STATE) << "Checking that calorimeter cells are displayed in eve for an AOD event.";

  // Arrange
  auto& eventManager = EventManager::getInstance();

  // Act
  // Event with AOD data only
  eventManager.GotoEvent(140);
  TEveElementList* calo = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Calo"));

  // Assert
  BOOST_CHECK_NE(calo, nullptr);
}

BOOST_AUTO_TEST_CASE(Should_DisplayCalo3DTowers_When_CorrectAODEvent_Test)
{
  LOG(STATE) << "Checking that calorimeter towers are displayed in eve for an AOD event.";

  // Arrange
  auto& eventManager = EventManager::getInstance();

  // Act
  // Event with AOD data only
  eventManager.GotoEvent(140);
  TEveElementList* calo = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Calo"));
  TEveElementList* histList = (TEveElementList*)(calo->FindChild("3D Histogram"));
  BOOST_CHECK_NE(histList, nullptr);
  TEveCalo3D* calo3d = (TEveCalo3D*)histList->FirstChild();
  BOOST_CHECK_NE(calo3d, nullptr);
  TEveCaloDataHist* data = (TEveCaloDataHist*)calo3d->GetData();
  BOOST_CHECK_NE(data, nullptr);
}

BOOST_AUTO_TEST_CASE(Should_DisplayEMCALCaloCells_When_CorrectAODEvent_Test)
{
  LOG(STATE) << "Checking that EMCAL calorimeter cells are displayed in eve for an AOD event.";

  // Arrange
  auto& eventManager = EventManager::getInstance();

  // Act
  // Event with AOD data only
  eventManager.GotoEvent(140);
  TEveElementList* calo = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Calo"));
  TEveElementList* EMCALList = (TEveElementList*)(calo->FindChild("EMCAL"));

  // Assert
  BOOST_CHECK_NE(EMCALList, nullptr);
}

BOOST_AUTO_TEST_CASE(Should_DisplayPHOSCaloCells_When_CorrectAODEvent_Test)
{
  LOG(STATE) << "Checking that PHOS calorimeter cells are displayed in eve for an AOD event.";

  // Arrange
  auto& eventManager = EventManager::getInstance();

  // Act
  // Event with AOD data only
  eventManager.GotoEvent(140);
  TEveElementList* calo = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Calo"));
  TEveElementList* PHOSList = (TEveElementList*)(calo->FindChild("PHOS"));

  // Assert
  BOOST_CHECK_NE(PHOSList, nullptr);
}

BOOST_AUTO_TEST_CASE(Should_DisplayMuonTracks_When_CorrectAODEvent_Test)
{
  LOG(STATE) << "Checking that muon tracks are displayed in eve for an AOD event.";

  // Arrange
  auto& eventManager = EventManager::getInstance();

  // Act
  // Event with AOD data only
  eventManager.GotoEvent(140);
  TEveElementList* muonList = (TEveElementList*)(gEve->GetCurrentEvent()->FindChild("Muon"));

  // Assert
  BOOST_CHECK_NE(muonList, nullptr);
}

} // namespace event_visualisation
} // namespace o2
