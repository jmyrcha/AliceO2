// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#define BOOST_TEST_MODULE Test EventVisualisation Detectors
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "EventVisualisationBase/DataReader.h"
#include "EventVisualisationDetectors/DataReaderAOD.h"
#include "EventVisualisationDetectors/DataReaderITS.h"
#include "EventVisualisationDetectors/DataReaderTPC.h"
#include "EventVisualisationBase/VisualisationConstants.h"
#include "EventVisualisationDataConverter/VisualisationEvent.h"
#include "EventVisualisationDataConverter/VisualisationCaloCell.h"
#include "EventVisualisationDataConverter/VisualisationTrack.h"
#include "EventVisualisationDataConverter/VisualisationCluster.h"
#include "EventVisualisationBase/DataInterpreter.h"
#include "EventVisualisationDetectors/DataInterpreterAOD.h"
#include "EventVisualisationDetectors/DataInterpreterITS.h"
#include "EventVisualisationDetectors/DataInterpreterTPC.h"

#include "FairLogger.h"

#include <TApplication.h>
#include <TEveManager.h>
#include <TString.h>
#include <TFile.h>

#include <memory>

namespace o2
{
namespace event_visualisation
{
/// What is done before and after all tests
// Needed for ITS - it uses geometry
struct Fixture {
  Fixture()
  {
    TApplication* app = new TApplication("o2eve", nullptr, nullptr);
    app->Connect("TEveBrowser", "CloseWindow()", "TApplication", app, "Terminate()");
    if (!TEveManager::Create(kTRUE, "FI")) {
      LOG(FATAL) << "Could not create TEveManager!";
    }
  }
  ~Fixture()
  {
    TEveManager::Terminate();
    gApplication->Terminate(0);
  }
};
BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_AUTO_TEST_CASE(Should_SetProperEventCount_When_DataFileIsCorrect_Test)
{
  LOG(STATE) << "Checking that each DataReader sets proper event count when data file is correct.";

  // Arrange
  DataReader* readers[] = {new DataReaderAOD(), new DataReaderITS(), new DataReaderTPC};
  int eventCounts[] = {150, 10, 10, 3865};

  // Act
  for (int i = 0; i < 4; i++) {
    readers[i]->open();
  }

  // Assert
  for (int i = 0; i < 4; i++) {
    BOOST_CHECK_EQUAL(readers[i]->getEventCount(), eventCounts[i]);
  }

  for (int i = 0; i < 4; i++) {
    delete readers[i];
  }
}

BOOST_AUTO_TEST_CASE(Should_ReturnEventData_When_EventNumberInRange_Test)
{
  LOG(STATE) << "Checking that each DataReader returns event data for the event specified.";

  // Arrange
  DataReader* readers[] = {new DataReaderAOD(), new DataReaderITS(), new DataReaderTPC};
  int eventCounts[] = {150, 10, 10, 3865};

  // Act
  for (int i = 0; i < 4; i++) {
    readers[i]->open();
  }

  // Assert
  for (int i = 0; i < 4; i++) {
    BOOST_CHECK_NE(readers[i]->getEventData(0), nullptr);
  }

  for (int i = 0; i < 4; i++) {
    delete readers[i];
  }
}

BOOST_AUTO_TEST_CASE(Should_ReturnNull_When_EventNumberOutsideRange_Test)
{
  LOG(STATE) << "Checking that each DataReader returns nullptr when event required is outside range.";

  // Arrange
  DataReader* readers[] = {new DataReaderAOD(), new DataReaderITS(), new DataReaderTPC};
  int eventCounts[] = {150, 10, 10, 3865};

  // Act
  for (int i = 0; i < 4; i++) {
    readers[i]->open();
  }

  // Assert
  for (int i = 0; i < 4; i++) {
    BOOST_CHECK_EQUAL(readers[i]->getEventData(eventCounts[i]), nullptr);
  }

  for (int i = 0; i < 4; i++) {
    delete readers[i];
  }
}

BOOST_AUTO_TEST_CASE(Should_InterpreEventTracks_When_CorrectEvent_Test)
{
  LOG(STATE) << "Checking that each DataInterpreter processes all tracks in an event.";

  // Arrange
  DataReader* readers[] = {new DataReaderAOD(), new DataReaderITS(), new DataReaderTPC};
  DataInterpreter* interpreters[] = {new DataInterpreterAOD(), new DataInterpreterITS(), new DataInterpreterTPC};
  int testedEvents[] = {140, 5, 5, 5};
  int trackCounts[] = {5946, 10, 33, 2};

  std::vector<std::unique_ptr<VisualisationEvent>> events(3);
  for (int i = 0; i < 3; i++) {
    events[i] = std::make_unique<VisualisationEvent>(0, 0, 0, 0, "", 0);
  }

  // Act
  for (int i = 0; i < 3; i++) {
    readers[i]->open();
    TObject* data = readers[i]->getEventData(testedEvents[i]);
    interpreters[i]->interpretDataForType(data, EVisualisationDataType::Tracks, *(events[i]));
  }

  // Assert
  for (int i = 0; i < 3; i++) {
    BOOST_CHECK_EQUAL(events[i]->getTrackCount(), trackCounts[i]);
  }

  for (int i = 0; i < 3; i++) {
    delete readers[i];
    delete interpreters[i];
  }
}

BOOST_AUTO_TEST_CASE(Should_InterpreEventClusters_When_CorrectEvent_Test)
{
  LOG(STATE) << "Checking that each DataInterpreter processes all clusters in an event.";

  // Arrange
  DataReader* readers[] = {new DataReaderAOD(), new DataReaderITS(), new DataReaderTPC};
  DataInterpreter* interpreters[] = {new DataInterpreterAOD(), new DataInterpreterITS(), new DataInterpreterTPC};
  int testedEvents[] = {140, 5, 5, 5};
  int clusterCounts[] = {0, 1393, 5472, 1};

  std::vector<std::unique_ptr<VisualisationEvent>> events(4);
  for (int i = 0; i < 3; i++) {
    events[i] = std::make_unique<VisualisationEvent>(0, 0, 0, 0, "", 0);
  }

  // Act
  for (int i = 0; i < 3; i++) {
    readers[i]->open();
    TObject* data = readers[i]->getEventData(testedEvents[i]);
    interpreters[i]->interpretDataForType(data, EVisualisationDataType::Clusters, *(events[i]));
  }

  // Assert
  for (int i = 0; i < 3; i++) {
    BOOST_CHECK_EQUAL(events[i]->getClusterCount(), clusterCounts[i]);
  }

  for (int i = 0; i < 3; i++) {
    delete readers[i];
    delete interpreters[i];
  }
}

BOOST_AUTO_TEST_CASE(Should_InterpretAODEventMuonTracks_When_CorrectAODEvent_Test)
{
  LOG(STATE) << "Checking that AOD DataInterpreter processes all muon tracks in an event.";

  // Arrange
  DataReader* reader = new DataReaderAOD();
  DataInterpreter* interpreter = new DataInterpreterAOD();
  std::unique_ptr<VisualisationEvent> event = std::make_unique<VisualisationEvent>(0, 0, 0, 0, "", 0);

  // Act
  reader->open();
  TObject* data = reader->getEventData(140);
  interpreter->interpretDataForType(data, EVisualisationDataType::Muon, *event);

  // Assert
  BOOST_CHECK_EQUAL(event->getMuonTrackCount(), 6);

  delete reader;
  delete interpreter;
}

BOOST_AUTO_TEST_CASE(Should_NotInterpretEventTracks_When_NotAODEvent_Test)
{
  LOG(STATE) << "Checking that each DataInterpreter except for AOD does not get muon tracks in an event.";

  // Arrange
  DataReader* readers[] = {new DataReaderITS(), new DataReaderTPC};
  DataInterpreter* interpreters[] = {new DataInterpreterITS(), new DataInterpreterTPC};

  std::vector<std::unique_ptr<VisualisationEvent>> events(4);
  for (int i = 0; i < 3; i++) {
    events[i] = std::make_unique<VisualisationEvent>(0, 0, 0, 0, "", 0);
  }

  // Act
  for (int i = 0; i < 2; i++) {
    readers[i]->open();
    TObject* data = readers[i]->getEventData(5);
    interpreters[i]->interpretDataForType(data, EVisualisationDataType::Muon, *(events[i]));
  }

  // Assert
  for (int i = 0; i < 2; i++) {
    BOOST_CHECK_EQUAL(events[i]->getMuonTrackCount(), 0);
  }

  for (int i = 0; i < 2; i++) {
    delete readers[i];
    delete interpreters[i];
  }
}

BOOST_AUTO_TEST_CASE(Should_InterpretAODEventCaloCells_When_CorrectAODEvent_Test)
{
  LOG(STATE) << "Checking that AOD DataInterpreter processes all calorimeter cells in an event.";

  // Arrange
  DataReader* reader = new DataReaderAOD();
  DataInterpreter* interpreter = new DataInterpreterAOD();
  std::unique_ptr<VisualisationEvent> event = std::make_unique<VisualisationEvent>(0, 0, 0, 0, "", 0);

  // Act
  reader->open();
  TObject* data = reader->getEventData(140);
  interpreter->interpretDataForType(data, EVisualisationDataType::Calo, *event);

  // Assert
  BOOST_CHECK_EQUAL(event->getCaloCellsCount(), 558);

  delete reader;
  delete interpreter;
}

BOOST_AUTO_TEST_CASE(Should_NotInterpretEventCaloCells_When_NotAODEvent_Test)
{
  LOG(STATE) << "Checking that each DataInterpreter except for AOD does not get calorimeter cells in an event.";

  // Arrange
  DataReader* readers[] = {new DataReaderITS(), new DataReaderTPC};
  DataInterpreter* interpreters[] = {new DataInterpreterITS(), new DataInterpreterTPC};

  std::vector<std::unique_ptr<VisualisationEvent>> events(4);
  for (int i = 0; i < 3; i++) {
    events[i] = std::make_unique<VisualisationEvent>(0, 0, 0, 0, "", 0);
  }

  // Act
  for (int i = 0; i < 2; i++) {
    readers[i]->open();
    TObject* data = readers[i]->getEventData(5);
    interpreters[i]->interpretDataForType(data, EVisualisationDataType::Calo, *(events[i]));
  }

  // Assert
  for (int i = 0; i < 2; i++) {
    BOOST_CHECK_EQUAL(events[i]->getCaloCellsCount(), 0);
  }

  for (int i = 0; i < 2; i++) {
    delete readers[i];
    delete interpreters[i];
  }
}

} // namespace event_visualisation
} // namespace o2
