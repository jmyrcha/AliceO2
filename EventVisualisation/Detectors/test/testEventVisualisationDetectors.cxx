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
#include "EventVisualisationDetectors/DataReaderVSD.h"
#include "EventVisualisationBase/VisualisationConstants.h"
#include "EventVisualisationDataConverter/VisualisationEvent.h"
#include "EventVisualisationDataConverter/VisualisationCaloCell.h"
#include "EventVisualisationDataConverter/VisualisationTrack.h"
#include "EventVisualisationDataConverter/VisualisationCluster.h"
#include "EventVisualisationBase/DataInterpreter.h"
#include "EventVisualisationDetectors/DataInterpreterAOD.h"
#include "EventVisualisationDetectors/DataInterpreterITS.h"
#include "EventVisualisationDetectors/DataInterpreterTPC.h"
#include "EventVisualisationDetectors/DataInterpreterVSD.h"

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
    LOG(INFO) << "Initializing TEveManager";
    if (!TEveManager::Create(kTRUE, "FI")) {
      LOG(FATAL) << "Could not create TEveManager!";
      exit(1);
    }
  }
  ~Fixture()
  {
    TEveManager::Terminate();
    gApplication->Terminate(0);
  }
};
BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_AUTO_TEST_CASE(eventCountTest)
{
  DataReader* readers[] = { new DataReaderAOD(), new DataReaderITS(), new DataReaderTPC, new DataReaderVSD };
  int eventCounts[] = {150, 10, 10, 3865};
  for (int i = 0; i < 4; i++) {
    readers[i]->open();
    BOOST_CHECK_EQUAL(readers[i]->GetEventCount(), eventCounts[i]);
    delete readers[i];
  }
}

BOOST_AUTO_TEST_CASE(getEventDataTest)
{
  DataReader* readers[] = { new DataReaderAOD(), new DataReaderITS(), new DataReaderTPC, new DataReaderVSD };
  int eventCounts[] = {150, 10, 10, 3865};
  for (int i = 0; i < 4; i++) {
    readers[i]->open();
    BOOST_CHECK_NE(readers[i]->getEventData(0), nullptr);
    BOOST_CHECK_EQUAL(readers[i]->getEventData(eventCounts[i]), nullptr);
    delete readers[i];
  }
}

BOOST_AUTO_TEST_CASE(interpretDataForTypeTest)
{
  // FIXME: VSD breaks, where is any correct code??
  DataReader* readers[] = { new DataReaderAOD(), new DataReaderITS(), new DataReaderTPC, new DataReaderVSD };
  DataInterpreter* interpreters[] = { new DataInterpreterAOD(), new DataInterpreterITS(), new DataInterpreterTPC, new DataInterpreterVSD };
  int testedEvents[] = { 140, 5, 5, 5 };

  int trackCounts[] = { 5946, 10, 33, 2 };
  int clusterCounts[] = { 0, 1393, 5472, 1 };
  int muonCounts[] = { 6, 0, 0, 0 };
  int caloCounts[] = { 558, 0, 0, 0 };

  for (int i = 0; i < 3; i++) {
    readers[i]->open();
    TObject* data = readers[i]->getEventData(testedEvents[i]);

    std::unique_ptr<VisualisationEvent> event = std::make_unique<VisualisationEvent>(0, 0, 0, 0, "", 0);
    for (int dataType = 0; dataType < EVisualisationDataType::NdataTypes; ++dataType) {
      interpreters[i]->interpretDataForType(data, (EVisualisationDataType)dataType, *event);
    }

    BOOST_CHECK_EQUAL(event->getTrackCount(), trackCounts[i]);
    BOOST_CHECK_EQUAL(event->getClusterCount(), clusterCounts[i]);
    BOOST_CHECK_EQUAL(event->getMuonTrackCount(), muonCounts[i]);
    BOOST_CHECK_EQUAL(event->getCaloCellsCount(), caloCounts[i]);

    delete readers[i];
    delete interpreters[i];
  }
}


} // namespace event_visualisation
} // namespace o2