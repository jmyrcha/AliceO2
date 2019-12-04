// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#define BOOST_TEST_MODULE Test EventVisualisation Base
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "EventVisualisationBase/GeometryManager.h"
#include "EventVisualisationBase/VisualisationConstants.h"

#include "FairLogger.h"

#include <TApplication.h>
#include <TEveManager.h>

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
    }
    // TODO: Check if geometries work without GeometryManager.cxx:80, then no event manager needed
    //  Otherwise it crashes the tests
  }
  ~Fixture()
  {
    TEveManager::Terminate();
    gApplication->Terminate(0);
  }
};
BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_AUTO_TEST_CASE(drawO2GeometryTest)
{
  auto& geomManager = GeometryManager::getInstance();
  const std::string O2Geometries[] = {"EMC", "HMP", "ITS", "MCH", "MID", "PHS", "TOF", "TPC", "TRD"};
  for (auto detName : O2Geometries) {
    TEveGeoShape* shape = geomManager.getGeometryForDetector(detName);
    BOOST_CHECK_NE(shape, nullptr);
  }
}

BOOST_AUTO_TEST_CASE(dontFailOnAbsentGeometryTest)
{
  auto& geomManager = GeometryManager::getInstance();
  BOOST_CHECK_NO_THROW(geomManager.getGeometryForDetector("ACO"));
}
} // namespace event_visualisation
} // namespace o2