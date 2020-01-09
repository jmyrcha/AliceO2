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
/// \file    testEventVisualisationGeometryManager.cxx
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#define BOOST_TEST_MODULE Test EventVisualisation Base
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "EventVisualisationBase/GeometryManager.h"
#include "EventVisualisationBase/VisualisationConstants.h"

#include "FairLogger.h"

#include <TApplication.h>
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveViewer.h>

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
    gEve->AddEvent(new TEveEventManager("Event", ""));
    TEveViewer* viewer = gEve->SpawnNewViewer("3D View", "");
    viewer->AddScene(gEve->GetEventScene());
  }
  ~Fixture()
  {
    TEveManager::Terminate();
    gApplication->Terminate(0);
  }
};
BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_AUTO_TEST_CASE(Should_LoadDetectorGeometry_When_GeometryIsO2_Test)
{
  LOG(STATE) << "Checking that each O2 geometry is loaded properly.";

  // Arrange
  auto& geomManager = GeometryManager::getInstance();
  const std::string O2Geometries[] = {"EMC", "HMP", "ITS", "MCH", "MID", "PHS", "TOF", "TPC", "TRD"};
  int geometriesCount = 9;
  TEveGeoShape* shapes[geometriesCount];
  for (int i = 0; i < geometriesCount; i++) {
    shapes[i] = nullptr;
  }

  // Act
  for (int i = 0; i < geometriesCount; i++) {
    shapes[i] = geomManager.getGeometryForDetector(O2Geometries[i]);
  }

  // Assert
  for (int i = 0; i < geometriesCount; i++) {
    BOOST_CHECK_NE(shapes[i], nullptr);
  }
}

BOOST_AUTO_TEST_CASE(Should_NotLoad_And_NotTerminate_When_GeometryDoesNotExist_Test)
{
  LOG(STATE) << "Checking that wrong geometry specifier raises a non-fatal error.";

  // Arrange
  auto& geomManager = GeometryManager::getInstance();

  // Act, assert
  BOOST_CHECK_NO_THROW(geomManager.getGeometryForDetector("ACO"));
}
} // namespace event_visualisation
} // namespace o2
