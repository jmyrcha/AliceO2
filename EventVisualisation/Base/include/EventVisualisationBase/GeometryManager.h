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
/// \file    GeometryManager.h
/// \author  Jeremi Niedziela
/// \author  julian.myrcha@cern.ch
///

#ifndef ALICE_O2_EVENTVISUALISATION_BASE_GEOMETRYMANAGER_H
#define ALICE_O2_EVENTVISUALISATION_BASE_GEOMETRYMANAGER_H

#include <TEveGeoShape.h>

#include <string>

namespace o2
{
namespace event_visualisation
{

/// GeometryManager allows access to geometries of detectors.
///
/// GeometryManager is a singleton class which opens ROOT files with
/// simplified geometries, reads drawing parameters (such as color or transparency)
/// from the config file, prepares and return ready-to-draw volumes.

class GeometryManager
{
 public:
  /// Returns an instance of GeometryManager
  static GeometryManager& getInstance();

  /// Returns ROOT shapes describing simplified geometry of given detector
  TEveGeoShape* getGeometryForDetector(std::string detectorName);

  /// Sets geometry to R2/R3
  void setR2Geometry(bool R2Geometry) { this->mR2Geometry = R2Geometry; }
  bool getR2Geometry() const { return this->mR2Geometry; }

 private:
  static GeometryManager* sInstance;

  /// If using R2 geometry
  bool mR2Geometry;

  /// Goes through all children nodes of geometry shape and sets drawing options
  void drawDeep(TEveGeoShape* geomShape, Color_t color, char transparency, Color_t lineColor);

  /// Default constructor
  GeometryManager();
  /// Default destructor
  ~GeometryManager() = default;
  /// Deleted copy constructor
  GeometryManager(GeometryManager const&) = delete;
  /// Deleted assignment operator
  void operator=(GeometryManager const&) = delete;
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_BASE_GEOMETRYMANAGER_H
