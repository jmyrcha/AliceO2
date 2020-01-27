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
/// \file    DataInterpreterAOD.h
/// \brief   Converting AOD data to Event Visualisation primitives
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#ifndef ALICE_O2_EVENTVISUALISATION_DETECTORS_DATAINTERPRETERAOD_H
#define ALICE_O2_EVENTVISUALISATION_DETECTORS_DATAINTERPRETERAOD_H

/// This class overrides DataInterpreter and implements method
/// returning visualisation objects representing data from AOD file.

#include "EventVisualisationBase/DataInterpreter.h"
#include "EventVisualisationBase/VisualisationConstants.h"
#include "EventVisualisationDataConverter/VisualisationCaloCell.h"
#include "EMCALBase/Geometry.h"
#include "PHOSBase/Geometry.h"

#include <TTree.h>

namespace o2
{
namespace event_visualisation
{
class DataInterpreterAOD : public DataInterpreter
{
 public:
  // Default constructor
  DataInterpreterAOD();

  // Default destructor
  ~DataInterpreterAOD() final = default;

  void interpretDataForType(TObject* data, EVisualisationDataType type, VisualisationEvent& event) final;

 private:
  void interpretTracks(TFile* AODFile, int eventId, VisualisationEvent& event) final;

  void interpretAODCaloCells(TFile* AODFile, int eventId, VisualisationEvent& event);
  void interpretEMCALCell(int absID, float amplitude, VisualisationEvent& event);
  void interpretPHOSCell(int absID, float amplitude, VisualisationEvent& event);

  void interpretMuonTracks(TFile* AODFile, int eventId, VisualisationEvent& event);

  bool cutCell(std::vector<float> caloAmplitudes, std::vector<int> caloAbsIds, int caloCount,
               float amplitude, char caloType, float maxCellEnergy, int maxCellEnergyAbsId);
  std::tuple<int, int, int> getModuleNumberColAndRow(int absId, char caloType);
  float getECross(std::vector<float>& caloAmplitudes, std::vector<int>& caloAbsIds,
                  char caloType, int imod, int icol, int irow);
  float getCaloCellAmplitude(std::vector<float> caloAmplitudes, std::vector<int> caloAbsIds, int absId);

  // Calorimeter geometries
  o2::emcal::Geometry* mEMCALGeom = nullptr;
  o2::phos::Geometry* mPHOSGeom = nullptr;

  // From AliRoot emcal_esdclustercells.C macro
  // Cluster cuts          { PHOS, EMCAL }
  int mNumMinCellsCut[2] = {3, 2};    /// Number of cells in cluster must be larger than this value.
  int mNumMaxCellsCut[2] = {60, 30};  /// Number of cells in cluster must be smaller than this value.
  float mEnergyCut[2] = {0.30, 0.30}; /// Cluster energy must be larger than this value.
  float mM02LowCut[2] = {0.70, 0.10}; /// Cluster shower shape major axis must be larger than this value.
  float mM02HigCut[2] = {7.00, 7.00}; /// Cluster shower shape major axis must be smaller than this value.
  float mM20LowCut[2] = {0.50, -1.0}; /// Cluster shower shape lower axis must be larger than this value.
  float mExoCut = 0.95;               /// Reject clusters with this exoticity value.
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_DETECTORS_DATAINTERPRETERAOD_H
