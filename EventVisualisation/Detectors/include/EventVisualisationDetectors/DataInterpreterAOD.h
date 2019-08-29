// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DataInterpreterAOD.h
/// \brief converting AOD data to Event Visualisation primitives
/// \author Maja Kabus

#ifndef ALICE_O2_EVENTVISUALISATION_DETECTORS_DATAINTERPRETERAOD_H
#define ALICE_O2_EVENTVISUALISATION_DETECTORS_DATAINTERPRETERAOD_H

/// This class overrides DataInterpreter and implements method
/// returning visualisation objects representing data from AOD file.

#include "EventVisualisationBase/DataInterpreter.h"
#include "EventVisualisationBase/VisualisationConstants.h"
#include "EMCALBase/Geometry.h"
#include "PHOSBase/Geometry.h"

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
  ~DataInterpreterAOD() final;

  std::unique_ptr<VisualisationEvent> interpretDataForType(TObject* data, EVisualisationDataType type) final;

 private:
  std::unique_ptr<VisualisationEvent> interpretAODTracks(TFile* AODFile, Int_t event);

  std::unique_ptr<VisualisationEvent> interpretAODCaloCells(TFile* AODFile, Int_t eventID);
  std::unique_ptr<VisualisationEvent> interpretEMCALCell(Int_t absID, Float_t amplitude, std::unique_ptr<VisualisationEvent> event);
  std::unique_ptr<VisualisationEvent> interpretPHOSCell(Int_t absID, Float_t amplitude, std::unique_ptr<VisualisationEvent> event);

  std::unique_ptr<VisualisationEvent> interpretMuonTracks(TFile* AODFile, Int_t eventID);

  Bool_t cutCell(Int_t absID, Float_t amplitude, Int_t type);

  // Calorimeter geometries
  o2::emcal::Geometry* mEMCALGeom = nullptr;
  o2::phos::Geometry* mPHOSGeom = nullptr;

  // From AliRoot emcal_esdclustercells.C macro
  // Cluster cuts          { PHOS, EMCAL }
  Int_t mNumMinCellsCut[2] = { 3, 2 };    /// Number of cells in cluster must be larger than this value.
  Int_t mNumMaxCellsCut[2] = { 60, 30 };  /// Number of cells in cluster must be smaller than this value.
  Float_t mEnergyCut[2] = { 0.30, 0.30 }; /// Cluster energy must be larger than this value.
  Float_t mM02LowCut[2] = { 0.70, 0.10 }; /// Cluster shower shape major axis must be larger than this value.
  Float_t mM02HigCut[2] = { 7.00, 7.00 }; /// Cluster shower shape major axis must be smaller than this value.
  Float_t mM20LowCut[2] = { 0.50, -1.0 }; /// Cluster shower shape lower axis must be larger than this value.
  Float_t mExoCut = 0.95;                 /// Reject clusters with this exoticity value.
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_DETECTORS_DATAINTERPRETERAOD_H
