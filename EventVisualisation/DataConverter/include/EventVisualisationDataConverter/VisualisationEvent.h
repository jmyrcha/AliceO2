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
/// \file    VisualisationEvent.h
/// \author  Jeremi Niedziela
/// \author  Maciej Grochowicz
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#ifndef ALICE_O2_EVENTVISUALISATION_DATACONVERTER_VISUALISATIONEVENT_H
#define ALICE_O2_EVENTVISUALISATION_DATACONVERTER_VISUALISATIONEVENT_H

#include "EventVisualisationDataConverter/VisualisationCluster.h"
#include "EventVisualisationDataConverter/VisualisationTrack.h"
#include "EventVisualisationDataConverter/VisualisationCaloCell.h"

#include <TEveTrack.h>

#include <vector>
#include <ctime>
#include <string>

namespace o2
{
namespace event_visualisation
{

/// Minimalistic description of an event
///
/// This class is used mainly for visualisation purposes.
/// It stores simple information about tracks, V0s, kinks, cascades,
/// clusters and calorimeter towers, which can be used for visualisation
/// or exported for external applications.

class VisualisationEvent
{
 public:
  // Default constructor
  VisualisationEvent(int eventNumber, int runNumber, double energy, int multiplicity, std::string collidingSystem, time_t timeStamp);

  // Adds minimalistic track to a minimalistic event
  void addTrack(const VisualisationTrack& track) { mTracks.push_back(track); }
  // Adds minimalistic muon track to a minimalistic event
  void addMuonTrack(const VisualisationTrack& track) { mMuonTracks.push_back(track); }
  // Adds minimalistic cluster to a minimalistic event
  void addCluster(const VisualisationCluster& cluster) { mClusters.push_back(cluster); }
  // Adds minimalistic calorimeter cell to a minimalistic event
  void addCaloCell(const VisualisationCaloCell& caloCell) { mCaloCells.push_back(caloCell); }

  // Event number getter
  inline int getEventNumber() { return mEventNumber; }
  // Run number getter
  inline int getRunNumber() { return mRunNumber; }
  // Energy getter
  inline double getEnergy() { return mEnergy; }
  // Multiplicity getter
  inline int getMultiplicity() { return mMultiplicity; }
  // Colliding system getter
  inline std::string getCollidingSystem() { return mCollidingSystem; }
  // TimeStamp getter
  inline std::time_t getTimeStamp() { return mTimeStamp; }

  // Returns track with index i
  const VisualisationTrack& getTrack(int i) const;
  // Returns number of tracks
  size_t getTrackCount() const { return mTracks.size(); }

  // Returns muon track with index i
  const VisualisationTrack& getMuonTrack(int i) const;
  // Returns number of muon tracks
  size_t getMuonTrackCount() const { return mMuonTracks.size(); }

  // Returns cluster with index i
  const VisualisationCluster& getCluster(int i) const;
  // Returns number of clusters
  size_t getClusterCount() const { return mClusters.size(); }

  // Returns calorimeter cell with index i
  const VisualisationCaloCell& getCaloCell(int i) const;
  // Returns number of calorimeter cells
  size_t getCaloCellsCount() const { return mCaloCells.size(); }

 private:
  int mEventNumber;                              /// event number in file
  int mRunNumber;                                /// run number
  double mEnergy;                                /// energy of the collision
  int mMultiplicity;                             /// number of particles reconstructed
  std::string mCollidingSystem;                  /// colliding system (e.g. proton-proton)
  std::time_t mTimeStamp;                        /// collision timestamp
  std::vector<VisualisationTrack> mTracks;       /// an array of minimalistic tracks
  std::vector<VisualisationTrack> mMuonTracks;   /// an array of minimalistic muon tracks
  std::vector<VisualisationCaloCell> mCaloCells; /// an array of minimalistic calorimeter cells
  std::vector<VisualisationCluster> mClusters;   /// an array of minimalistic clusters
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_DATACONVERTER_VISUALISATIONEVENT_H
