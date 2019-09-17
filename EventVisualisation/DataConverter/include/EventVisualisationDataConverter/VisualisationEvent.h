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
/// \author julian.myrcha@cern.ch
/// \author p.nowakowski@cern.ch
///

#ifndef ALICE_O2_DATACONVERTER_VISUALISATIONEVENT_H
#define ALICE_O2_DATACONVERTER_VISUALISATIONEVENT_H

#include "EventVisualisationDataConverter/VisualisationCluster.h"
#include "EventVisualisationDataConverter/VisualisationTrack.h"

#include <vector>
#include <ctime>

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

  // Adds minimalistic track inside minimalistic event
  void addTrack(const VisualisationTrack& track) { mTracks.push_back(track); }
  // Adds minimalistic cluster inside minimalistic event
  void addCluster(const VisualisationCluster& cluster) { mClusters.push_back(cluster); }

  // Multiplicity getter
  inline int GetMultiplicity() { return mMultiplicity; }

  // Returns track with index i
  const VisualisationTrack& getTrack(int i) const;
  // Returns number of tracks
  size_t getTrackCount() const { return mTracks.size(); }

  // Returns cluster with index i
  const VisualisationCluster& getCluster(int i) const;
  // Returns number of clusters
  size_t getClusterCount() const { return mClusters.size(); }

 private:
  int mEventNumber;                            /// event number in file
  int mRunNumber;                              /// run number
  double mEnergy;                              /// energy of the collision
  int mMultiplicity;                           /// number of particles reconstructed
  std::string mCollidingSystem;                /// colliding system (e.g. proton-proton)
  std::time_t mTimeStamp;                      /// collision timestamp
  std::vector<VisualisationTrack> mTracks;     /// an array of minimalistic tracks
  std::vector<VisualisationCluster> mClusters; /// an array of minimalistic clusters
};

} // namespace event_visualisation
} // namespace o2

#endif
