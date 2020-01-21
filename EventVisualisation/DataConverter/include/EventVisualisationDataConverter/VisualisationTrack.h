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
/// \file    VisualisationTrack.h
/// \author  Jeremi Niedziela
/// \author  Maciej Grochowicz
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#ifndef ALICE_O2_EVENTVISUALISATION_DATACONVERTER_VISUALISATIONTRACK_H
#define ALICE_O2_EVENTVISUALISATION_DATACONVERTER_VISUALISATIONTRACK_H

#include "ConversionConstants.h"
#include "ReconstructionDataFormats/PID.h"

#include <iosfwd>
#include <string>
#include <vector>
#include <array>
#include <cmath>

namespace o2
{
namespace event_visualisation
{

/// Minimalistic description of a particle's track
///
/// This class is used mainly for visualisation purpose.
/// It keeps basic information about a track, such as its vertex,
/// momentum, PID, phi and theta or helix curvature.

class VisualisationTrack
{
 public:
  // Default constructor
  VisualisationTrack();

  // Constructor with properties initialisation
  VisualisationTrack(
    int ID,
    int type,
    int charge,
    double energy,
    int parentID,
    o2::track::PID PID,
    double signedPT,
    double mass,
    double pxpypz[],
    double startXYZ[],
    double endXYZ[],
    double helixCurvature,
    double theta,
    double phi,
    float C1Pt21Pt2,
    unsigned long long flags);

  // Add child particle (coming from decay of this particle)
  void addChild(int childID);
  // Add xyz coordinates of the point along the track
  void addPolyPoint(double x, double y, double z);
  // Add xyz coordinates of the point along the track
  void addPolyPoint(double xyz[3]);
  // Track type setter (standard track, V0, kink, cascade)
  void setTrackType(ETrackType type);

  // ID getter
  inline int getID() const { return mID; }
  // Type getter
  inline std::string getType() const { return mType; }
  // Charge getter
  inline int getCharge() const { return mCharge; }
  // Energy getter
  inline double getEnergy() const { return mEnergy; }
  // Parent ID getter
  inline int getParentID() const { return mParentID; }
  // PID (particle identification code) getter
  inline o2::track::PID getPID() const { return mPID; }
  // Signed transverse momentum getter
  inline double getSignedPt() const { return mSignedPT; }
  // Mass getter
  inline double getMass() const { return mMass; }
  // Momentum vector getter
  inline double* getMomentum() { return mMomentum; }
  // Vertex getter
  inline double* getVertex() { return mStartCoordinates; }
  // Helix curvature getter
  inline double getHelixCurvature() const { return mHelixCurvature; }
  // Theta getter
  inline double getTheta() const { return mTheta; }
  // Phi getter
  inline double getPhi() const { return mPhi; }
  // Momentum sigma getter
  inline float getC1Pt21Pt2() const { return mC1Pt21Pt2; }
  // Reconstruction flags getter
  inline unsigned long long getFlags() const { return mFlags; }

  // Beta (velocity) getter
  inline double getBeta() const { return sqrt(1 - std::pow(mMass / mEnergy, 2)); }
  // Checks if there was specific flag set for reconstruction
  inline bool isRecoFlagSet(unsigned long long mask) const { return (mFlags & mask & 0xffff) > 0; }

  size_t getPointCount() { return mPolyX.size(); }
  std::array<double, 3> getPoint(size_t i) { return std::array<double, 3>{mPolyX[i], mPolyY[i], mPolyZ[i]}; }

 private:
  // Set coordinates of the beginning of the track
  void setStartCoordinates(double xyz[3]);
  // Set coordinates of the end of the track
  void setEndCoordinates(double xyz[3]);
  /// Set momentum vector
  void setMomentum(double pxpypz[3]);

  int mID;                     /// Unique identifier of the track
  std::string mType;           /// Type (standard, V0 mother, daughter etc.)
  int mCharge;                 /// Charge of the particle
  double mEnergy;              /// Energy of the particle
  int mParentID;               /// ID of the parent-track (-1 means no parent)
  o2::track::PID mPID;         /// Particle ID
  double mSignedPT;            /// Signed transverse momentum
  double mMass;                /// Mass of the particle
  double mMomentum[3];         /// Momentum vector
  double mStartCoordinates[3]; /// Vector of track's start coordinates
  double mEndCoordinates[3];   /// Vector of track's end coordinates
  double mHelixCurvature;      /// Helix curvature of the trajectory
  double mTheta;               /// An angle from Z-axis to the radius vector pointing to the particle
  double mPhi;                 /// An angle from X-axis to the radius vector pointing to the particle
  float mC1Pt21Pt2;            /// Momentum sigma, 14th element of covariance matrix
  unsigned long long mFlags;   /// Reconstruction status flags - track quality parameters

  std::vector<int> mChildrenIDs; /// Uniqe IDs of children particles

  /// Polylines -- array of points along the trajectory of the track
  std::vector<double> mPolyX;
  std::vector<double> mPolyY;
  std::vector<double> mPolyZ;
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_DATACONVERTER_VISUALISATIONTRACK_H
