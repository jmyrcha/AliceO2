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
/// \file    DataReaderAOD.h
/// \brief   AOD detector-specific reading from file(s)
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#ifndef ALICE_O2_EVENTVISUALISATION_DETECTORS_DATAREADERAOD_H
#define ALICE_O2_EVENTVISUALISATION_DETECTORS_DATAREADERAOD_H

#include "EventVisualisationBase/DataReader.h"
#include "ReconstructionDataFormats/PID.h"

#include <TFile.h>
#include <TList.h>

namespace o2
{
namespace event_visualisation
{

class AODTrack : public TObject
{
 public:
  Int_t mId;          // The index of the collision vertex, to which the track is attached
  Float_t mX;         // X coordinate for the point of parametrisation
  Float_t mAlpha;     // Local <--> global coor.system rotation angle
  Float_t mY;         // fP[0] local Y-coordinate of a track (cm)
  Float_t mZ;         // fP[1] local Z-coordinate of a track (cm)
  Float_t mSnp;       // fP[2] local sine of the track momentum azimuthal angle
  Float_t mTgl;       // fP[3] tangent of the track momentum dip angle
  Float_t mSigned1Pt; // fP[4] 1/pt (1/(GeV/c))
  o2::track::PID mPID;
  ULong64_t mFlags; // Reconstruction status flags - track quality parameters

  AODTrack() = default;
  AODTrack(AODTrack& other) = default;
  AODTrack(Int_t id, Float_t x, Float_t alpha, Float_t y, Float_t z, Float_t snp, Float_t tgl, Float_t signed1Pt, o2::track::PID pid, ULong64_t flags) : mId(id), mX(x), mAlpha(alpha), mY(y), mZ(z), mSnp(snp), mTgl(tgl), mSigned1Pt(signed1Pt), mPID(pid), mFlags(flags) {}
};

class AODCalo : public TObject
{
 public:
  Int_t mId;           // The index of the collision vertex, to which the cell is attached
  Short_t mCellNumber; // Cell absolute Id. number
  Float_t mAmplitude;  // Cell amplitude (= energy!)
  Char_t mType;        // Cell type (-1 is undefined, 0 is PHOS, 1 is EMCAL)

  AODCalo() = default;
  AODCalo(AODCalo& other) = default;
  AODCalo(Int_t id, Short_t cellNumber, Float_t amplitude, Char_t type) : mId(id), mCellNumber(cellNumber), mAmplitude(amplitude), mType(type) {}
};

class AODMuonTrack : public TObject
{
 public:
  Int_t mId;                       // The index of the collision vertex, to which the muon is attached
  Float_t mInverseBendingMomentum; // Inverse bending momentum (GeV/c ** -1) times the charge
  Float_t mThetaX;                 // Angle of track at vertex in X direction (rad)
  Float_t mThetaY;                 // Angle of track at vertex in Y direction (rad)
  Float_t mZ;                      // Z coordinate (cm)
  Float_t mBendingCoor;            // bending coordinate (cm)
  Float_t mNonBendingCoor;         // non bending coordinate (cm)

  AODMuonTrack() = default;
  AODMuonTrack(AODMuonTrack& other) = default;
  AODMuonTrack(Int_t id, Float_t inverseBendingMomentum, Float_t thetaX, Float_t thetaY, Float_t z, Float_t bendingCoor, Float_t nonBendingCoor) : mId(id), mInverseBendingMomentum(inverseBendingMomentum), mThetaX(thetaX), mThetaY(thetaY), mZ(z), mBendingCoor(bendingCoor), mNonBendingCoor(nonBendingCoor) {}
};

class DataReaderAOD : public DataReader
{
 private:
  TFile* mAODFile;

  TList* mTracks;
  TList* mCaloCells;
  TList* mMuonTracks;

 public:
  DataReaderAOD();
  void open() final;

  void setOnlineEventData(TList* data, EVisualisationDataType type) final;

  TObject* getEventData(int eventNumber, EVisualisationDataType dataType, EDataSource source) final;
  TObject* getTracks(int eventNumber);
  TObject* getCaloCells(int eventNumber);
  TObject* getMuonTracks(int eventNumber);
  TObject* getOnlineTracks(int eventNumber);
  TObject* getOnlineCaloCells(int eventNumber);
  TObject* getOnlineMuonTracks(int eventNumber);
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_DETECTORS_DATAREADERAOD_H
