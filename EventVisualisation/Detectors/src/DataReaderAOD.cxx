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
/// \file    DataReaderAOD.cxx
/// \brief   AOD detector-specific reading from file(s)
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#include "EventVisualisationDetectors/DataReaderAOD.h"

#include "FairLogger.h"

#include <TSystem.h>
#include <TTree.h>
#include <TVector2.h>

namespace o2
{
namespace event_visualisation
{

DataReaderAOD::DataReaderAOD()
{
  mTracks = new TList();
	mCaloCells = new TList();
	mMuonTracks = new TList();
}

void DataReaderAOD::open()
{
  TString file = "O2aod.root";

  this->mAODFile = TFile::Open(file);
  if (!this->mAODFile) {
    LOG(FATAL) << "There is no " << file.Data() << " file in current directory!";
  }

  TTree* trec = dynamic_cast<TTree*>(this->mAODFile->Get("O2tracks"));
  TTree* calo = dynamic_cast<TTree*>(this->mAODFile->Get("O2calo"));
  TTree* muon = dynamic_cast<TTree*>(this->mAODFile->Get("O2mu"));

  if (!trec || !calo || !muon) {
    LOG(FATAL) << "Incorrect AOD file format, branch missing!";
  }

  Int_t trackEventID;
  Int_t caloEventID;
  Int_t muonEventID;
  trec->SetBranchAddress("fID4Tracks", &trackEventID);
  calo->SetBranchAddress("fID4Calo", &caloEventID);
  muon->SetBranchAddress("fID4mu", &muonEventID);
  Int_t maxEventID = 0;

  int numTrackEvents = trec->GetEntriesFast();
  int numCaloEvents = calo->GetEntriesFast();
  int numMuonEvents = muon->GetEntriesFast();

  for (int i = numTrackEvents - 1; i >= 0; i--) {
    trec->GetEntry(i);
    if (trackEventID > maxEventID)
      maxEventID = trackEventID;
    // Some heuristic - assuming that events are in some order, perhaps a bit intertwined
    if (trackEventID < maxEventID - 10)
      break;
  }

  for (int i = numCaloEvents - 1; i >= 0; i--) {
    calo->GetEntry(i);
    if (caloEventID > maxEventID)
      maxEventID = caloEventID;
    // Some heuristic - assuming that events are in some order, perhaps a bit intertwined
    if (caloEventID < maxEventID - 10)
      break;
  }

  for (int i = numMuonEvents - 1; i >= 0; i--) {
    muon->GetEntry(i);
    if (muonEventID > maxEventID)
      maxEventID = muonEventID;
    // Some heuristic - assuming that events are in some order, perhaps a bit intertwined
    if (muonEventID < maxEventID - 10)
      break;
  }

  // TODO: Slow method. How to get events number?
  LOG(INFO) << "Setting max ev to: " << maxEventID + 1;
  setEventCount(maxEventID + 1);
}

TObject* DataReaderAOD::getEventData(int eventNumber, EVisualisationDataType dataType, EDataSource source)
{
  if (!this->hasEventData(eventNumber)) {
    return new TList();
  }

  if (dataType == Tracks) {
    if (source == EDataSource::SourceOffline) {
      return getTracks(eventNumber);
    }
    if (source == EDataSource::SourceOnline) {
      return getOnlineTracks(eventNumber);
    }
  }
  if (dataType == Calo) {
    if (source == EDataSource::SourceOffline) {
      return getCaloCells(eventNumber);
    }
    if (source == EDataSource::SourceOnline) {
      return getOnlineCaloCells(eventNumber);
    }
  }
  if (dataType == Muon) {
    if (source == EDataSource::SourceOffline) {
      return getMuonTracks(eventNumber);
    }
    if (source == EDataSource::SourceOnline) {
      return getOnlineMuonTracks(eventNumber);
    }
  }

  return new TList();
}

TObject* DataReaderAOD::getTracks(int eventId)
{
  TTree* tracks = (TTree*)mAODFile->Get("O2tracks");

  // Read all tracks parameters to buffers
  AODTrack fileTrack;
  tracks->SetBranchAddress("fID4Tracks", &fileTrack.mId);
  tracks->SetBranchAddress("fX", &fileTrack.mX);
  tracks->SetBranchAddress("fAlpha", &fileTrack.mAlpha);
  tracks->SetBranchAddress("fY", &fileTrack.mY);
  tracks->SetBranchAddress("fZ", &fileTrack.mZ);
  tracks->SetBranchAddress("fSnp", &fileTrack.mSnp);
  tracks->SetBranchAddress("fTgl", &fileTrack.mTgl);
  tracks->SetBranchAddress("fSigned1Pt", &fileTrack.mSigned1Pt);
  tracks->SetBranchAddress("fPID", &fileTrack.mPID);
  tracks->SetBranchAddress("fFlags", &fileTrack.mFlags);

  TList* fileTracks = new TList();
  int tracksCount = tracks->GetEntriesFast();
  for (int i = 0; i < tracksCount; i++) {
    tracks->GetEntry(i);

    // TODO: This is provisional.
    //  AOD is supposed to contain references from primary vertex to elements from given collision.
    if (fileTrack.mId > eventId + 10)
      break; // End event search
    if (fileTrack.mId != eventId)
      continue;

    fileTracks->Add(new AODTrack(fileTrack));
  }

  return fileTracks;
}

TObject* DataReaderAOD::getCaloCells(int eventId)
{
  TTree* calo = (TTree*)mAODFile->Get("O2calo");

  AODCalo fileCalo;
  calo->SetBranchAddress("fID4Calo", &fileCalo.mId);
  calo->SetBranchAddress("fCellNumber", &fileCalo.mCellNumber);
  calo->SetBranchAddress("fAmplitude", &fileCalo.mAmplitude);
  calo->SetBranchAddress("fType", &fileCalo.mType);

  TList* fileCaloList = new TList();
  int caloCount = calo->GetEntriesFast();
  for (int i = 0; i < caloCount; i++) {
    calo->GetEntry(i);

    // TODO: This is provisional.
    //  AOD is supposed to contain references from primary vertex to elements from given collision.
    if (fileCalo.mId > eventId + 10)
      break; // End event search
    if (fileCalo.mId != eventId)
      continue;

    fileCaloList->Add(new AODCalo(fileCalo));
  }

  return fileCaloList;
}

TObject* DataReaderAOD::getMuonTracks(int eventId)
{
  TTree* muon = (TTree*)mAODFile->Get("O2mu");

  AODMuonTrack fileMuon;
  muon->SetBranchAddress("fID4mu", &fileMuon.mId);
  muon->SetBranchAddress("fInverseBendingMomentum", &fileMuon.mInverseBendingMomentum);
  muon->SetBranchAddress("fThetaX", &fileMuon.mThetaX);
  muon->SetBranchAddress("fThetaY", &fileMuon.mThetaY);
  muon->SetBranchAddress("fZ", &fileMuon.mZ);
  muon->SetBranchAddress("fBendingCoor", &fileMuon.mBendingCoor);
  muon->SetBranchAddress("fNonBendingCoor", &fileMuon.mNonBendingCoor);

  TList* fileMuonList = new TList();
  int muonCount = muon->GetEntriesFast();
  for (int i = 0; i < muonCount; i++) {
    muon->GetEntry(i);

    // TODO: This is provisional.
    //  AOD is supposed to contain references from primary vertex to elements from given collision.
    if (fileMuon.mId > eventId + 10)
      break; // End event search
    if (fileMuon.mId != eventId)
      continue;

    fileMuonList->Add(new AODMuonTrack(fileMuon));
  }

  return fileMuonList;
}

TObject* DataReaderAOD::getOnlineTracks(int eventNumber)
{
  TList* onlineTracks = new TList();

  for (auto obj : *mTracks) {
    AODTrack* track = (AODTrack*)obj;
    if (track->mId > eventNumber + 10) {
      break; // End event search
    }

    if (track->mId == eventNumber)
      onlineTracks->Add(track);
  }

  return onlineTracks;
}

TObject* DataReaderAOD::getOnlineCaloCells(int eventNumber)
{
  TList* onlineCalo = new TList();

  for (auto obj : *mCaloCells) {
    AODCalo* calo = (AODCalo*)obj;
    if (calo->mId > eventNumber + 10)
      break; // End event search
    if (calo->mId != eventNumber)
      continue;

    onlineCalo->Add(calo);
  }

  return onlineCalo;
}

TObject* DataReaderAOD::getOnlineMuonTracks(int eventNumber)
{
  TList* onlineMuonTracks = new TList();

  for (auto obj : *mMuonTracks) {
    AODMuonTrack* muon = (AODMuonTrack*)obj;
    if (muon->mId > eventNumber + 10)
      break; // End event search
    if (muon->mId != eventNumber)
      continue;

    onlineMuonTracks->Add(muon);
  }

  return onlineMuonTracks;
}

void DataReaderAOD::setOnlineEventData(TList* data, EVisualisationDataType type)
{
  LOG(INFO) << "Setting online data of size: " << data->GetEntries();
  if (type == Tracks) {
    mTracks = data;
  } else if (type == Calo) {
    mCaloCells = data;
  } else if (type == Muon) {
    mMuonTracks = data;
  }
}

} // namespace event_visualisation
} // namespace o2
