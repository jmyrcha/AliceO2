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
/// \file   DataReaderAOD.cxx
/// \brief  AOD Detector-specific reading from file(s)
/// \author Maja Kabus

#include "EventVisualisationDetectors/DataReaderAOD.h"
#include "ReconstructionDataFormats/Track.h"

#include <TTree.h>
#include <TVector2.h>
#include <TError.h>

namespace o2
{
namespace event_visualisation
{

DataReaderAOD::DataReaderAOD() = default;

void DataReaderAOD::open()
{
  TString file = "O2aod.root";
  this->mAODFile = TFile::Open(file);

  TTree* trec = static_cast<TTree*>(this->mAODFile->Get("O2tracks"));
  TTree* calo = static_cast<TTree*>(this->mAODFile->Get("O2calo"));
  TTree* muon = static_cast<TTree*>(this->mAODFile->Get("O2mu"));

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

  for (int i = numTrackEvents - 1; i > 0; i--) {
    trec->GetEntry(i);
    if (trackEventID > maxEventID)
      maxEventID = trackEventID;
    // Some heuristic - assuming that events are in some order, perhaps a bit intertwined
    if (trackEventID < maxEventID - 10)
      break;
  }

  for (int i = numCaloEvents - 1; i > 0; i--) {
    calo->GetEntry(i);
    if (caloEventID > maxEventID)
      maxEventID = caloEventID;
    // Some heuristic - assuming that events are in some order, perhaps a bit intertwined
    if (caloEventID < maxEventID - 10)
      break;
  }

  for (int i = numMuonEvents - 1; i > 0; i--) {
    muon->GetEntry(i);
    if (muonEventID > maxEventID)
      maxEventID = muonEventID;
    // Some heuristic - assuming that events are in some order, perhaps a bit intertwined
    if (muonEventID < maxEventID - 10)
      break;
  }

  // TODO: Slow method. How to get events number?
  std::cout << "Setting max ev to: " << maxEventID << std::endl;
  fMaxEv = maxEventID;
}

Int_t DataReaderAOD::GetEventCount()
{
  return fMaxEv;
}

TObject* DataReaderAOD::getEventData(int no)
{
  /// FIXME: Redesign the data reader class
  TList* list = new TList();
  list->Add(this->mAODFile);
  TVector2* v = new TVector2(no, 0);
  list->Add(v);
  return list;
}
} // namespace event_visualisation
} // namespace o2
