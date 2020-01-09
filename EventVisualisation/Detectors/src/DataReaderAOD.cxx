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

TObject* DataReaderAOD::getEventData(int eventNumber, EVisualisationDataType dataType)
{
  if (!this->hasEventData(eventNumber)) {
    return new TList();
  }

  /// FIXME: Redesign the data reader class
  TList* list = new TList();
  list->Add(this->mAODFile);
  TVector2* v = new TVector2(eventNumber, 0);
  list->Add(v);
  return list;
}
} // namespace event_visualisation
} // namespace o2
