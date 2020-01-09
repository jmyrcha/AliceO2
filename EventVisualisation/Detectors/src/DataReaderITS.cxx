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
/// \file    DataReaderITS.cxx
/// \brief   ITS detector-specific reading from file(s)
/// \author  julian.myrcha@cern.ch
/// \author  p.nowakowski@cern.ch
///

#include "DataFormatsITSMFT/ROFRecord.h"
#include "EventVisualisationDetectors/DataReaderITS.h"

#include "FairLogger.h"

#include <TSystem.h>
#include <TTree.h>
#include <TVector2.h>

namespace o2
{
namespace event_visualisation
{

void DataReaderITS::open()
{
  TString clusterFile = "o2clus_its.root";
  TString trackFile = "o2trac_its.root";

  this->mTracFile = TFile::Open(trackFile);
  if (!this->mTracFile) {
    LOG(FATAL) << "There is no " << trackFile.Data() << " file in current directory!";
  }
  this->mClusFile = TFile::Open(clusterFile);
  if (!this->mClusFile) {
    LOG(FATAL) << "There is no " << clusterFile.Data() << " file in current directory!";
  }

  TTree* tracksRof = dynamic_cast<TTree*>(this->mTracFile->Get("ITSTracksROF"));
  TTree* clustersRof = dynamic_cast<TTree*>(this->mClusFile->Get("ITSClustersROF"));

  if (!tracksRof || !clustersRof) {
    LOG(FATAL) << "Incorrect ITS file format, branch missing!";
  }

  // Read all track RO frames to a buffer to count number of elements
  std::vector<o2::itsmft::ROFRecord>* trackROFrames = nullptr;
  tracksRof->SetBranchAddress("ITSTracksROF", &trackROFrames);
  tracksRof->GetEntry(0);

  // Read all cluster RO frames to a buffer
  std::vector<o2::itsmft::ROFRecord>* clusterROFrames = nullptr;
  clustersRof->SetBranchAddress("ITSClustersROF", &clusterROFrames);
  clustersRof->GetEntry(0);

  if (trackROFrames->size() == clusterROFrames->size()) {
    setEventCount(trackROFrames->size());
  } else {
    LOG(FATAL) << "DataReaderITS: Inconsistent number of readout frames in files";
  }
}

TObject* DataReaderITS::getEventData(int eventNumber, EVisualisationDataType dataType)
{
  if (!this->hasEventData(eventNumber)) {
    return new TList();
  }

  /// FIXME: Redesign the data reader class
  TList* list = new TList();
  list->Add(this->mTracFile);
  list->Add(this->mClusFile);
  TVector2* v = new TVector2(eventNumber, 0);
  list->Add(v);
  return list;
}
} // namespace event_visualisation
} // namespace o2
