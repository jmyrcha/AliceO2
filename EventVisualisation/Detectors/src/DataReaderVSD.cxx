// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DataReaderVSD.cxx
/// \brief VSD specific reading from file(s) (Visualisation Summary Data)
/// \author julian.myrcha@cern.ch
/// \author p.nowakowski@cern.ch

#include "EventVisualisationDetectors/DataReaderVSD.h"

#include "FairLogger.h"

#include <TSystem.h>
#include <TFile.h>
#include <TPRegexp.h>
#include <TEveEventManager.h>

namespace o2
{
namespace event_visualisation
{

DataReaderVSD::DataReaderVSD()
  : DataReader(),
    mFile(nullptr),
    mMaxEv(-1),
    mCurEv(-1)
{
}

DataReaderVSD::~DataReaderVSD()
{
  delete mFile;
}

void DataReaderVSD::open()
{
  TString ESDFileName = "events_0.root";
  Warning("GotoEvent", "OPEN");
  mMaxEv = -1;
  mCurEv = -1;
  mFile = TFile::Open(ESDFileName);
  if (!mFile) {
    LOG(FATAL) << "There is no " << ESDFileName.Data() << " file in current directory!";
  }

  assert(mEvDirKeys.size() == 0);

  TPMERegexp name_re("Event\\d+");
  TObjLink* lnk = mFile->GetListOfKeys()->FirstLink();
  while (lnk) {
    if (name_re.Match(lnk->GetObject()->GetName())) {
      mEvDirKeys.push_back((TKey*)lnk->GetObject());
    }
    lnk = lnk->Next();
  }

  mMaxEv = mEvDirKeys.size();
  if (mMaxEv == 0) {
    LOG(FATAL) << "No VSD events to show.";
  }
}

TObject* DataReaderVSD::getEventData(int ev)
{
  if (ev < 0 || ev >= this->mMaxEv) {
    LOG(WARNING) << "Invalid VSD event id " << ev;
    return nullptr;
  }
  this->mCurEv = ev;
  return this->mEvDirKeys[this->mCurEv]->ReadObj();
}

} // namespace event_visualisation
} // namespace o2
