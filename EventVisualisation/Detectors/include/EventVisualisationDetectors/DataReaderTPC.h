// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DataReaderTPC.h
/// \brief TPC Detector-specific reading from file(s)
/// \author julian.myrcha@cern.ch

#ifndef ALICE_O2_EVENTVISUALISATION_DETECTORS_DATAREADERTPC_H
#define ALICE_O2_EVENTVISUALISATION_DETECTORS_DATAREADERTPC_H

#include "EventVisualisationBase/DataReader.h"

#include <TFile.h>

namespace o2
{
namespace event_visualisation
{

class DataReaderTPC : public DataReader
{
 private:
  Int_t mMaxEv;
  TFile* mClusFile;
  TFile* mTracFile;
  Int_t mTPCReadoutCycle = 100; // ms, provisional

 public:
  DataReaderTPC();
  void open() override;
  int getEventCount() const override;
  TObject* getEventData(int eventNumber) override;
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_DETECTORS_DATAREADERTPC_H
