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
/// \file    DataReaderTPC.h
/// \brief   TPC detector-specific reading from file(s)
/// \author  julian.myrcha@cern.ch
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

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
  TFile* mClusFile;
  TFile* mTracFile;
  int mTPCReadoutCycle = 100; // ms, provisional

 public:
  DataReaderTPC() = default;
  void open() final;
  TObject* getEventData(int eventNumber, EVisualisationDataType dataType) final;
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_DETECTORS_DATAREADERTPC_H
