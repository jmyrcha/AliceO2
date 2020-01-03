// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DataInterpreterTPC.h
/// \brief Converting TPC data to Event Visualisation primitives
/// \author julian.myrcha@cern.ch
/// \author p.nowakowski@cern.ch

#ifndef ALICE_O2_EVENTVISUALISATION_DETECTORS_DATAINTERPRETERTPC_H
#define ALICE_O2_EVENTVISUALISATION_DETECTORS_DATAINTERPRETERTPC_H

///
/// This class overrides DataInterpreter and implements method
/// returning visualisation objects representing data from TPC file.

#include "EventVisualisationBase/DataInterpreter.h"
#include "EventVisualisationBase/VisualisationConstants.h"

namespace o2
{
namespace event_visualisation
{

class DataInterpreterTPC : public DataInterpreter
{
 private:
  Int_t mTPCReadoutCycle = 100; // ms, provisional

 public:
  // Default constructor
  DataInterpreterTPC();

  // Default destructor
  ~DataInterpreterTPC() final;

  void interpretDataForType(TObject* data, EVisualisationDataType type, VisualisationEvent& event) final;
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_DETECTORS_DATAINTERPRETERTPC_H
