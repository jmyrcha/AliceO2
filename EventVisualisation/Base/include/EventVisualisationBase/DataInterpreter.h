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
/// \file    DataInterpreter.h
/// \author  Jeremi Niedziela
/// \author  julian.myrcha@cern.ch
/// \author  p.nowakowski@cern.ch
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#ifndef ALICE_O2_EVENTVISUALISATION_BASE_DATAINTERPRETER_H
#define ALICE_O2_EVENTVISUALISATION_BASE_DATAINTERPRETER_H

#include "EventVisualisationBase/VisualisationConstants.h"
#include "EventVisualisationDataConverter/VisualisationEvent.h"

#include <TEveElement.h>
#include <TFile.h>

namespace o2
{
namespace event_visualisation
{

/// DataInterpreter is a template class for detector-specific visualisation classes.
///
/// Each detector should override this template class, providing method(s)
/// to interpret detector-specific data (such as hits, digits, clusters)
/// as a set of visualisation objects (points, lines, boxes).

class DataInterpreter
{
 public:
  // Default constructor
  DataInterpreter() = default;
  // Virtual destructor
  virtual ~DataInterpreter() = default;

  // Should return visualisation objects for required data type
  virtual void interpretDataForType(TObject* data, EVisualisationDataType type, VisualisationEvent& event)
  {
    TList* list = (TList*)data;
    Int_t eventId = ((TVector2*)list->At(0))->X();

    if (type == Tracks) {
      TFile* file = (TFile*)list->At(1);
      interpretTracks(file, eventId, event);
    } else if (type == Clusters) {
      TFile* file = (TFile*)list->At(2);
      interpretClusters(file, eventId, event);
    }
  }

 private:
  virtual void interpretTracks(TFile* file, Int_t eventId, VisualisationEvent& event){};
  virtual void interpretClusters(TFile* file, Int_t eventId, VisualisationEvent& event){};
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_BASE_DATAINTERPRETER_H
