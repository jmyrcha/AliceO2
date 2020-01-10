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
/// \file    DataReader.h
/// \brief   Abstract base class for detector-specific data reading
/// \author  julian.myrcha@cern.ch
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#ifndef ALICE_O2_EVENTVISUALISATION_BASE_DATAREADER_H
#define ALICE_O2_EVENTVISUALISATION_BASE_DATAREADER_H

#include "EventVisualisationBase/VisualisationConstants.h"

class TObject;
class TList;

namespace o2
{
namespace event_visualisation
{

class DataReader
{
 public:
  virtual ~DataReader() = default;
  virtual TObject* getEventData(int eventNumber, EVisualisationDataType dataType, EDataSource source) = 0;
  virtual void setOnlineEventData(TList* data, EVisualisationDataType type) = 0;
  virtual void open() = 0;

  int getEventCount() { return mEventCount; };
  void setEventCount(int eventCount) { mEventCount = eventCount; };
  bool hasEventData(int eventNumber)
  {
    return (eventNumber >= 0 && eventNumber < mEventCount);
  }

 private:
  int mEventCount;
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_BASE_DATAREADER_H
