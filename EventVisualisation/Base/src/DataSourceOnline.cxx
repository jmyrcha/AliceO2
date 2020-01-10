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
/// \file    DataSourceOnline.cxx
/// \brief   Reading event data from DPL / Framework
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#include "EventVisualisationBase/DataSourceOnline.h"

namespace o2
{
namespace event_visualisation
{

DataReader* DataSourceOnline::sInstance[EVisualisationGroup::NvisualisationGroups];

TObject* DataSourceOnline::getEventData(int no, EVisualisationGroup purpose, EVisualisationDataType dataType)
{
  if (sInstance[purpose] == nullptr) {
    return new TList();
  }
  return sInstance[purpose]->getEventData(no, dataType, SourceOnline);
}

bool DataSourceOnline::hasEventData(int eventNumber, EVisualisationGroup detector)
{
  if (sInstance[detector] == nullptr)
    return false;
  return sInstance[detector]->hasEventData(eventNumber);
}

int DataSourceOnline::getEventCount()
{
  int eventCount = 0;
  for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; i++) {
    if (sInstance[i] != nullptr && sInstance[i]->getEventCount() > eventCount) {
      eventCount = sInstance[i]->getEventCount();
    }
  }
  return eventCount;
};

} // namespace event_visualisation
} // namespace o2
