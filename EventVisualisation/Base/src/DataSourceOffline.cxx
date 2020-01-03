// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DataSourceOffline.h
/// \brief Grouping reading from file(s)
/// \author julian.myrcha@cern.ch
/// \author p.nowakowski@cern.ch

#include "EventVisualisationBase/DataSourceOffline.h"

namespace o2
{
namespace event_visualisation
{

DataReader* DataSourceOffline::sInstance[EVisualisationGroup::NvisualisationGroups];

TObject* DataSourceOffline::getEventData(int no, EVisualisationGroup purpose)
{
  if (sInstance[purpose] == nullptr)
    return nullptr;
  return sInstance[purpose]->getEventData(no);
}

int DataSourceOffline::getEventCount()
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
