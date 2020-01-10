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
/// \file    DataSourceOnline.h
/// \brief   Reading event data from DPL / Framework
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#ifndef ALICE_O2_EVENTVISUALISATION_BASE_DATASOURCEONLINE_H
#define ALICE_O2_EVENTVISUALISATION_BASE_DATASOURCEONLINE_H

#include "EventVisualisationBase/DataSource.h"

#include <TList.h>
#include <TObject.h>

namespace o2
{
namespace event_visualisation
{

class DataSourceOnline : public DataSource
{
 protected:
  static DataReader* sInstance[EVisualisationGroup::NvisualisationGroups];

 public:
  DataSourceOnline() = default;

  ~DataSourceOnline() override = default;
  DataSourceOnline(DataSourceOnline const&) = delete;

  /// Deleted assigment operator
  void operator=(DataSourceOnline const&) = delete;

  int getEventCount() override;
  void setEventCount(int eventCount, EVisualisationGroup detector)
  {
    sInstance[detector]->setEventCount(eventCount);
  }

  void registerReader(DataReader* reader, EVisualisationGroup purpose) override
  {
    sInstance[purpose] = reader;
  }

  bool hasEventData(int eventNumber, EVisualisationGroup detector) override;
  TObject* getEventData(int no, EVisualisationGroup detector, EVisualisationDataType dataType) override;
  void setOnlineEventData(TList* data, EVisualisationGroup detector, EVisualisationDataType type)
  {
    if (sInstance[detector]) {
      sInstance[detector]->setOnlineEventData(data, type);
    }
  }
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_BASE_DATASOURCEONLINE_H
