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
/// \file    DataSourceOffline.h
/// \brief   Group reading from file(s)
/// \author  julian.myrcha@cern.ch
/// \author  p.nowakowski@cern.ch
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#ifndef ALICE_O2_EVENTVISUALISATION_BASE_DATASOURCEOFFLINE_H
#define ALICE_O2_EVENTVISUALISATION_BASE_DATASOURCEOFFLINE_H

#include "EventVisualisationBase/DataSource.h"

#include <TObject.h>

namespace o2
{
namespace event_visualisation
{

class DataSourceOffline : public DataSource
{
 protected:
  static DataReader* sInstance[EVisualisationGroup::NvisualisationGroups];

 public:
  DataSourceOffline() = default;

  ~DataSourceOffline() override = default;
  DataSourceOffline(DataSourceOffline const&) = delete;

  /// Deleted assignment operator
  void operator=(DataSourceOffline const&) = delete;

  int getEventCount() override;

  void registerReader(DataReader* reader, EVisualisationGroup detector) override
  {
    sInstance[detector] = reader;
  }

  bool hasEventData(int eventNumber, EVisualisationGroup detector) override;
  TObject* getEventData(int no, EVisualisationGroup detector, EVisualisationDataType dataType) override;
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_BASE_DATASOURCEOFFLINE_H
