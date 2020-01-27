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
/// \file    DataSource.h
/// \brief   Group data reading
/// \author  julian.myrcha@cern.ch
/// \author  p.nowakowski@cern.ch
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#ifndef ALICE_O2_EVENTVISUALISATION_BASE_DATASOURCE_H
#define ALICE_O2_EVENTVISUALISATION_BASE_DATASOURCE_H

#include "EventVisualisationBase/VisualisationConstants.h"
#include "EventVisualisationBase/DataReader.h"

class TObject;

namespace o2
{
namespace event_visualisation
{

class DataSource
{
 public:
  virtual TObject* getEventData(int no, EVisualisationGroup detector, EVisualisationDataType dataType) = 0;
  virtual bool hasEventData(int eventNumber, EVisualisationGroup detector) { return false; }
  virtual int getEventCount() = 0;
  virtual void registerReader(DataReader* reader, EVisualisationGroup detector) = 0;

  DataSource() = default;

  /// Default destructor
  virtual ~DataSource() = default;

  /// Deleted copy constructor
  DataSource(DataSource const&) = delete;

  /// Deleted assignemt operator
  void operator=(DataSource const&) = delete;
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_BASE_DATASOURCE_H
