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
/// \file    EventManager.h
/// \author  Jeremi Niedziela
/// \author julian.myrcha@cern.ch
/// \author p.nowakowski@cern.ch

#ifndef ALICE_O2_EVENTVISUALISATION_VIEW_EVENTMANAGER_H
#define ALICE_O2_EVENTVISUALISATION_VIEW_EVENTMANAGER_H

#include "EventVisualisationDataConverter/VisualisationEvent.h"
#include "EventVisualisationBase/DataInterpreter.h"
#include "EventVisualisationBase/DataReader.h"

#include "CCDB/Manager.h"

#include <TEveElement.h>
#include <TEveEventManager.h>
#include <TGeoMatrix.h>
#include <TEveQuadSet.h>
#include <TEveTrackPropagator.h>

#include <string>

namespace o2
{
namespace event_visualisation
{

/// EventManager is a singleton class managing event loading.
///
/// This class is a hub for data visualisation classes, providing them with objects of requested type
/// (Raw data, hits, digits, clusters, ESDs, AODs...). It is a role of detector-specific data macros to
/// interpret data from different formats as visualisation objects (points, lines...) and register them
/// for drawing in the MultiView.

class DataSource;

class EventManager : public TEveEventManager
{
 public:
  enum EDataSource {
    SourceOnline,  ///< Online reconstruction is a source of events
    SourceOffline, ///< Local files are the source of events
    SourceHLT      ///< HLT reconstruction is a source of events
  };

  /// Returns an instance of EventManager
  static EventManager& getInstance();

  /// Setter of the current data source type
  inline void setDataSourceType(EDataSource source) { mCurrentDataSourceType = source; }

  /// Sets the CDB path in CCDB Manager
  inline void setCdbPath(const TString& path)
  {
    o2::ccdb::Manager::Instance()->setDefaultStorage(path.Data());
  }

  Int_t getCurrentEvent() const { return mCurrentEvent; }
  DataSource* getDataSource() const { return mDataSource; }

  void setDataSource(DataSource* dataSource) { this->mDataSource = dataSource; }

  void Open() override;
  void GotoEvent(Int_t /*event*/) override;
  void NextEvent() override;
  void PrevEvent() override;
  void Close() override;

  void AfterNewEventLoaded() override;

  void AddNewEventCommand(const TString& cmd) override;
  void RemoveNewEventCommand(const TString& cmd) override;
  void ClearNewEventCommands() override;

  void registerDetector(DataReader* reader, DataInterpreter* interpreter, EVisualisationGroup type);

 private:
  static EventManager* mInstance;
  DataInterpreter* mDataInterpreters[EVisualisationGroup::NvisualisationGroups];
  DataReader* mDataReaders[EVisualisationGroup::NvisualisationGroups];

  /// store lists of visualisation element in current event (row, clusters ...)
  TEveElementList* mDataTypeLists[EVisualisationDataType::NdataTypes];

  EDataSource mCurrentDataSourceType = EDataSource::SourceOffline;
  DataSource* mDataSource = nullptr;
  Int_t mCurrentEvent = 0;

  Width_t mWidth;

  /// Default constructor
  EventManager();
  /// Default destructor
  ~EventManager() final;
  /// Deleted copy constructor
  EventManager(EventManager const&) = delete;
  /// Deleted assignemt operator
  void operator=(EventManager const&) = delete;

  //Display visualisation event
  void displayVisualisationEvent(VisualisationEvent& event, const std::string& detectorName);

  void displayTracks(VisualisationEvent& event, const std::string& detectorName);
  void displayTracksByPt(VisualisationEvent& event, const std::string& detectorName);

  void displayMuonTracks(VisualisationEvent& event);
  void setupMuonTrackPropagator(TEveTrackPropagator* prop, Bool_t tracker, Bool_t trigger);
  void displayCalo(VisualisationEvent& event);
  void setCaloQuadSet(const Float_t quadSize, const TGeoHMatrix* matrix, TEveQuadSet* quadSet);

  void displayClusters(VisualisationEvent& event, const std::string& detectorName);
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_VIEW_EVENTMANAGER_H
