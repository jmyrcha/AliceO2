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
/// \file    EventManager.cxx
/// \author  Jeremi Niedziela
/// \author julian.myrcha@cern.ch
/// \author p.nowakowski@cern.ch

#include "EventVisualisationView/EventManager.h"
#include "EventVisualisationDataConverter/VisualisationEvent.h"
#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationView/MultiView.h"
#include "EventVisualisationBase/DataSource.h"
#include "EventVisualisationBase/DataInterpreter.h"
#include "EventVisualisationBase/DataSourceOffline.h"
#include "EventVisualisationDetectors/DataReaderVSD.h"
#include "EventVisualisationDetectors/DataReaderITS.h"

#include <TEveManager.h>
#include <TEveProjectionManager.h>
#include <TEveTrackPropagator.h>
#include <TSystem.h>
#include <TEnv.h>
#include <TEveElement.h>
#include <TGListTree.h>

#include <iostream>

using namespace std;

namespace o2
{
namespace event_visualisation
{

EventManager* EventManager::mInstance = nullptr;

EventManager& EventManager::getInstance()
{
  if (mInstance == nullptr) {
    mInstance = new EventManager();
  }
  return *mInstance;
}

EventManager::EventManager() : TEveEventManager("Event", "")
{
  for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; i++) {
    mDataInterpreters[i] = nullptr;
    mDataReaders[i] = nullptr;
  }
}

void EventManager::Open()
{
  switch (mCurrentDataSourceType) {
    case SourceOnline:
      break;
    case SourceOffline: {
      DataSourceOffline* source = new DataSourceOffline();
      for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; i++) {
        if (mDataInterpreters[i] != nullptr) {
          mDataReaders[i]->open();
          source->registerReader(mDataReaders[i], static_cast<EVisualisationGroup>(i));
        }
      }
      setDataSource(source);
    } break;
    case SourceHLT:
      break;
  }
}

void EventManager::GotoEvent(Int_t no)
{
  //-1 means last event
  if (no == -1) {
    no = getDataSource()->GetEventCount() - 1;
  }

  this->mCurrentEvent = no;

  MultiView::getInstance()->destroyAllEvents();

  for (int i = 0; i < EVisualisationDataType::NdataTypes; ++i) {
    mDataTypeLists[i] = new TEveElementList(gDataTypeNames[i].c_str());
  }

  for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; ++i) {
    for (int dataType = 0; dataType < EVisualisationDataType::NdataTypes; ++dataType) {
      DataInterpreter* interpreter = mDataInterpreters[i];
      if (interpreter) {
        TObject* data = getDataSource()->getEventData(no, (EVisualisationGroup)i);
        std::unique_ptr<VisualisationEvent> event = interpreter->interpretDataForType(data, (EVisualisationDataType)dataType);
        displayVisualisationEvent(*event, gVisualisationGroupName[i]);
      }
    }
  }

  for (int i = 0; i < EVisualisationDataType::NdataTypes; ++i) {
    MultiView::getInstance()->registerElement(mDataTypeLists[i]);
  }

  MultiView::getInstance()->redraw3D();
}

void EventManager::NextEvent()
{
  Int_t event = (this->mCurrentEvent + 1) % getDataSource()->GetEventCount();
  GotoEvent(event);
}

void EventManager::PrevEvent()
{
  GotoEvent(this->mCurrentEvent - 1);
}

void EventManager::Close()
{
  delete this->mDataSource;
  this->mDataSource = nullptr;
}

void EventManager::AfterNewEventLoaded()
{
  TEveEventManager::AfterNewEventLoaded();
}

void EventManager::AddNewEventCommand(const TString& cmd)
{
  TEveEventManager::AddNewEventCommand(cmd);
}

void EventManager::RemoveNewEventCommand(const TString& cmd)
{
  TEveEventManager::RemoveNewEventCommand(cmd);
}

void EventManager::ClearNewEventCommands()
{
  TEveEventManager::ClearNewEventCommands();
}

EventManager::~EventManager()
{
  for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; i++) {
    if (mDataInterpreters[i] != nullptr) {
      delete mDataInterpreters[i];
      delete mDataReaders[i];

      mDataInterpreters[i] = nullptr;
      mDataReaders[i] = nullptr;
    }
  }
  mInstance = nullptr;
}

void EventManager::displayVisualisationEvent(VisualisationEvent& event, const std::string& detectorName)
{
  size_t trackCount = event.getTrackCount();

  auto* list = new TEveTrackList(detectorName.c_str());
  list->IncDenyDestroy();

  for (size_t i = 0; i < trackCount; ++i) {
    VisualisationTrack track = event.getTrack(i);
    TEveRecTrackD t;
    double* p = track.getMomentum();
    t.fP = { p[0], p[1], p[2] };
    t.fSign = track.getCharge() > 0 ? 1 : -1;
    auto* vistrack = new TEveTrack(&t, &TEveTrackPropagator::fgDefault);
    vistrack->SetLineColor(kMagenta);
    size_t pointCount = track.getPointCount();
    vistrack->Reset(pointCount);

    for (size_t j = 0; j < pointCount; ++j) {
      auto point = track.getPoint(j);
      vistrack->SetNextPoint(point[0], point[1], point[2]);
    }
    list->AddElement(vistrack);
  }

  if (trackCount != 0) {
    mDataTypeLists[EVisualisationDataType::ESD]->AddElement(list);
  }

  size_t clusterCount = event.getClusterCount();
  auto* point_list = new TEvePointSet(detectorName.c_str());
  point_list->IncDenyDestroy();
  point_list->SetMarkerColor(kBlue);

  for (size_t i = 0; i < clusterCount; ++i) {
    VisualisationCluster cluster = event.getCluster(i);
    point_list->SetNextPoint(cluster.X(), cluster.Y(), cluster.Z());
  }

  if (clusterCount != 0) {
    mDataTypeLists[EVisualisationDataType::Clusters]->AddElement(point_list);
  }
}

void EventManager::registerDetector(DataReader* reader, DataInterpreter* interpreter, EVisualisationGroup type)
{
  mDataReaders[type] = reader;
  mDataInterpreters[type] = interpreter;
}

} // namespace event_visualisation
} // namespace o2
