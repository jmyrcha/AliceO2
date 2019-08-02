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
#include "EventVisualisationBase/DataSource.h"
#include "EventVisualisationBase/DataInterpreter.h"
#include <EventVisualisationBase/DataSourceOffline.h>
#include <EventVisualisationDetectors/DataReaderVSD.h>
#include <EventVisualisationDetectors/DataReaderITS.h>

#include <TEveManager.h>
#include <TEveProjectionManager.h>
#include <TEveTrackPropagator.h>
#include <TSystem.h>
#include <TEnv.h>
#include <TEveElement.h>
#include <TGListTree.h>

#include <iostream>
#include <EventVisualisationView/MultiView.h>

using namespace std;

namespace o2
{
namespace event_visualisation
{

EventManager* EventManager::instance = nullptr;

EventManager& EventManager::getInstance()
{
  if (instance == nullptr) {
    instance = new EventManager();
  }
  return *instance;
}

EventManager::EventManager() : TEveEventManager("Event", "")
{
  for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; i++) {
    dataInterpreters[i] = nullptr;
    dataReaders[i] = nullptr;
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
        if (dataInterpreters[i] != nullptr) {
          dataReaders[i]->open();
          source->registerReader(dataReaders[i], static_cast<EVisualisationGroup>(i));
        }
      }
      setDataSource(source);
    }
      break;
    case SourceHLT:
      break;
  }
}

void EventManager::GotoEvent(Int_t no) {
    //-1 means last event
    if(no == -1) {
        no = getDataSource()->GetEventCount()-1;
    }

    this->currentEvent = no;

    MultiView::getInstance()->destroyAllEvents();

    for (int i = 0; i < EVisualisationDataType::NdataTypes; ++i)
    {
        dataTypeLists[i] = new TEveElementList(gDataTypeNames[i].c_str());
    }

    for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; ++i) {
      for(int dataType = 0; dataType < EVisualisationDataType::NdataTypes; ++dataType) {
        DataInterpreter* interpreter = dataInterpreters[i];
        if(interpreter) {
            TObject *data = getDataSource()->getEventData(no, (EVisualisationGroup)i);
            std::unique_ptr<VisualisationEvent> event = interpreter->interpretDataForType(data, (EVisualisationDataType)dataType);
            displayVisualisationEvent(*event, gVisualisationGroupName[i]);
        }
      }
    }

    for (int i = 0; i < EVisualisationDataType::NdataTypes; ++i)
    {
        MultiView::getInstance()->registerElement(dataTypeLists[i]);
    }

    MultiView::getInstance()->redraw3D();
}

void EventManager::NextEvent()
{
  Int_t event = (this->currentEvent + 1) % getDataSource()->GetEventCount();
  GotoEvent(event);
}

void EventManager::PrevEvent()
{
  GotoEvent(this->currentEvent - 1);
}

void EventManager::Close()
{
  delete this->dataSource;
  this->dataSource = nullptr;
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


EventManager::~EventManager() {
    for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; i++)
    {
        if(dataInterpreters[i] != nullptr) {
            delete dataInterpreters[i];
            assert(dataReaders[i] != nullptr);
            delete dataReaders[i];

      dataInterpreters[i] = nullptr;
      dataReaders[i] = nullptr;
    }
  }
  instance = nullptr;
}


void EventManager::displayVisualisationEvent(VisualisationEvent &event, const std::string &detectorName) {
  size_t trackCount = event.getTrackCount();

  auto *list = new TEveTrackList(detectorName.c_str());
  list->IncDenyDestroy();

  for(size_t i = 0; i < trackCount; ++i) {
    VisualisationTrack track = event.getTrack(i);
    TEveRecTrackD t;
    double *p = track.getMomentum();
    t.fP = { p[0], p[1], p[2] };
    t.fSign = track.getCharge() > 0 ? 1 : -1;
    auto* vistrack = new TEveTrack(&t, &TEveTrackPropagator::fgDefault);
    vistrack->SetLineColor(kMagenta);
    size_t pointCount = track.getPointCount();
    vistrack->Reset(pointCount);

    for(size_t j = 0; j < pointCount; ++j) {
      auto point = track.getPoint(j);
      vistrack->SetNextPoint(point[0], point[1], point[2]);
    }
    list->AddElement(vistrack);
  }

  if(trackCount != 0)
  {
    dataTypeLists[EVisualisationDataType::ESD]->AddElement(list);
  }

  size_t clusterCount = event.getClusterCount();
  auto *point_list = new TEvePointSet(detectorName.c_str());
  point_list->IncDenyDestroy();
  point_list->SetMarkerColor(kBlue);

  for(size_t i = 0; i < clusterCount; ++i) {
    VisualisationCluster cluster = event.getCluster(i);
    point_list->SetNextPoint(cluster.X(), cluster.Y(), cluster.Z());
  }

  if(clusterCount != 0)
  {
    dataTypeLists[EVisualisationDataType::Clusters]->AddElement(point_list);
  }
}

void
EventManager::registerDetector(DataReader* reader, DataInterpreter* interpreter, EVisualisationGroup type)
{
  dataReaders[type] = reader;
  dataInterpreters[type] = interpreter;
}

}
}

