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

namespace o2  {
namespace event_visualisation {

EventManager *EventManager::instance = nullptr;

EventManager& EventManager::getInstance()
{
  if( instance == nullptr)
      instance = new EventManager();
  return *instance;
}

EventManager::EventManager() : TEveEventManager("Event", "") {
}

void EventManager::Open() {
    switch(mCurrentDataSourceType)
    {
        case SourceOnline:
            break;
        case SourceOffline: {
              DataSourceOffline *source = new DataSourceOffline();
              if(DataInterpreter::getInstance(EVisualisationGroup::VSD)) {
                  DataReader *vsd = new DataReaderVSD();
                  vsd->open();
                  source->registerReader(vsd, EVisualisationGroup::VSD);
              }
              if(DataInterpreter::getInstance(EVisualisationGroup::ITS)) {
                  DataReader *its = new DataReaderITS();
                  its->open();
                  source->registerReader(its, EVisualisationGroup::ITS);
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
    for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; i++) {
      DataInterpreter* interpreter = DataInterpreter::getInstance((EVisualisationGroup)i);
      if(interpreter) {
        TObject *data = getDataSource()->getEventData(no, (EVisualisationGroup)i);
        std::unique_ptr<VisualisationEvent> event = interpreter->interpretDataForType(data, ESD);
        displayVisualisationEvent(*event);
      }
    }
}

void EventManager::NextEvent() {
    Int_t event = (this->currentEvent + 1) % getDataSource()->GetEventCount();
    GotoEvent(event);
}

void EventManager::PrevEvent() {
    GotoEvent(this->currentEvent - 1);
}

void EventManager::Close() {
    delete this->dataSource;
    this->dataSource = nullptr;
}

void EventManager::AfterNewEventLoaded() {
    TEveEventManager::AfterNewEventLoaded();
}

void EventManager::AddNewEventCommand(const TString &cmd) {
    TEveEventManager::AddNewEventCommand(cmd);
}

void EventManager::RemoveNewEventCommand(const TString &cmd) {
    TEveEventManager::RemoveNewEventCommand(cmd);
}

void EventManager::ClearNewEventCommands() {
    TEveEventManager::ClearNewEventCommands();
}

EventManager::~EventManager() {
    instance = nullptr;
}

void EventManager::DropEvent() {
  DestroyElements();
}

void EventManager::displayVisualisationEvent(VisualisationEvent &event) {
    size_t trackCount = event.getTrackCount();

    auto *list = new TEveTrackList;
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
    MultiView::getInstance()->registerElement(list);
}

}
}

