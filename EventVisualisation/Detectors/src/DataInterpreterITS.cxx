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
/// \file    DataInterpreterITS.cxx
/// \author  Julian Myrcha

#include "EventVisualisationDetectors/DataInterpreterITS.h"

#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationBase/Track.h"
#include "EventVisualisationDataConverter/MinimalisticEvent.h"
#include "EventVisualisationDetectors/DataSourceOfflineITS.h"

#include <TEveManager.h>
#include <TEveTrackPropagator.h>
#include <TGListTree.h>

#include <iostream>

using namespace std;

namespace o2 {
namespace event_visualisation {

DataInterpreterITS::DataInterpreterITS() = default;

DataInterpreterITS::~DataInterpreterITS() = default;

TEveElement* DataInterpreterITS::interpretDataForType(EDataType type) {
  // Assuming offline - where online vs offline should be recognized?
  DataSourceOfflineITS dataSource;

  switch(type) {
    case EDataType::Raw:
      dataSource.OpenRawFile();
      break;
    case EDataType::Hits:
      break;
    case EDataType::Digits:
      break;
    case EDataType::Clusters:
      dataSource.OpenClustersFile();
      break;
    case EDataType::ESD:
      break;
    case EDataType::AOD:
      dataSource.OpenClustersFile();
      dataSource.OpenTracksFile();
      break;
    case EDataType::NoData:
    default:
      break;
  }

  if(!dataSource.GotoEvent(0))
  {
    std::cout << "Data source could not go to event 0" << std::endl;
    return nullptr;
  }

  TEveElement* event;
  switch(type) {
    case EDataType::Raw:
      break;
    case EDataType::Hits:
      break;
    case EDataType::Digits:
      break;
    case EDataType::Clusters:
      break;
    case EDataType::ESD:
      break;
    case EDataType::AOD:
      event = loadTracks();
      break;
    case EDataType::NoData:
    default:
      break;
  }

  return event;
}

TEveElement* DataInterpreterITS::loadTracks()
{
  auto gman = o2::its::GeometryTGeo::Instance();

  TEvePointSet* clusters = new TEvePointSet("clusters");
  clusters->SetMarkerColor(kBlue);
  for (const auto& c : mClusters) {
    const auto& gloC = c.getXYZGloRot(*gman);
    clusters->SetNextPoint(gloC.X(), gloC.Y(), gloC.Z());
  }
  std::cout << "Got eve clusters" << '\n';

  TEveTrackList* tracks = new TEveTrackList("tracks");
  auto prop = tracks->GetPropagator();
  prop->SetMagField(0.5);
  prop->SetMaxR(50.);
  for (const auto& rec : mTracks) {
    std::array<float, 3> p;
    rec.getPxPyPzGlo(p);
    TEveRecTrackD t;
    t.fP = { p[0], p[1], p[2] };
    t.fSign = (rec.getSign() < 0) ? -1 : 1;
    TEveTrack* track = new TEveTrack(&t, prop);
    track->SetLineColor(kMagenta);
    tracks->AddElement(track);

    if (mClusters.empty())
      continue;
    TEvePointSet* tpoints = new TEvePointSet("tclusters");
    tpoints->SetMarkerColor(kGreen);
    int nc = rec.getNumberOfClusters();
    int idxRef = rec.getFirstClusterEntry();
    while (nc--) {
      Int_t idx = (*mClIdxBuffer)[idxRef + nc];
      const Cluster& c = mClusters[idx];
      const auto& gloC = c.getXYZGloRot(*gman);
      tpoints->SetNextPoint(gloC.X(), gloC.Y(), gloC.Z());
    }
    track->AddElement(tpoints);
  }
  tracks->MakeTracks();
  std::cout << "Got eve tracks" << '\n';

  TEveElement* event = new TEveElement();
  event->AddElement(clusters);
  event->AddElement(tracks);
  std::cout << "Added elements" << '\n';
  return event;
}

}
}
