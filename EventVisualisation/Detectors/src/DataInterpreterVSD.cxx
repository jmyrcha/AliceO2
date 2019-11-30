// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DataInterpreterVSD.cxx
/// \brief Converting VSD data to Event Visualisation primitives
/// \author julian.myrcha@cern.ch
/// \author p.nowakowski@cern.ch

#include "EventVisualisationDetectors/DataInterpreterVSD.h"
#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationDataConverter/VisualisationEvent.h"

#include <TEveManager.h>
#include <TEveTrackPropagator.h>
#include <TGListTree.h>

#include <iostream>

using namespace std;

namespace o2
{
namespace event_visualisation
{

DataInterpreterVSD::~DataInterpreterVSD()
{
  //this->DropEvent();
  if (mVSD) {
    delete mVSD;
    mVSD = nullptr;
  }
}

void DataInterpreterVSD::interpretDataForType(TObject* data, EVisualisationDataType type, VisualisationEvent& event)
{
  if (mVSD == nullptr)
    mVSD = new TEveVSD;
  this->DropEvent();

  // Connect to new event-data.
  this->mDirectory = dynamic_cast<TDirectory*>(data);
  this->mVSD->SetDirectory(this->mDirectory);

  this->AttachEvent();

  // Load event data into visualization structures.

  //        this->LoadClusters(this->fITSClusters, "ITS", 0);
  //        this->LoadClusters(this->fTPCClusters, "TPC", 1);
  //        this->LoadClusters(this->fTRDClusters, "TRD", 2);
  //        this->LoadClusters(this->fTOFClusters, "TOF", 3);
  if (type == Tracks) {
    LoadEsdTracks(event);
  }
}

void DataInterpreterVSD::LoadClusters(TEvePointSet*& ps, const TString& det_name, Int_t det_id)
{
  if (ps == nullptr) {
    ps = new TEvePointSet(det_name);
    ps->SetMainColor((Color_t)(det_id + 2));
    ps->SetMarkerSize(0.5);
    ps->SetMarkerStyle(2);
    ps->IncDenyDestroy();
  } else {
    ps->Reset();
  }

  //TEvePointSelector ss(fVSD->fTreeC, ps, "fV.fX:fV.fY:fV.fZ", TString::Format("fDetId==%d", det_id));
  //ss.Select();
  ps->SetTitle(TString::Format("N=%d", ps->Size()));

  gEve->AddElement(ps);
}

void DataInterpreterVSD::AttachEvent()
{
  // Attach event data from current directory.
  mVSD->LoadTrees();
  mVSD->SetBranchAddresses();
}

void DataInterpreterVSD::DropEvent()
{
  assert(mVSD != nullptr);
  // Drop currently held event data, release current directory.
  // Drop old visualization structures.
  this->mViewers = gEve->GetViewers();
  this->mViewers->DeleteAnnotations();
  //TEveEventManager *manager = gEve->GetCurrentEvent();
  //assert(manager != nullptr);
  //manager->DestroyElements();

  // Drop old event data.
  mVSD->DeleteTrees();
  delete mDirectory;
  mDirectory = nullptr;
}

void DataInterpreterVSD::LoadEsdTracks(VisualisationEvent& event)
{
  TEveTrackList* list = new TEveTrackList();
  TEveTrackPropagator* trkProp = list->GetPropagator();
  trkProp->SetMagField(0.5);
  trkProp->SetStepper(TEveTrackPropagator::kRungeKutta);

  Int_t nTracks = mVSD->fTreeR->GetEntries();
  for (Int_t n = 0; n < nTracks; n++) {
    mVSD->fTreeR->GetEntry(n);

    auto* eve_track = new TEveTrack(&mVSD->fR, trkProp);
    eve_track->MakeTrack();

    auto start = eve_track->GetLineStart();
    auto end = eve_track->GetLineEnd();
    auto p = eve_track->GetMomentum();

    double track_start[3] = { start.fX, start.fY, start.fZ };
    double track_end[3] = { end.fX, end.fY, end.fZ };
    double track_p[3] = { p[0], p[1], p[2] };

    VisualisationTrack track(eve_track->GetCharge(), 0.0, 0, 0, 0.0, 0.0, track_start, track_end, track_p, 0, 0.0, 0.0, 0.0, 0, 0);

    for (Int_t i = 0; i < eve_track->GetN(); ++i) {
      Float_t x, y, z;
      eve_track->GetPoint(i, x, y, z);
      track.addPolyPoint(x, y, z);
    }
    delete eve_track;

    event.addTrack(track);
  }
  delete list;
}

} // namespace event_visualisation
} // namespace o2
