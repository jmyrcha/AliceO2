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
/// \file    DataInterpreterTPC.cxx
/// \brief   Converting TPC data to Event Visualisation primitives
/// \author  julian.myrcha@cern.ch
/// \author  p.nowakowski@cern.ch
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#include "EventVisualisationDetectors/DataInterpreterTPC.h"
#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationDataConverter/VisualisationEvent.h"

#include "TPCBase/Mapper.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "DataFormatsTPC/ClusterNative.h"
#include "DataFormatsTPC/ClusterNativeHelper.h"

#include <TEveManager.h>
#include <TEveTrackPropagator.h>
#include <TEveTrack.h>
#include <TGListTree.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector2.h>
#include <TROOT.h>

#include <iostream>
#include <gsl/span>

namespace o2
{
namespace event_visualisation
{

void DataInterpreterTPC::interpretTracks(TFile* file, Int_t eventId, VisualisationEvent& event)
{
  TTree* tracks = (TTree*)file->Get("tpcrec");

  // Read all tracks to a buffer
  std::vector<tpc::TrackTPC>* trkArr = nullptr;
  tracks->SetBranchAddress("TPCTracks", &trkArr);
  tracks->GetEntry(0);

  TEveTrackList* trackList = new TEveTrackList("tracks");
  trackList->IncDenyDestroy();
  auto prop = trackList->GetPropagator();
  prop->SetMagField(0.5); // 0.1 * field in simulation workflow (different units)
  prop->SetMaxR(250.0);

  // Tracks are not in order
  int startTime = 2 * mTPCReadoutCycle * eventId;
  int endTime = startTime + 2 * mTPCReadoutCycle;

  int trkInd = 0;
  for (int i = 0; i < trkArr->size(); i++) {
    const auto& rec = (*trkArr)[i];

    if (rec.getTime0() < startTime || rec.getTime0() > endTime) {
      continue;
    }

    std::array<float, 3> p;
    rec.getPxPyPzGlo(p);

    TEveRecTrackD t;
    t.fP = {p[0], p[1], p[2]};
    t.fSign = (rec.getSign() < 0) ? -1 : 1;
    TEveTrack* eve_track = new TEveTrack(&t, prop);
    eve_track->MakeTrack();

    auto start = eve_track->GetLineStart();
    auto end = eve_track->GetLineEnd();
    double track_start[3] = {start.fX, start.fY, start.fZ};
    double track_end[3] = {end.fX, end.fY, end.fZ};
    double track_p[3] = {p[0], p[1], p[2]};
    o2::track::PID pid;
    double mass = pid.getMass();
    double momentum = rec.getP();
    double energy = TMath::Sqrt(momentum * momentum + mass * mass);
    // FIXME: Temporarily harcoded magnetic field
    float bz = 5.0f;

    // int ID, int type, int charge, double energy, int parentID, o2::track::PID PID, double signedPT, double mass, double pxpypz[], double startXYZ[], double endXYZ[], double helixCurvature, double theta, double phi, float C1Pt21Pt2, unsigned long long flags
    VisualisationTrack track(trkInd, ETrackType::Standard, eve_track->GetCharge(), energy, -1, pid, 1.0 / rec.getQ2Pt(), mass, track_p, track_start, track_end, rec.getCurvature(bz), rec.getTheta(), rec.getPhi(), rec.getSigma1Pt2(), 0);

    for (int j = 0; j < eve_track->GetN(); j++) {
      float x, y, z;
      eve_track->GetPoint(j, x, y, z);
      track.addPolyPoint(x, y, z);
    }
    delete eve_track;

    event.addTrack(track);
    trkInd++;
  }
  delete trackList;
}

void DataInterpreterTPC::interpretClusters(TFile* file, Int_t eventId, VisualisationEvent& event)
{
  //Why cannot TPC clusters be read like other clusters?
  auto reader = new tpc::ClusterNativeHelper::Reader();
  reader->init(file->GetName());
  reader->read(0);

  auto access = std::make_unique<o2::tpc::ClusterNativeAccess>();
  std::unique_ptr<tpc::ClusterNative[]> clusterBuffer;
  tpc::MCLabelContainer clusterMCBuffer;

  reader->fillIndex(*access, clusterBuffer, clusterMCBuffer);

  const auto& mapper = tpc::Mapper::instance();
  const auto& clusterRefs = access->clusters;

  for (int sector = 0; sector < o2::tpc::Constants::MAXSECTOR; sector++) {
    for (int row = 0; row < o2::tpc::Constants::MAXGLOBALPADROW; row++) {
      const auto& c = clusterRefs[sector][row];

      const auto pad = mapper.globalPadNumber(tpc::PadPos(row, c->getPad()));
      const tpc::LocalPosition3D localXYZ(mapper.padCentre(pad).X(), mapper.padCentre(pad).Y(), c->getTime());
      const auto globalXYZ = mapper.LocalToGlobal(localXYZ, sector);
      double xyz[3] = {globalXYZ.X(), globalXYZ.Y(), globalXYZ.Z()};

      VisualisationCluster cluster(xyz);
      event.addCluster(cluster);
    }
  }
}

} // namespace event_visualisation
} // namespace o2
