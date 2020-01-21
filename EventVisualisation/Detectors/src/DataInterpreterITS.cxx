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
/// \brief   Converting ITS data to Event Visualisation primitives
/// \author  julian.myrcha@cern.ch
/// \author  p.nowakowski@cern.ch
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#include "EventVisualisationDetectors/DataInterpreterITS.h"
#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationDataConverter/VisualisationEvent.h"
#include "EventVisualisationDataConverter/ConversionConstants.h"

#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/Cluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITSBase/GeometryTGeo.h"

#include <TEveManager.h>
#include <TEveTrackPropagator.h>
#include <TEveTrack.h>
#include <TGListTree.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector2.h>

#include <iostream>
#include <gsl/span>

namespace o2
{
namespace event_visualisation
{

DataInterpreterITS::DataInterpreterITS()
{
  //Prepare coordinate translator
  base::GeometryManager::loadGeometry("O2geometry.root", "FAIRGeom");
  its::GeometryTGeo* gman = its::GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::utils::bit2Mask(o2::TransformType::T2GRot));
}

void DataInterpreterITS::interpretDataForType(TObject* data, EVisualisationDataType type, VisualisationEvent& event)
{
  TList* list = (TList*)data;
  Int_t eventId = ((TVector2*)list->At(2))->X();

  if (type == Clusters) {
    its::GeometryTGeo* gman = its::GeometryTGeo::Instance();

    TFile* clustFile = (TFile*)list->At(1);
    TTree* clusters = (TTree*)clustFile->Get("o2sim");
    TTree* clustersRof = (TTree*)clustFile->Get("ITSClustersROF");

    //Read all clusters to a buffer
    std::vector<itsmft::Cluster>* clusArr = nullptr;
    clusters->SetBranchAddress("ITSCluster", &clusArr);
    clusters->GetEntry(0);

    //Read all cluster RO frames to a buffer
    std::vector<itsmft::ROFRecord>* clusterROFrames = nullptr;
    clustersRof->SetBranchAddress("ITSClustersROF", &clusterROFrames);
    clustersRof->GetEntry(0);
    auto currentClusterROF = clusterROFrames->at(eventId);

    int first, last;
    first = currentClusterROF.getROFEntry().getIndex();
    last = first + currentClusterROF.getNROFEntries();

    gsl::span<itsmft::Cluster> mClusters = gsl::make_span(&(*clusArr)[first], last - first);

    for (const auto& c : mClusters) {
      const auto& gloC = c.getXYZGloRot(*gman);
      double xyz[3] = {gloC.X(), gloC.Y(), gloC.Z()};
      VisualisationCluster cluster(xyz);
      event.addCluster(cluster);
    }
  } else if (type == Tracks) {
    TFile* trackFile = (TFile*)list->At(0);
    TFile* clustFile = (TFile*)list->At(1);

    TTree* tracksTree = (TTree*)trackFile->Get("o2sim");
    TTree* tracksRofTree = (TTree*)trackFile->Get("ITSTracksROF");

    TTree* clustersTree = (TTree*)clustFile->Get("o2sim");
    TTree* clustersRofTree = (TTree*)clustFile->Get("ITSClustersROF");

    //Read all tracks to a buffer
    std::vector<its::TrackITS>* trkArr = nullptr;
    tracksTree->SetBranchAddress("ITSTrack", &trkArr);
    tracksTree->GetEntry(0);

    //Read all track RO frames to a buffer
    std::vector<itsmft::ROFRecord>* trackROFrames = nullptr;
    tracksRofTree->SetBranchAddress("ITSTracksROF", &trackROFrames);
    tracksRofTree->GetEntry(0);

    //Read all clusters to a buffer
    std::vector<itsmft::Cluster>* clusArr = nullptr;
    clustersTree->SetBranchAddress("ITSCluster", &clusArr);
    clustersTree->GetEntry(0);

    //Read all cluster RO frames to a buffer
    std::vector<itsmft::ROFRecord>* clusterROFrames = nullptr;
    clustersRofTree->SetBranchAddress("ITSClustersROF", &clusterROFrames);
    clustersRofTree->GetEntry(0);

    TEveTrackList* trackList = new TEveTrackList("tracks");
    trackList->IncDenyDestroy();
    auto prop = trackList->GetPropagator();
    prop->SetMagField(0.5);
    prop->SetMaxR(50.);

    auto currentTrackROF = trackROFrames->at(eventId);

    int first, last;
    first = currentTrackROF.getROFEntry().getIndex();
    last = first + currentTrackROF.getNROFEntries();

    gsl::span<its::TrackITS> tracks = gsl::make_span(&(*trkArr)[first], last - first);

    int trkInd = 0;
    for (const auto& rec : tracks) {
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

      //                    TEvePointSet* tpoints = new TEvePointSet("tclusters");
      //                    int nc = rec.getNumberOfClusters();
      //                    while (nc--) {
      //                        Int_t idx = rec.getClusterEntry(nc);
      //                        itsmft::Cluster& c = (*clusArr)[idx];
      //                        const auto& gloC = c.getXYZGloRot(*gman);
      //                        tpoints->SetNextPoint(gloC.X(), gloC.Y(), gloC.Z());
      //                    }

      event.addTrack(track);
      trkInd++;
    }
    delete trackList;
  }
}

} // namespace event_visualisation
} // namespace o2
