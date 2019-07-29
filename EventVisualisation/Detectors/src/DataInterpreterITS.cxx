// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DataInterpreterITS.cxx
/// \brief converting ITS data to Event Visualisation primitives
/// \author julian.myrcha@cern.ch
/// \author p.nowakowski@cern.ch

#include "EventVisualisationDetectors/DataInterpreterITS.h"

#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationBase/Track.h"

#include "EventVisualisationDataConverter/MinimalisticEvent.h"

#include <TEveManager.h>
#include <TEveTrackPropagator.h>
#include <TGListTree.h>
#include <TFile.h>
#include <TTree.h>

#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/Cluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITSBase/GeometryTGeo.h"

#include <iostream>
#include <gsl/span>

using namespace std;

namespace o2 {
namespace event_visualisation {

DataInterpreterITS::DataInterpreterITS() = default;

DataInterpreterITS::~DataInterpreterITS() = default;

TEveElement* DataInterpreterITS::interpretDataForType(TObject* data, EDataType type) {
    TList *list = (TList*)data;
    Int_t event = ((TVector2*)list->At(2))->X();

    //Prepare coordinate translator
    base::GeometryManager::loadGeometry("O2geometry.root", "FAIRGeom");
    its::GeometryTGeo* gman = its::GeometryTGeo::Instance();
    gman->fillMatrixCache(o2::utils::bit2Mask(o2::TransformType::T2GRot));

    if(type == EDataType::Clusters) {
        TFile *clustFile = (TFile*)list->At(1);
        TTree *clusters = (TTree*)clustFile->Get("o2sim");
        TTree *clustersRof = (TTree*)clustFile->Get("ITSClustersROF");

        //Read all clusters to a buffer
        std::vector<itsmft::Cluster>* clusArr = nullptr;
        clusters->SetBranchAddress("ITSCluster", &clusArr);
        clusters->GetEntry(0);

        //Read all cluster RO frames to a buffer
        std::vector<itsmft::ROFRecord> *clusterROFrames = nullptr;
        clustersRof->SetBranchAddress("ITSClustersROF", &clusterROFrames);
        clustersRof->GetEntry(0);

        auto currentClusterROF = clusterROFrames->at(event);

        int first, last;
        first = currentClusterROF.getROFEntry().getIndex();
        last = first + currentClusterROF.getNROFEntries();

        gsl::span<itsmft::Cluster> mClusters = gsl::make_span(&(*clusArr)[first], last - first);

        TEvePointSet *cPointSet = new TEvePointSet("clusters");
        cPointSet->SetMarkerColor(kBlue);

        for (const auto& c : mClusters) {
            const auto& gloC = c.getXYZGloRot(*gman);
            cPointSet->SetNextPoint(gloC.X(), gloC.Y(), gloC.Z());
        }

        return cPointSet;
    }
    else if(type == EDataType::ESD) {
        TFile *trackFile = (TFile*)list->At(0);
        TFile *clustFile = (TFile*)list->At(1);

        TTree *tracks = (TTree*)trackFile->Get("o2sim");
        TTree *tracksRof = (TTree*)trackFile->Get("ITSTracksROF");

        TTree *clusters = (TTree*)clustFile->Get("o2sim");
        TTree *clustersRof = (TTree*)clustFile->Get("ITSClustersROF");

        //Read all tracks to a buffer
        std::vector<its::TrackITS>* trkArr = nullptr;
        tracks->SetBranchAddress("ITSTrack", &trkArr);
        tracks->GetEntry(0);

        //Read all track RO frames to a buffer
        std::vector<itsmft::ROFRecord> *trackROFrames = nullptr;
        tracksRof->SetBranchAddress("ITSTracksROF", &trackROFrames);
        tracksRof->GetEntry(0);

        //Read all clusters to a buffer
        std::vector<itsmft::Cluster>* clusArr = nullptr;
        clusters->SetBranchAddress("ITSCluster", &clusArr);
        clusters->GetEntry(0);

        //Read all cluster RO frames to a buffer
        std::vector<itsmft::ROFRecord> *clusterROFrames = nullptr;
        clustersRof->SetBranchAddress("ITSClustersROF", &clusterROFrames);
        clustersRof->GetEntry(0);

        TEveTrackList* trackList = new TEveTrackList("tracks");
        trackList->IncDenyDestroy();
        auto prop = trackList->GetPropagator();
        prop->SetMagField(0.5);
        prop->SetMaxR(50.);

        auto currentTrackROF = trackROFrames->at(event);

        int first, last;
        first = currentTrackROF.getROFEntry().getIndex();
        last = first + currentTrackROF.getNROFEntries();

        gsl::span<o2::its::TrackITS> mTracks = gsl::make_span(&(*trkArr)[first], last - first);

        struct TrackletTest {
            std::array<float, 3> mP;
            Int_t sign;
            std::vector<double> mPolyX;
            std::vector<double> mPolyY;
            std::vector<double> mPolyZ;
        };

        std::vector<TrackletTest> points;

        for (const auto& rec : mTracks) {
            std::array<float, 3> p;
            rec.getPxPyPzGlo(p);
            TEveRecTrackD t;
            t.fP = { p[0], p[1], p[2] };
            t.fSign = (rec.getSign() < 0) ? -1 : 1;
            TEveTrack* track = new TEveTrack(&t, prop);
            track->MakeTrack();

            TrackletTest test;

            test.mP = { p[0], p[1], p[2] };
            test.sign = t.fSign;

//            TEvePointSet* tpoints = new TEvePointSet("tclusters");
//            int nc = rec.getNumberOfClusters();
//            while (nc--) {
//                Int_t idx = rec.getClusterEntry(nc);
//                itsmft::Cluster& c = (*clusArr)[idx];
//                const auto& gloC = c.getXYZGloRot(*gman);
//                tpoints->SetNextPoint(gloC.X(), gloC.Y(), gloC.Z());
//            }

            for(Int_t i = 0; i < track->GetN(); ++i) {
                Float_t x,y,z;
                track->GetPoint(i, x, y, z);
                test.mPolyX.push_back(x);
                test.mPolyY.push_back(y);
                test.mPolyZ.push_back(z);
            }

            points.push_back(test);
        }

        for(const auto& track: points) {
            TEveRecTrackD t;
            t.fP = { track.mP[0], track.mP[1], track.mP[2] };
            t.fSign = track.sign;
            TEveTrack* vistrack = new TEveTrack(&t, &TEveTrackPropagator::fgDefault);
            vistrack->SetLineColor(kMagenta);
            vistrack->Reset(track.mPolyZ.size());

            for(size_t i = 0; i < track.mPolyZ.size(); ++i) {
                vistrack->SetNextPoint(track.mPolyX[i], track.mPolyY[i], track.mPolyZ[i]);
            }
            trackList->AddElement(vistrack);
        }

        return trackList;
    }
}

}
}
