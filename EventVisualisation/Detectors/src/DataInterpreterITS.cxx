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

#include <TEveManager.h>
#include <TEveTrackPropagator.h>
#include <TGListTree.h>
#include <TFile.h>
#include <TTree.h>

#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/Cluster.h"
#include "ITSBase/GeometryTGeo.h"

#include <iostream>

using namespace std;

namespace o2 {
namespace event_visualisation {

DataInterpreterITS::DataInterpreterITS() = default;

DataInterpreterITS::~DataInterpreterITS() = default;

TEveElement* DataInterpreterITS::interpretDataForType(TObject* data, EDataType type) {
    auto *file = (TFile *) data;
    TTree* tree = (TTree*)file->Get("o2sim");

    std::cerr << tree->GetName() << std::endl;

    std::vector<its::TrackITS>* trkArr = nullptr;
    tree->SetBranchAddress("ITSTrack", &trkArr);
    tree->GetEntry(0);

    auto *ff = TFile::Open("o2clus_its.root");

    TTree* tree2 = (TTree*)ff->Get("o2sim");
    std::cerr << tree2->GetName() << std::endl;

    std::vector<itsmft::Cluster>* clusArr = nullptr;
    tree2->SetBranchAddress("ITSCluster", &clusArr);

    std::cerr << trkArr << " " << clusArr << std::endl;

    TEveTrackList* tracks = new TEveTrackList("tracks");
    auto prop = tracks->GetPropagator();
    prop->SetMagField(0.5);
    prop->SetMaxR(50.);

    auto *points = new TEvePointSet("event_track");

    o2::base::GeometryManager::loadGeometry("O2geometry.root", "FAIRGeom");
    o2::its::GeometryTGeo* gman = its::GeometryTGeo::Instance();
    gman->fillMatrixCache(o2::utils::bit2Mask(o2::TransformType::T2GRot));

    for (const auto& rec : *trkArr) {
        std::array<float, 3> p;
        rec.getPxPyPzGlo(p);
        TEveRecTrackD t;
        t.fP = { p[0], p[1], p[2] };
        t.fSign = (rec.getSign() < 0) ? -1 : 1;
        TEveTrack* track = new TEveTrack(&t, prop);
        track->SetLineColor(kMagenta);
        tracks->AddElement(track);
        TEvePointSet* tpoints = new TEvePointSet("tclusters");
        tpoints->SetMarkerColor(kGreen);
        int nc = rec.getNumberOfClusters();
        while (nc--) {
            Int_t idx = rec.getClusterEntry(nc); //Is this the correct function??
            itsmft::Cluster& c = (*clusArr)[idx];
            const auto& gloC = c.getXYZGloRot(*gman);
            tpoints->SetNextPoint(gloC.X(), gloC.Y(), gloC.Z());
        }
        track->AddElement(tpoints);
    }
    tracks->MakeTracks();

    return tracks;
}

}
}
