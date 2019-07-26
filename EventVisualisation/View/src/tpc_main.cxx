//
// Created by jmy on 10.07.19.
//

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
/// \file    main.cxx
/// \author  Jeremi Niedziela
///

#include "EventVisualisationView/Initializer.h"
#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationBase/VisualisationConstants.h"

#include <TApplication.h>
#include <TEveBrowser.h>
#include <TEveManager.h>

#include <ctime>
#include <iostream>

using namespace std;
using namespace o2::event_visualisation;

#include <iostream>
#include <array>
#include <algorithm>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TEveManager.h>
#include <TEveBrowser.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGFrame.h>
#include <TGTab.h>
#include <TGLCameraOverlay.h>
#include <TEveFrameBox.h>
#include <TEveQuadSet.h>
#include <TEveTrans.h>
#include <TEvePointSet.h>
#include <TEveTrackPropagator.h>
#include <TEveTrack.h>
#include <Rtypes.h>
#include <gsl/span>
#include <DetectorsBase/GeometryManager.h>

#include "EventVisualisationView/MultiView.h"

//#include "ITSMFTReconstruction/ChipMappingITS.h"
//#include "ITSMFTReconstruction/DigitPixelReader.h"
//#include "ITSMFTReconstruction/RawPixelReader.h"
//#include "ITSMFTBase/SegmentationAlpide.h"
//#include "TPCBase/Digit.h"
#include "TPCBase/Mapper.h"
//#include "ITSBase/GeometryTGeo.h"
//#include "DataFormatsTPC/ClusterNative.h"
//#include "DataFormatsTPC/ClusterNativeHelper.h"
//#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsTPC/TrackTPC.h"
//#include "DetectorsCommonDataFormats/DetID.h"
//#include "CommonDataFormat/InteractionRecord.h"
#include "TGenericClassInfo.h"

using namespace o2::tpc;

extern TEveManager* gEve;

static TGNumberEntry* gEntry;

class Data
{
public:
    void loadData(int entry);
    void displayData(int entry);
    int getLastEvent() const { return mLastEvent; }
    //void setDigiTree(TTree* t) { mDigiTree = t; }
    //void setClusTree(TTree* t);
    void setTracTree(TTree* t);

private:
    // Data loading members
    //ClusterNativeHelper::Reader reader;
    int mLastEvent = 0;
    //std::vector<Digit> mDigits;
    //ClusterNativeAccess clusterIndex;
    //std::unique_ptr<ClusterNative[]> mClusterBuffer;
    //MCLabelContainer mcBuffer;
    //std::vector<ClusterNative>* mClusterBuffer = nullptr;
    //gsl::span<ClusterNative> mClusters;
    //std::unique_ptr<ClusterNativeAccess> mClusters;
    //vector<TrackTPC> tracks;
    //MCLabelContainer tracksMC;
    std::vector<TrackTPC>* mTrackBuffer = nullptr;
    gsl::span<TrackTPC> mTracks;
    //void loadDigits();
    //void loadDigits(int entry);
    //void loadClusters(int entry);
    void loadTracks(int entry);

    //TTree* mDigiTree = nullptr;
    //TTree* mClusTree = nullptr;
    TTree* mTracTree = nullptr;

    // TEve-related members
    TEveElementList* mEvent = nullptr;
    //TEveElement* getEveClusters();
    TEveElement* getEveTracks();
} evdata;

//void Data::setClusTree(TTree* tree)
//{
//    if (tree == nullptr) {
//        std::cerr << "No tree for clusters !\n";
//        return;
//    }
//    tree->SetBranchAddress("TPCCluster", &mClusterBuffer);
//    mClusTree = tree;
//}

//void Data::loadClusters(int entry)
//{
//    static int lastLoaded = -1;
//
//    if (mClusTree == nullptr)
//        return;
//
//    auto event = entry;
//    if ((event < 0) || (event >= mClusTree->GetEntries())) {
//        std::cerr << "Clusters: Out of event range ! " << event << '\n';
//        return;
//    }
//    if (event != lastLoaded) {
//        mClusterBuffer->clear();
//        mClusTree->GetEntry(event);
//        lastLoaded = event;
//    }
//
//    int first = 0, last = mClusterBuffer->size();
//    mClusters = gsl::make_span(&(*mClusterBuffer)[first], last - first);
//
//    std::cout << "Number of TPC Clusters: " << mClusters.size() << '\n';
//}

void Data::setTracTree(TTree* tree)
{
    if (tree == nullptr) {
        std::cerr << "No tree for tracks !\n";
        return;
    }
    tree->SetBranchAddress("TPCTracks", &mTrackBuffer);
    mTracTree = tree;
}

void Data::loadTracks(int entry)
{
    static int lastLoaded = -1;

    if (mTracTree == nullptr)
        return;

    auto event = entry;
    if ((event < 0) || (event >= mTracTree->GetEntries())) {
        std::cerr << "Tracks: Out of event range ! " << event << '\n';
        return;
    }
    if (event != lastLoaded) {
        mTrackBuffer->clear();
        mTracTree->GetEntry(event);
        lastLoaded = event;
    }

    int first = 0, last = mTrackBuffer->size();
    mTracks = gsl::make_span(&(*mTrackBuffer)[first], last - first);

    std::cout << "Number of TPC Tracks: " << mTracks.size() << '\n';
}

void Data::loadData(int entry)
{
    //loadDigits(entry);
    //loadClusters(entry);
    loadTracks(entry);
}

//TEveElement* Data::getEveClusters()
//{
//    const auto& mapper = Mapper::instance();
//    TEvePointSet* clusters = new TEvePointSet("clusters");
//    clusters->SetMarkerColor(kBlue);
//    for (const auto& c : mClusters) {
//        const auto pad = mapper.globalPadNumber(PadPos(c.getRow(),
//                                                       c.getPad()));
//        const auto& localXYZ = mapper.padCentre(pad);
//        const auto globalXYZ = mapper.LocalToGlobal(localXYZ,
//                                                    CRU(c.getCRU()).sector());
//        clusters->SetNextPoint(globalXYZ.X(), globalXYZ.Y(), globalXYZ.Z());
//    }
//    return clusters;
//}

TEveElement* Data::getEveTracks()
{
    TEveTrackList* tracks = new TEveTrackList("tracks");
    auto prop = tracks->GetPropagator();
    prop->SetMagField(0.5);
    //prop->SetMaxR(50.);
    for (const auto& rec : mTracks) {
        std::array<float, 3> p;
        rec.getPxPyPzGlo(p);
        TEveRecTrackD t;
        t.fP = { p[0], p[1], p[2] };
        t.fSign = (rec.getSign() < 0) ? -1 : 1;
        TEveTrack* track = new TEveTrack(&t, prop);
        track->SetLineColor(kMagenta);
        tracks->AddElement(track);

//        TEvePointSet* tpoints = new TEvePointSet("tclusters");
//        tpoints->SetMarkerColor(kGreen);
//        int nc = rec.getNClusterReferences();
//        while (nc--) {
//            uint8_t sector, row;
//            uint32_t clusterIndexInRow;
//            rec.getClusterReference(nc, sector, row, clusterIndexInRow);
//            const ClusterNative& cl = rec.getCluster(nc, *clusters, sector, row);
//            const ClusterNative& clLast = rec.getCluster(0, *clusters);
//            const auto& gloC = c.getXYZGloRot(*gman);
//            tpoints->SetNextPoint(gloC.X(), gloC.Y(), gloC.Z());
//        }
//        track->AddElement(tpoints);
    }
    tracks->MakeTracks();

    return tracks;
}

void Data::displayData(int entry)
{
    std::string ename("Event #");
    ename += std::to_string(entry);

    // Event display
    //auto clusters = getEveClusters();
    auto tracks = getEveTracks();
    delete mEvent;
    mEvent = new TEveElementList(ename.c_str());
    //mEvent->AddElement(clusters);
    mEvent->AddElement(tracks);
    auto multi = o2::event_visualisation::MultiView::getInstance();
    multi->registerElement(mEvent);

    gEve->Redraw3D(kFALSE);
}

void load(int entry)
{
    int lastEvent = evdata.getLastEvent();
    if (lastEvent > entry) {
        std::cerr << "\nERROR: Cannot stay or go back over events. Please increase the event number !\n\n";
        gEntry->SetIntNumber(lastEvent - 1);
        return;
    }

    gEntry->SetIntNumber(entry);

    std::cout << "\n*** Event #" << entry << " ***\n";
    evdata.loadData(entry);
    evdata.displayData(entry);
}


void load()
{
    auto event = gEntry->GetNumberEntry()->GetIntNumber();
    load(event);
}

void next()
{
    auto event = gEntry->GetNumberEntry()->GetIntNumber();
    event++;
    load(event);
}

void prev()
{
    auto event = gEntry->GetNumberEntry()->GetIntNumber();
    event--;
    load(event);
}

void setupGeometry()
{
  // read path to geometry files from config file
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);

  // get geometry from Geometry Manager and register in multiview
  auto multiView = MultiView::getInstance();

  for(int iDet=0;iDet<NvisualisationGroups;++iDet){
    EVisualisationGroup det = static_cast<EVisualisationGroup>(iDet);
    std::string detName = gVisualisationGroupName[det];
    if(settings.GetValue((detName+".draw").c_str(), false))
    {
      if(   detName=="TPC" || detName=="MCH" || detName=="MTR"
            || detName=="MID" || detName=="MFT" || detName=="AD0"
            || detName=="FMD"){// don't load MUON+MFT and AD and standard TPC to R-Phi view

        multiView->drawGeometryForDetector(detName, false, true, false);
      }
      else if(detName=="RPH"){// special TPC geom from R-Phi view

        multiView->drawGeometryForDetector(detName, false, false, true, false);
      }
      else{// default
        multiView->drawGeometryForDetector(detName);
      }
    }
  }
}

int main(int argc, char **argv)
{
    // create ROOT application environment
    TApplication *app = new TApplication("o2-tpc-eve", &argc, argv);
    app->Connect("TEveBrowser", "CloseWindow()", "TApplication", app, "Terminate()");

    cout<<"Initializing TEveManager"<<endl;
    if(!(gEve=TEveManager::Create())){
        cout<<"FATAL -- Could not create TEveManager!!"<<endl;
        exit(0);
    }

    int entry = 0;
    //std::string digifile = "tpcdigits.root";
    //bool rawdata = false;
    std::string tracfile = "tpctracks.root";
    std::string inputGeom = "O2geometry.root";

    // Geometry
    o2::base::GeometryManager::loadGeometry(inputGeom, "FAIRGeom");
    TEveBrowser* browser = gEve->GetBrowser();

    // Event View
    std::cout << "Going to setup the geometry..." << std::endl;
    setupGeometry();

    // Event navigation
    browser->StartEmbedding(TRootBrowser::kBottom);
    auto frame = new TGMainFrame(gClient->GetRoot(), 1000, 600, kVerticalFrame);

    auto h = new TGHorizontalFrame(frame);
    auto b = new TGTextButton(h, "PrevEvnt", "prev()");
    h->AddFrame(b);
    gEntry = new TGNumberEntry(h, 0, 5, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 0, 10000);
    gEntry->Connect("ValueSet(Long_t)", 0, 0, "load()");
    h->AddFrame(gEntry);
    b = new TGTextButton(h, "NextEvnt", "next()");
    h->AddFrame(b);

    frame->AddFrame(h);

    frame->MapSubwindows();
    frame->MapWindow();
    browser->StopEmbedding("Navigator");

    TFile* file;

    // Data sources
//    if (rawdata) {
//        std::ifstream* rawfile = new std::ifstream(digifile.data(), std::ifstream::binary);
//        if (rawfile->good()) {
//            delete rawfile;
//            std::cout << "Running with raw digits...\n";
//            evdata.setRawPixelReader(digifile.data());
//        } else
//            std::cerr << "\nERROR: Cannot open file: " << digifile << "\n\n";
//    } else {
//        file = TFile::Open(digifile.data());
//        if (file && gFile->IsOpen()) {
//            file->Close();
//            std::cout << "Running with MC digits...\n";
//            evdata.setDigitPixelReader(digifile.data());
//            //evdata.setDigiTree((TTree*)gFile->Get("o2sim"));
//        } else
//            std::cerr << "\nERROR: Cannot open file: " << digifile << "\n\n";
//    }

    file = TFile::Open(tracfile.data());
    if (file && gFile->IsOpen()) {
      evdata.setTracTree((TTree*)gFile->Get("tpcrec"));
    }
    else
        std::cerr << "\nERROR: Cannot open file: " << tracfile << "\n\n";

//    std::cout << "\n **** Navigation over events ****\n";
//    std::cout << " load(event) \t jump to the specified event \n";
//    std::cout << " next() \t\t load next event \n";
//    std::cout << " prev() \t\t load previous event \n";

    // Manually adding an event
    gEve->AddEvent(new TEveEventManager("Event", "ALICE TPC Event"));

    load(entry);

    // Start the application
    app->Run(kTRUE);

    // Terminate application
    TEveManager::Terminate();
    app->Terminate();

    return 0;
}

void tpc_main()
{
  // A dummy function with the same name as this macro
}