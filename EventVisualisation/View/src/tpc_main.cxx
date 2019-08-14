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

#include "TPCSimulation/Point.h" // for o2::tpc::HitGroup
#include "TPCBase/Digit.h"
#include "TPCBase/Mapper.h"
#include "TPCReconstruction/RawReader.h"
#include "DataFormatsTPC/ClusterNative.h"
#include "DataFormatsTPC/ClusterNativeHelper.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "TPCSimulation/SAMPAProcessing.h"
#include "TGenericClassInfo.h"
#include "TPCBase/CDBInterface.h"

using namespace o2::tpc;

extern TEveManager* gEve;

static TGNumberEntry* gEntry;

// TODO: Any TPC official constant?
const int tpcReadoutCycle = 100; // in ms

class Data
{
public:
    void loadData(int entry);
    void displayData(int entry);
    int getLastEvent() const { return mLastEvent; }
    void setHitsTree(TTree* t);
    void setDigiTree(TTree* t);
    //void setClusTree(TTree* t);
    void setTracTree(TTree* t);

    void setRawReader(std::string infile)
    {
      auto reader = new RawReader();
      if(!reader->addInputFile(-1, -1, -1, infile, -1)) {
        std::cout << "Could not add file for raw reader!" << std::endl;
        return;
      }
      mRawReader = reader;
      mRawReader->loadEvent(0);
    }
    void setClusReader(std::string infile) {
      auto reader = new ClusterNativeHelper::Reader();
      reader->init(infile.c_str());
      mClusterIndex = std::make_unique<o2::tpc::ClusterNativeAccess>();
      mClusReader = reader;
      std::cout << "Cluster reader set" << std::endl;
    }
    bool getRawData() { return mRawData; }
    void setRawData(bool rawData) { mRawData = rawData; }

private:
    // Data loading members
    ClusterNativeHelper::Reader* mClusReader = nullptr;
    bool mRawData = false;
    RawReader* mRawReader = nullptr;
    int mLastEvent = 0;
    int mEventsCount = 0;

    std::vector<std::vector<HitGroup>*> mHitsBuffer;
    std::vector<gsl::span<HitGroup>> mHits;
    std::vector<std::vector<Digit>*> mDigitBuffer;
    std::vector<gsl::span<Digit>> mDigits;
    std::unique_ptr<ClusterNative[]> mClusterBuffer;
    MCLabelContainer mClusterMCBuffer;
    std::unique_ptr<ClusterNativeAccess> mClusterIndex;
    const ClusterNativeAccess* mClusterIndexStruct = nullptr;
    //std::vector<ClusterNative>* mClusterBuffer = nullptr;
    //gsl::span<ClusterNative> mClusters;
    std::vector<TrackTPC>* mTrackBuffer = nullptr;
    gsl::span<TrackTPC> mTracks;
    void loadHits(int entry);
    void loadDigits();
    void loadDigits(int entry);
    void loadRawDigits();
    void loadRawDigits(int entry);
    void loadClusters(int entry);
    void loadTracks(int entry);

    TTree* mHitsTree = nullptr;
    TTree* mDigiTree = nullptr;
    //TTree* mClusTree = nullptr;
    TTree* mTracTree = nullptr;

    // TEve-related members
    TEveElementList* mEvent = nullptr;
    TEveElement* getEveHits();
    TEveElement* getEveDigits();
    TEveElement* getEveClusters();
    TEveElement* getEveTracks();
} evdata;

void Data::setHitsTree(TTree* tree)
{
  if (tree == nullptr) {
    std::cerr << "No tree for hits!\n";
    return;
  }
  mHitsBuffer.resize(o2::tpc::Constants::MAXSECTOR);
  for(int sector = 0; sector < o2::tpc::Constants::MAXSECTOR; sector++) {
    std::stringstream hitsStr;
    hitsStr << "TPCHitsShiftedSector" << sector;
    tree->SetBranchAddress(hitsStr.str().c_str(), &(mHitsBuffer[sector]));
  }

  mHitsTree = tree;
}

// Based on: macro/analyzeHits::analyseTPC()
void Data::loadHits(int entry)
{
  static int lastLoaded = -1;

  if (mHitsTree == nullptr)
    return;

  int eventsCount = mHitsTree->GetBranch("TPCHitsShiftedSector0")->GetEntries();
  std::cout << "Hits: Number of events: " << eventsCount << std::endl;
  if ((entry < 0) || (entry >= eventsCount)) {
    std::cerr << "Hits: Out of event range ! " << entry << '\n';
    return;
  }
  if (entry != lastLoaded) {
    for(int sector = 0; sector < o2::tpc::Constants::MAXSECTOR; sector++) {
      mHitsBuffer[sector]->clear();
    }
    mHitsTree->GetEntry(entry);
    lastLoaded = entry;
  }

  int size = 0;
  int first = 0;
  mHits.resize(o2::tpc::Constants::MAXSECTOR);
  for(int sector = 0; sector < o2::tpc::Constants::MAXSECTOR; sector++) {
    int last = mHitsBuffer[sector]->size();
    mHits[sector] = gsl::make_span(&(*mHitsBuffer[sector])[first], last - first);
    size += mHits[sector].size();
  }

  std::cout << "Number of TPC Hits: " << size << '\n';
}

void Data::setDigiTree(TTree* tree)
{
    if (tree == nullptr) {
        std::cerr << "No tree for digits!\n";
        return;
    }

    mDigitBuffer.resize(o2::tpc::Constants::MAXSECTOR);
    for(int sector = 0; sector < o2::tpc::Constants::MAXSECTOR; sector++) {
      std::stringstream digiStr;
      digiStr << "TPCDigit_" << sector;
      tree->SetBranchAddress(digiStr.str().c_str(), &(mDigitBuffer[sector]));
    }

    mDigiTree = tree;

    int time = 0;
    mDigiTree->GetEntry(0); // Assuming that there will be 1 entry for all data
    for(int sector = 0; sector < o2::tpc::Constants::MAXSECTOR; sector++) {
      for(int i = 0; i < mDigitBuffer[sector]->size(); i++) {
        std::cout << "Digits: Sector: " << sector << " number of digits: " << mDigitBuffer[sector]->size() << '\n';
        int digitTime = (*mDigitBuffer[sector])[i].getTimeStamp();
        if(digitTime > time) time = digitTime;
        std::cout << "Digits: Sector: " << sector << " last digit time: " << digitTime << " current max time: " << time << '\n';
      }
    }
    int eventCount = time / (2 * tpcReadoutCycle) + (time % (2 * tpcReadoutCycle) > 0 ? 1 : 0);
    std::cout << "Digits: Current number of events: " << mEventsCount << " digits events: " << eventCount << '\n';
    if(eventCount > mEventsCount) mEventsCount = eventCount;
}

// Based on: CalibRawBase::processEventRawReader()
void Data::loadRawDigits()
{
  mDigitBuffer.resize(1);

  int size = 0;
  o2::tpc::PadPos padPos;
  while (std::shared_ptr<std::vector<uint16_t>> data = mRawReader->getNextData(padPos)) {
    if (!data)
      continue;

    CRU cru(mRawReader->getRegion());

    // row is local in region (CRU)
    const int row = padPos.getRow();
    const int pad = padPos.getPad();
    if (row == 255 || pad == 255)
      continue;

    int timeBin = 0;
    for (const auto& signalI : *data) {
      const float signal = float(signalI);
      mDigitBuffer[0]->emplace_back(cru, signal, row, pad, timeBin);
      size++;
      ++timeBin;
    }
    std::cout << "timeBin counts: " << timeBin << std::endl;
  }

  int first = 0;
  int last = mDigitBuffer[0]->size();
  mDigits[0] = gsl::make_span(&(*mDigitBuffer[0])[first], last - first);

  std::cout << "Number of TPC Digits: " << size << '\n';
}


void Data::loadRawDigits(int entry)
{
  if (mRawReader == nullptr)
    return;

  if ((entry < 0) || (entry >= mRawReader->getEventNumber())) {
    std::cerr << "Raw digits: Out of event range ! " << entry << '\n';
    return;
  }

  int eventId = mRawReader->loadEvent(entry);

  mLastEvent++;
  loadRawDigits();
}

void Data::loadDigits(int entry)
{
  static int lastLoaded = -1;

  if (mDigiTree == nullptr)
    return;

  if ((entry < 0) || (entry >= mEventsCount)) {
    std::cerr << "Digits: Out of event range ! " << entry << '\n';
    return;
  }
  if (entry != lastLoaded) {
//    for(int i = 0; i < o2::tpc::Constants::MAXSECTOR; i++) {
//      mDigitBuffer[i]->clear();
//    }
    //mDigiTree->GetEntry(entry);
    lastLoaded = entry;
  }

  int size = 0;
  int startTime = 2 * tpcReadoutCycle * entry;
  int endTime = startTime + 2 * tpcReadoutCycle;
  std::cout << "Digit time bounds: " << startTime << ", " << endTime << '\n';
  bool findFirst = true;
  mDigits.resize(o2::tpc::Constants::MAXSECTOR);
  for(int sector = 0; sector < o2::tpc::Constants::MAXSECTOR; sector++) {
    int first = 0;
    int last = 0;
    for(int i = 0; i < mDigitBuffer[sector]->size(); i++) {
      int digitTime = (*mDigitBuffer[sector])[i].getTimeStamp();
      if(findFirst && digitTime >= startTime && digitTime <= endTime) {
        first = i;
        findFirst = false;
      }
      else if(!findFirst && digitTime > endTime) {
        last = i;
        break;
      }
    }
    mDigits[sector] = gsl::make_span(&(*mDigitBuffer[sector])[first], last - first);
    size += mDigits[sector].size();
    std::cout << "Digits: Sector: " << sector << " number of digits: " << mDigits[sector].size() << " indices: " << first << ", " << last << '\n';
  }

  std::cout << "Number of TPC Digits: " << size << '\n';
}

//void Data::setClusTree(TTree* tree)
//{
//    if (tree == nullptr) {
//        std::cerr << "No tree for clusters !\n";
//        return;
//    }
//    tree->SetBranchAddress("TPCCluster", &mClusterBuffer);
//    mClusTree = tree;
//}

void Data::loadClusters(int entry)
{
    static int lastLoaded = -1;

    if (mClusReader == nullptr)
        return;

    if ((entry < 0) || (entry >= mClusReader->getTreeSize())) {
        std::cerr << "Clusters: Out of event range ! " << entry << '\n';
        return;
    }
    if (entry != lastLoaded) {
        //mClusterBuffer->clear();
        mClusReader->read(entry);
        //mClusTree->GetEntry(event);
        lastLoaded = entry;
    }

    //int first = 0, last = mClusterBuffer->size();
    //mClusters = gsl::make_span(&(*mClusterBuffer)[first], last - first);

    // Based on: MatchTPCITS::loadTPCClustersChunk()
    mClusReader->fillIndex(*mClusterIndex.get(), mClusterBuffer, mClusterMCBuffer);
    mClusterIndexStruct = mClusterIndex.get();

    std::cout << "Number of TPC Clusters: " << mClusterIndexStruct->nClustersTotal << '\n';
}

void Data::setTracTree(TTree* tree)
{
    if (tree == nullptr) {
        std::cerr << "No tree for tracks !\n";
        return;
    }
    tree->SetBranchAddress("TPCTracks", &mTrackBuffer);
    mTracTree = tree;

    int time = 0;
    mTracTree->GetEntry(0); // Assuming that there will be 1 entry for all data
    for(int i = 0; i < mTrackBuffer->size(); i++) {
      std::cout << "Tracks: Number of tracks: " << mTrackBuffer->size() << '\n';
      int trackTime = (*mTrackBuffer)[i].getTime0();
      if(trackTime > time) time = trackTime;
      std::cout << "Tracks: Last track time: " << trackTime << " current max time: " << time << '\n';
    }
    int eventCount = time / (2 * tpcReadoutCycle) + (time % (2 * tpcReadoutCycle) > 0 ? 1 : 0);
    std::cout << "Tracks: Current number of events: " << mEventsCount << " tracks events: " << eventCount << '\n';
    if(eventCount > mEventsCount) mEventsCount = eventCount;
}

void Data::loadTracks(int entry)
{
    static int lastLoaded = -1;

    if (mTracTree == nullptr)
      return;

    std::cout << "Tracks: Number of events: " << mEventsCount << std::endl;
    if ((entry < 0) || (entry >= mEventsCount)) {
      std::cerr << "Tracks: Out of event range ! " << entry << '\n';
      return;
    }
    if (entry != lastLoaded) {
  //        mTrackBuffer->clear();
  //        mTracTree->GetEntry(entry);
      lastLoaded = entry;
    }

    int size = 0;
    int startTime = 2 * tpcReadoutCycle * entry;
    int endTime = startTime + 2 * tpcReadoutCycle;
    std::cout << "Track time bounds: " << startTime << ", " << endTime << '\n';
    bool findFirst = true;
    int first = 0;
    int last = 0;
    for (int i = 0; i < mTrackBuffer->size(); i++) {
      int trackTime = (*mTrackBuffer)[i].getTime0();
      if (findFirst && trackTime >= startTime && trackTime <= endTime) {
        first = i;
        findFirst = false;
      } else if (!findFirst && trackTime > endTime) {
        last = i;
        break;
      }
    }
    mTracks = gsl::make_span(&(*mTrackBuffer)[first], last - first);
    size += mTracks.size();
    std::cout << "Tracks: Number of tracks: " << mTracks.size() << " indices: " << first << ", " << last << '\n';

    std::cout << "Number of TPC Tracks: " << mTracks.size() << '\n';
}

void Data::loadData(int entry)
{
    loadHits(entry);
    if(evdata.getRawData()) {
      loadRawDigits(entry);
    } else {
      loadDigits(entry);
    }
    loadClusters(entry);
    loadTracks(entry);
}

TEveElement* Data::getEveHits()
{
  const auto& mapper = Mapper::instance();
  TEvePointSet* hits = new TEvePointSet("hits");
  hits->SetMarkerColor(kYellow);
  hits->SetMarkerSize(0.2);

  for(int i = 0; i < mHits.size(); i++) {
    for (const auto& hv : mHits[i]) {
      for(int j = 0; j < hv.getSize(); j++) {
        const auto& h = hv.getHit(j);
        hits->SetNextPoint(h.GetX(), h.GetY(), h.GetZ());
      }
    }
  }
  return hits;
}

TEveElement* Data::getEveDigits()
{
    CDBInterface::instance().setUseDefaults();
    const auto& mapper = Mapper::instance();
    SAMPAProcessing& sampa = SAMPAProcessing::instance();
    TEvePointSet* digits = new TEvePointSet("digits");
    digits->SetMarkerColor(kBlue);
    digits->SetMarkerSize(0.2);

    for(int i = 0; i < mDigits.size(); i++) {
      for (const auto& d : mDigits[i]) {
        const auto pad = mapper.globalPadNumber(PadPos(d.getRow(),
                                                       d.getPad()));
        const CRU cru(d.getCRU());
        float z = sampa.getZfromTimeBin(d.getTimeStamp(), cru.side());

        const auto& localXYZ = mapper.padCentre(pad);
        const auto globalXYZ = mapper.LocalToGlobal(localXYZ, cru.sector());
        // TODO: One needs event time0 to get proper z-coordinate
        digits->SetNextPoint(globalXYZ.X(), globalXYZ.Y(), z);
      }
    }
    return digits;
}

TEveElement* Data::getEveClusters()
{
    CDBInterface::instance().setUseDefaults();
    SAMPAProcessing& sampa = SAMPAProcessing::instance();
    const auto& mapper = Mapper::instance();
    TEvePointSet* clusters = new TEvePointSet("clusters");
    clusters->SetMarkerColor(kRed);
    clusters->SetMarkerSize(0.2);

    const auto& clusterRefs = mClusterIndexStruct->clusters;
    for(int sector = 0; sector < o2::tpc::Constants::MAXSECTOR; sector++) {
      Sector sec(sector);
      for(int row = 0; row < o2::tpc::Constants::MAXGLOBALPADROW; row++) {
        const auto& c = clusterRefs[sector][row];
        float z = sampa.getZfromTimeBin(c->getTime(), sec.side());

        const auto pad = mapper.globalPadNumber(PadPos(row, c->getPad()));
        const auto& localXYZ = mapper.padCentre(pad);
        const auto globalXYZ = mapper.LocalToGlobal(localXYZ, sector);
        // TODO: One needs event time0 to get proper z-coordinate
        clusters->SetNextPoint(globalXYZ.X(), globalXYZ.Y(), z);
      }
    }
    return clusters;
}

TEveElement* Data::getEveTracks()
{
    CDBInterface::instance().setUseDefaults();
    SAMPAProcessing& sampa = SAMPAProcessing::instance();
    const auto& mapper = Mapper::instance();
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

        TEvePointSet* tpoints = new TEvePointSet("tclusters");
        tpoints->SetMarkerColor(kGreen);
        tpoints->SetMarkerSize(0.2);

        int nc = rec.getNClusterReferences();
        while (nc--) {
            uint8_t sector, row;
            uint32_t clusterIndexInRow;
            rec.getClusterReference(nc, sector, row, clusterIndexInRow);
            const auto& cl = rec.getCluster(nc, *mClusterIndexStruct, sector, row);
            Sector sec(sector);
            const auto pad = mapper.globalPadNumber(PadPos(row, cl.getPad()));
            float z = sampa.getZfromTimeBin(cl.getTime(), sec.side());
            z -= rec.getTime0();
            const auto& localXYZ = mapper.padCentre(pad);
            const auto globalXYZ = mapper.LocalToGlobal(localXYZ, sector);
            // TODO: One needs event time0 to get proper z-coordinate
            tpoints->SetNextPoint(globalXYZ.X(), globalXYZ.Y(), z);
        }
        track->AddElement(tpoints);
    }
    tracks->MakeTracks();

    return tracks;
}

void Data::displayData(int entry)
{
    std::string ename("Event #");
    ename += std::to_string(entry);

    // Event display
    auto hits = getEveHits();
    auto digits = getEveDigits();
    auto clusters = getEveClusters();
    auto tracks = getEveTracks();
    delete mEvent;
    mEvent = new TEveElementList(ename.c_str());
    mEvent->AddElement(hits);
    mEvent->AddElement(digits);
    mEvent->AddElement(clusters);
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
    std::string rawfile = "o2sim.root"; // should be e.g. GBTx0_Run005 - not a file per se?
    std::string hitsfile = "o2sim.root";
    std::string digifile = "tpcdigits.root";
    std::string clusfile = "tpc-native-clusters.root";
    evdata.setRawData(false);
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
    file = TFile::Open(hitsfile.data());
    if (file && gFile->IsOpen()) {
      evdata.setHitsTree((TTree*)gFile->Get("o2sim"));
    } else
      std::cerr << "\nERROR: Cannot open file: " << hitsfile << "\n\n";

    if (evdata.getRawData()) {
      // TODO: Another raw file for TPC? If at all?
        std::ifstream* rawfileStream = new std::ifstream(rawfile.data(), std::ifstream::binary);
        if (rawfileStream->good()) {
            delete rawfileStream;
            std::cout << "Running with raw digits...\n";
            evdata.setRawReader(rawfile.data());
        } else
            std::cerr << "\nERROR: Cannot open file: " << rawfile << "\n\n";
    } else {
        file = TFile::Open(digifile.data());
        if (file && gFile->IsOpen()) {
            std::cout << "Running with MC digits...\n";
            evdata.setDigiTree((TTree*)gFile->Get("o2sim"));
        } else
            std::cerr << "\nERROR: Cannot open file: " << digifile << "\n\n";
    }

    evdata.setClusReader(clusfile);

    file = TFile::Open(tracfile.data());
    if (file && gFile->IsOpen()) {
      evdata.setTracTree((TTree*)gFile->Get("tpcrec"));
    }
    else
        std::cerr << "\nERROR: Cannot open file: " << tracfile << "\n\n";

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