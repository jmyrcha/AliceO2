//
// Created by jmy on 09.07.19.
//

#include "EventVisualisationDetectors/DataSourceOfflineITS.h"

//#include "EventVisualisationView/MultiView.h"
//#include "ITSMFTReconstruction/ChipMappingITS.h"
//#include "ITSMFTReconstruction/DigitPixelReader.h"
//#include "ITSMFTReconstruction/RawPixelReader.h"
//#include "ITSMFTBase/SegmentationAlpide.h"
//#include "ITSMFTBase/Digit.h"
//#include "ITSBase/GeometryTGeo.h"
//#include "DataFormatsITSMFT/Cluster.h"
//#include "DataFormatsITSMFT/ROFRecord.h"
//#include "DataFormatsITS/TrackITS.h"
//#include "DetectorsCommonDataFormats/DetID.h"
//#include "CommonDataFormat/InteractionRecord.h"

//#include <iostream>
//#include <array>
//#include <algorithm>
//#include <fstream>

//#include <TFile.h>
//#include <TTree.h>
#include <TEveManager.h>
//#include <TEveBrowser.h>
//#include <TGButton.h>
//#include <TGNumberEntry.h>
//#include <TGFrame.h>
//#include <TGTab.h>
//#include <TGLCameraOverlay.h>
//#include <TEveFrameBox.h>
//#include <TEveQuadSet.h>
//#include <TEveTrans.h>
//#include <TEvePointSet.h>
//#include <TEveTrackPropagator.h>
//#include <TEveTrack.h>
//#include <TGenericClassInfo.h>
//#include <TEveElement.h>
//#include <Rtypes.h>
//#include <gsl/span>

extern TEveManager* gEve;

using namespace o2::itsmft;

namespace o2 {
namespace event_visualisation {

DataSourceOfflineITS::DataSourceOfflineITS():
DataSourceOffline()
{
  fgRawFileName = TString("itsdigits.root");
  fgDigitsFileName = TString("itsdigits.root");
  fgClustersFileName = TString("o2clus_its.root");
  fgTracksFileName = TString("o2trac_its.root");
}

void DataSourceOfflineITS::OpenRawFile()
{
  // Here should be check and search in default raw file paths
  // AliRoot/EVE/EveBase/AliEveDataSourceOffline.cxx lines 342-369

  std::ifstream* rawfile = new std::ifstream(fgRawFileName.Data(), std::ifstream::binary);
  if (rawfile->good()) {
    delete rawfile;
    std::cout << "Running with raw digits...\n";
    auto reader = new RawPixelReader<ChipMappingITS>();
    reader->openInput(fgRawFileName.Data());
    mPixelReader = reader;
    reader->getNextChipData(mChipData);
    mIR = mChipData.getInteractionRecord();
  } else
    std::cerr << "\nERROR: Cannot open file: " << fgRawFileName << "\n\n";
}

void DataSourceOfflineITS::OpenDigitsFile()
{
  // Here should be check and search in default digits file paths

  TFile* file = TFile::Open(fgDigitsFileName.Data());
  if (file && gFile->IsOpen()) {
    file->Close();
    std::cout << "Running with MC digits...\n";
    auto reader = new DigitPixelReader();
    reader->openInput(fgDigitsFileName.Data(), o2::detectors::DetID("ITS"));
    reader->init();
    reader->readNextEntry();
    mPixelReader = reader;
    mPixelReader->getNextChipData(mChipData);
    mIR = mChipData.getInteractionRecord();
  } else
    std::cerr << "\nERROR: Cannot open file: " << fgDigitsFileName << "\n\n";
}

void DataSourceOfflineITS::OpenClustersFile()
{
  // Here should be check and search in default clusters file paths

  TFile* file = TFile::Open(fgClustersFileName.Data());
  if (file && gFile->IsOpen())
  {
    TTree* tree = (TTree*)gFile->Get("o2sim");
    if (tree == nullptr) {
      std::cerr << "No tree for clusters !\n";
      return;
    }
    tree->SetBranchAddress("ITSCluster", &mClusterBuffer);
    mClusTree = tree;

    TTree* roft = (TTree*)gFile->Get("ITSClustersROF");
    if (roft != nullptr) {
      std::vector<o2::itsmft::ROFRecord>* roFrames = &mClustersROF;
      roft->SetBranchAddress("ITSClustersROF", &roFrames);
      roft->GetEntry(0);
    }
  }
  else
    std::cerr << "ERROR: Cannot open file: " << fgClustersFileName << "\n\n";
}

void DataSourceOfflineITS::OpenTracksFile()
{
  // Here should be check and search in default tracks file paths

  TFile* file = TFile::Open(fgTracksFileName.Data());
  if (file && gFile->IsOpen())
  {
    TTree* tree = (TTree*)gFile->Get("o2sim");
    if (tree == nullptr) {
      std::cerr << "No tree for tracks !\n";
      return;
    }
    tree->SetBranchAddress("ITSTrack", &mTrackBuffer);
    Info("OpenTracksFile", "Setting branch address for track clusters...");
    tree->SetBranchAddress("ITSTrackClusIdx", &mClIdxBuffer);
    mTracTree = tree;

    TTree* roft = (TTree*)gFile->Get("ITSTracksROF");
    if (roft != nullptr) {
      std::vector<o2::itsmft::ROFRecord>* roFrames = &mTracksROF;
      roft->SetBranchAddress("ITSTracksROF", &roFrames);
      roft->GetEntry(0);
    }
  }
  else
    std::cerr << "\nERROR: Cannot open file: " << fgTracksFileName << "\n\n";
}

void DataSourceOfflineITS::Open(TString fileName)
{
  DataSourceOffline::Open(fileName);
  std::cout << "DataSourceOfflineITS::Open()" << std::endl;
}

Int_t DataSourceOfflineITS::GotoEvent(Int_t ev)
{
  Warning("GotoEvent", "GOTOEVENT");
  if (ev < 0 || ev >= this->GetEventCount()) {
    Warning("GotoEvent", "Invalid event id %d.", ev);
    return kFALSE;
  }

  LoadDigits(ev);
  LoadClusters(ev);
  LoadTracks(ev);
  return kTRUE;
}

void DataSourceOfflineITS::LoadDigits(Int_t ev)
{
  if (mPixelReader == nullptr)
    return;

  for (; mLastEvent < ev; mLastEvent++) {
    auto ir = mChipData.getInteractionRecord();
    do {
      if (!mPixelReader->getNextChipData(mChipData))
        return;
      ir = mChipData.getInteractionRecord();
    } while (mIR == ir);
    mIR = ir;
  }
  mLastEvent++;
  LoadDigits();
}

void DataSourceOfflineITS::LoadDigits()
{
  auto ir = mChipData.getInteractionRecord();
  std::cout << "orbit/crossing: " << ' ' << ir.orbit << '/' << ir.bc << '\n';

  mDigits.clear();

  do {
    auto chipID = mChipData.getChipID();
    auto pixels = mChipData.getData();
    for (auto& pixel : pixels) {
      auto col = pixel.getCol();
      auto row = pixel.getRow();
      mDigits.emplace_back(chipID, 0, row, col);
    }
    if (!mPixelReader->getNextChipData(mChipData))
      return;
    ir = mChipData.getInteractionRecord();
  } while (mIR == ir);
  mIR = ir;

  std::cout << "Number of ITSDigits: " << mDigits.size() << '\n';
}

void DataSourceOfflineITS::LoadClusters(Int_t ev)
{
  static int lastLoaded = -1;

  if (mClusTree == nullptr)
    return;

  auto event = ev; // If no RO frame informaton available, assume one entry per a RO frame.
  if (!mClustersROF.empty()) {
    if ((event < 0) || (event >= (int)mClustersROF.size())) {
      std::cerr << "Clusters: Out of event range ! " << event << '\n';
      return;
    }
    auto rof = mClustersROF[ev];
    event = rof.getROFEntry().getEvent();
  }
  if ((event < 0) || (event >= mClusTree->GetEntries())) {
    std::cerr << "Clusters: Out of event range ! " << event << '\n';
    return;
  }
  if (event != lastLoaded) {
    mClusterBuffer->clear();
    mClusTree->GetEntry(event);
    lastLoaded = event;
  }

  int first = 0, last = mClusterBuffer->size();
  if (!mClustersROF.empty()) {
    auto rof = mClustersROF[ev];
    first = rof.getROFEntry().getIndex();
    last = first + rof.getNROFEntries();
  }
  mClusters = gsl::make_span(&(*mClusterBuffer)[first], last - first);

  std::cout << "Number of ITSClusters: " << mClusters.size() << '\n';
}

void DataSourceOfflineITS::LoadTracks(Int_t ev)
{
  static int lastLoaded = -1;

  if (mTracTree == nullptr)
    return;

  auto event = ev; // If no RO frame informaton available, assume one entry per a RO frame.
  if (!mTracksROF.empty()) {
    if ((event < 0) || (event >= (int)mTracksROF.size())) {
      std::cerr << "Clusters: Out of event range ! " << event << '\n';
      return;
    }
    auto rof = mTracksROF[ev];
    event = rof.getROFEntry().getEvent();
  }
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
  if (!mTracksROF.empty()) {
    auto rof = mTracksROF[ev];
    first = rof.getROFEntry().getIndex();
    last = first + rof.getNROFEntries();
  }
  mTracks = gsl::make_span(&(*mTrackBuffer)[first], last - first);

  std::cout << "Number of ITSTracks: " << mTracks.size() << '\n';
}

}
}
