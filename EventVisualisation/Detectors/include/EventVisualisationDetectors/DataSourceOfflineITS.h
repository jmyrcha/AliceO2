//
// Created by jmy on 09.07.19.
//

#ifndef O2EVE_DATASOURCEOFFLINEITS_H
#define O2EVE_DATASOURCEOFFLINEITS_H

#include "EventVisualisationBase/DataSourceOffline.h"

#include "ITSMFTReconstruction/DigitPixelReader.h"
#include "ITSMFTReconstruction/RawPixelReader.h"
//#include "ITSMFTReconstruction/ChipMappingITS.h"
//#include "ITSMFTBase/SegmentationAlpide.h"
#include "ITSMFTBase/Digit.h"
#include "ITSBase/GeometryTGeo.h"
#include "DataFormatsITSMFT/Cluster.h"
//#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITS/TrackITS.h"

#include <TFile.h>
#include <TTree.h>
#include <TEvePointSet.h>
#include <TEveTrackPropagator.h>
#include <TEveTrack.h>

using namespace o2::itsmft;

namespace o2 {
namespace event_visualisation {

class DataSourceOfflineITS : public DataSourceOffline {
public:

    DataSourceOfflineITS();

    void Open(TString fileName) override;

    void OpenRawFile() override;
    void OpenDigitsFile() override;
    void OpenClustersFile() override;
    void OpenTracksFile() override;

    Int_t GotoEvent(Int_t ev) override;
    Int_t GetEventCount() override { return 10; } // dummy value for now

    void LoadDigits(Int_t ev);
    void LoadDigits();
    void LoadClusters(Int_t ev);
    void LoadTracks(Int_t ev);

    gsl::span<Cluster> GetClusters() { return mClusters; }
    std::vector<Int_t>* GetClIdxBuffer() { return mClIdxBuffer; }
    gsl::span<its::TrackITS> GetTracks() { return mTracks; }

 private:
  PixelReader* mPixelReader = nullptr;
  ChipPixelData mChipData;
  o2::InteractionRecord mIR;
  TTree* mDigiTree = nullptr;
  TTree* mClusTree = nullptr;
  TTree* mTracTree = nullptr;

  TEveElementList* mEvent = nullptr;
  TEveElementList* mChip = nullptr;

  std::vector<Digit> mDigits;
  std::vector<Cluster>* mClusterBuffer = nullptr;
  gsl::span<Cluster> mClusters;
  std::vector<ROFRecord> mClustersROF;
  std::vector<its::TrackITS>* mTrackBuffer = nullptr;
  std::vector<Int_t>* mClIdxBuffer = nullptr;
  gsl::span<its::TrackITS> mTracks;
  std::vector<ROFRecord> mTracksROF;
};

}
}

#endif //O2EVE_DATASOURCEOFFLINEITS_H
