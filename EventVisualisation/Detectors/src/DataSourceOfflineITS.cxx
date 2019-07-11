//
// Created by jmy on 09.07.19.
//

//#include "EventVisualisationView/MultiView.h"
#include "EventVisualisationDetectors/DataSourceOfflineITS.h"
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
//#include <TEveManager.h>
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

using namespace o2::itsmft;

namespace o2 {
namespace event_visualisation {

DataSourceOfflineITS::DataSourceOfflineITS():
DataSourceOffline()
{
  fgRawFileName = TString("itsdigits.root");
}

void DataSourceOfflineITS::OpenRawFile()
{
  // Here should be check and search in default raw file paths
  // AliRoot/EVE/EveBase/AliEveDataSourceOffline.cxx lines 342-369

  TString digifile = fgRawFileName;
  std::ifstream* rawfile = new std::ifstream(digifile.Data(), std::ifstream::binary);
  if (rawfile->good()) {
    delete rawfile;
    std::cout << "Running with raw digits...\n";
    auto reader = new RawPixelReader<ChipMappingITS>();
    reader->openInput(digifile.Data());
    mPixelReader = reader;
    reader->getNextChipData(mChipData);
    mIR = mChipData.getInteractionRecord();
  } else
    std::cerr << "\nERROR: Cannot open file: " << digifile << "\n\n";
}

void DataSourceOfflineITS::Open(TString fileName)
{
  DataSourceOffline::Open(fileName);
}

Bool_t DataSourceOfflineITS::GotoEvent(Int_t ev)
{
  return kTRUE;
}

}
}
