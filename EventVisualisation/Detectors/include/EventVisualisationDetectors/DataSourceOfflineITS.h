//
// Created by jmy on 09.07.19.
//

#ifndef O2EVE_DATASOURCEOFFLINEITS_H
#define O2EVE_DATASOURCEOFFLINEITS_H

#include "EventVisualisationBase/DataSourceOffline.h"
#include "ITSMFTReconstruction/DigitPixelReader.h"
#include "ITSMFTReconstruction/RawPixelReader.h"

using namespace o2::itsmft;

namespace o2 {
namespace event_visualisation {

class DataSourceOfflineITS : public DataSourceOffline {
public:

    DataSourceOfflineITS();

    void Open(TString fileName) override;

    void OpenRawFile() override;

    Bool_t GotoEvent(Int_t ev);

 private:
  PixelReader* mPixelReader = nullptr;
  ChipPixelData mChipData;
  o2::InteractionRecord mIR;

};

}
}

#endif //O2EVE_DATASOURCEOFFLINEITS_H
