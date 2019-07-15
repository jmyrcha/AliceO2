//
// Created by jmy on 26.02.19.
//

#ifndef ALICE_O2_EVENTVISUALISATION_BASE_DATASOURCEOFFLINE_H
#define ALICE_O2_EVENTVISUALISATION_BASE_DATASOURCEOFFLINE_H

#include "EventVisualisationBase/DataSource.h"

#include "ReconstructionDataFormats/Track.h"

#include <TString.h>
//#include <gsl/span>

namespace o2  {
namespace event_visualisation {

class DataSourceOffline : public DataSource {
protected:
    TString fgESDFileName;
    TString  fgRawFileName;
    TString  fgDigitsFileName;
    TString  fgClustersFileName;
    TString fgTracksFileName;

    Bool_t fIsOpen = kFALSE;
    Int_t mLastEvent = 0;
public:
    virtual void Open(TString ESDFileName) override;
    virtual void OpenRawFile() {};
    virtual void OpenDigitsFile() {};
    virtual void OpenClustersFile() {};
    virtual void OpenTracksFile() {};

    Int_t GotoEvent(Int_t event) override {};
    void NextEvent() override {};

    /// Default constructor
    DataSourceOffline();

    ~DataSourceOffline() override {};

    /// Deleted copy constructor
    DataSourceOffline(DataSourceOffline const&) = delete;

    /// Deleted assigment operator
    void operator=(DataSourceOffline const&) = delete;
};

}
}


#endif //ALICE_O2_EVENTVISUALISATION_BASE_DATASOURCEOFFLINE_H
