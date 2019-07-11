//
// Created by jmy on 26.02.19.
//

#ifndef ALICE_O2_EVENTVISUALISATION_BASE_DATASOURCEOFFLINE_H
#define ALICE_O2_EVENTVISUALISATION_BASE_DATASOURCEOFFLINE_H

#include "EventVisualisationBase/DataSource.h"

#include <TString.h>

namespace o2  {
namespace event_visualisation {

class DataSourceOffline : public DataSource {
protected:
    TString fgESDFileName;
    TString  fgRawFileName;
    bool isOpen = kFALSE;
public:
    virtual void Open(TString ESDFileName) {
        this->fgESDFileName = ESDFileName;
    };
    virtual void OpenRawFile() {};
    Bool_t GotoEvent(Int_t event) override {};
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
