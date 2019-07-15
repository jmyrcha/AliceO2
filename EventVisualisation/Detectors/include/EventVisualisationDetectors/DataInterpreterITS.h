//
// Created by jmy on 25.06.19.
//

#ifndef ALICE_O2_EVENTVISUALISATION_BASE_DATAINTERPRETERITS_H
#define ALICE_O2_EVENTVISUALISATION_BASE_DATAINTERPRETERITS_H

/// DataInterpreterITS prepares random events
///
/// This class overrides DataInterpreter and implements method
/// returning visualisation objects representing data from ITS file
/// with tracks colored by PID only.

#include "EventVisualisationBase/DataInterpreter.h"
#include "EventVisualisationBase/EventManager.h"
#include "EventVisualisationBase/VisualisationConstants.h"

#include "EventVisualisationDetectors/DataSourceOfflineITS.h"

namespace o2 {
namespace event_visualisation {


class DataInterpreterITS : public DataInterpreter {
public:
    /// Default constructor
    DataInterpreterITS();

    /// Default destructor
    ~DataInterpreterITS() final;

    void Open(EDataType type) override;
    void Close() override { dataSource->Close(); }

    Int_t GetEventCount() override { return dataSource->GetEventCount(); };
    Int_t GotoEvent(Int_t ev) override;

    /// Returns a list of random tracks colored by PID
    TEveElement *interpretDataForType(EDataType type) final;

    TEveElement* loadTracks();

    DataSourceOfflineITS *getDataSource() { return dataSource; }
    void setDataSource(DataSourceOfflineITS *dataSource) { this->dataSource = dataSource; }

 protected:
    DataSourceOfflineITS *dataSource = nullptr;
};

}
}

#endif //ALICE_O2_EVENTVISUALISATION_BASE_DATAINTERPRETERITS_H
