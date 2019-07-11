//
// Created by jmy on 26.02.19.
//

#ifndef ALICE_O2_EVENTVISUALISATION_BASE_DATASOURCE_H
#define ALICE_O2_EVENTVISUALISATION_BASE_DATASOURCE_H

#include <TQObject.h>

namespace o2 {
namespace event_visualisation {

class DataSource : public TQObject {
public:
    virtual Bool_t GotoEvent(Int_t event) {};
    virtual void NextEvent() {};

    DataSource() = default;

    /// Default destructor
    virtual ~DataSource() = default;

    /// Deleted copy constructor
    DataSource(DataSource const &) = delete;

    /// Deleted assigment operator
    void operator=(DataSource const &) = delete;
};

}
}

#endif //ALICE_O2_EVENTVISUALISATION_BASE_DATASOURCE_H
