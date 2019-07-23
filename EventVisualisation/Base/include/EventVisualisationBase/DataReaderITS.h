//
// Created by jmy on 09.07.19.
//

#ifndef O2EVE_DATAREADERITS_H
#define O2EVE_DATAREADERITS_H

#include <EventVisualisationBase/DataReader.h>

namespace o2 {
namespace event_visualisation {


class DataReaderITS : public DataReader {
private:
public:
    DataReaderITS();
    void open() override;
    Bool_t GotoEvent(Int_t ev);

    Int_t GetEventCount() override;
    TObject* getEventData(int no) override;
};

}
}

#endif //O2EVE_DATAREADERITS_H
