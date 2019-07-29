// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DataReaderITS.h
/// \brief ITS Detector-specific reading from file(s)
/// \author julian.myrcha@cern.ch

#ifndef O2EVE_EVENTVISUALISATION_DETECTORS_DATAREADERITS_H
#define O2EVE_EVENTVISUALISATION_DETECTORS_DATAREADERITS_H

#include <TFile.h>
#include "EventVisualisationBase/DataReader.h"

namespace o2 {
    namespace event_visualisation {


        class DataReaderITS : public DataReader {
        private:
            Int_t fMaxEv;
            TFile *clusFile;
            TFile *tracFile;
        public:
            DataReaderITS();
            void open() override;
            Int_t GetEventCount() override;

            /*
            zwraca TLista
             */

            TObject* getEventData(int no) override;
        };

    }
}

#endif //O2EVE_DATAREADERITS_H
