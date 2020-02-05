#include "AODViewerSpec.h"

// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.


#include "Framework/DataProcessorSpec.h"
#include "Framework/DataRefUtils.h"
#include "Headers/DataHeader.h"
#include "Framework/Lifetime.h"
#include "Framework/ControlService.h"
#include <FairMQLogger.h>
#include "ReconstructionDataFormats/Vertex.h"

using namespace o2::framework;

namespace o2
{
  namespace viewer
  {

    DataProcessorSpec getAODViewerSpec()
    {
      // init function return a lambda taking a ProcessingContext
      auto doIt = [](ProcessingContext& pc) {
        // will be triggered whenever data arrives
        LOG(INFO) << "VIEWER GOT DATA -----------------";


        const auto secret = pc.inputs().get<int*>("input");

        LOG(INFO) << "RECEIVED " << *secret;
      };

      return DataProcessorSpec{
        /*ID*/ "AODViewer",
        /*INPUT CHANNELS*/ Inputs{InputSpec{"input", "AOD", "TRACKS", 0, Lifetime::Timeframe}},
        /*OUTPUT CHANNELS*/ Outputs{},
        /* ALGORITHM */  AlgorithmSpec(doIt),
        /* OPTIONS */Options{}};
    }
  } // namespace viewer
} // namespace o2
