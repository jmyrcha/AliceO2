
#include "AODReaderSpec.h"
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

    DataProcessorSpec getAODReaderSpec()
    {
      // init function return a lambda taking a ProcessingContext
      auto doIt = [](ProcessingContext& pc) {
       // user logic goes here
        LOG(INFO) << "HELLO FROM READER ====================================";

        o2::dataformats::Vertex vertex;
        vertex.setX(-11.);

        int secret = 11;
        // data gets send here
        pc.outputs().snapshot(OutputRef{"outtracks"},secret);

      };

      return DataProcessorSpec{
              /*ID*/ "AODReader",
              /*INPUT CHANNELS*/ Inputs{},
              /*OUTPUT CHANNELS*/ Outputs{{{"outtracks"}, "AOD", "TRACKS", 0, Lifetime::Timeframe}},
              /* ALGORITHM */  AlgorithmSpec(doIt),
              /* OPTIONS */Options{}};
    }
  } // namespace viewer
} // namespace o2