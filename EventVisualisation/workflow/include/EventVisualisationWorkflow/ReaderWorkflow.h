//
// Created by maja on 16.07.19.
//

#ifndef O2EVE_READERWORKFLOW_H
#define O2EVE_READERWORKFLOW_H

#include "Framework/WorkflowSpec.h"
#include <vector>
#include <string>

namespace o2
{
namespace event_visualisation
{
namespace reader_workflow
{

/// define input and output types of the workflow
enum struct InputType { Digitizer, // directly read digits from channel {TPC:DIGITS}
                        Digits,    // read digits from file
                        Raw,       // read hardware clusters in raw page format from file
                        Clusters,  // read native clusters from file
};

/// create the workflow for TPC reconstruction
framework::WorkflowSpec getWorkflow(std::vector<int> const& tpcSectors,           //
                                    std::vector<int> const& laneConfiguration,    //
                                    bool propagateMC = true, unsigned nLanes = 1, //
                                    std::string const& cfgInput = "digitizer",    //
                                    std::string const& cfgOutput = "tracks"       //
);

} // end namespace reader_workflow
} // end namespace event_visualisation
} // end namespace o2

#endif //O2EVE_READERWORKFLOW_H
