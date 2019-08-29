//
// Created by maja on 16.07.19.
//

#include "Framework/WorkflowSpec.h"
#include "Framework/ConfigParamSpec.h"
#include "Algorithm/RangeTokenizer.h"

#include <string>
#include <stdexcept>
#include <unordered_map>

#include "EventVisualisationWorkflow/WriterWorkflow.h"

// add workflow options, note that customization needs to be declared before
// including Framework/runDataProcessing
void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<o2::framework::ConfigParamSpec> options{
    { "input-type", o2::framework::VariantType::String, "digits", { "digitizer, digits, raw, clusters" } },
    { "output-type", o2::framework::VariantType::String, "tracks", { "digits, raw, clusters, tracks" } },
    { "disable-mc", o2::framework::VariantType::Bool, false, { "disable sending of MC information" } },
    { "tpc-sectors", o2::framework::VariantType::String, "0-35", { "TPC sector range, e.g. 5-7,8,9" } },
    { "tpc-lanes", o2::framework::VariantType::Int, 1, { "number of parallel lanes up to the tracker" } },
  };
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h" // the main driver

using namespace o2::framework;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Possibly read the configuration and set parameters
  auto tpcSectors = o2::RangeTokenizer::tokenize<int>(cfgc.options().get<std::string>("tpc-sectors"));
  // the lane configuration defines the subspecification ids to be distributed among the lanes.
  std::vector<int> laneConfiguration;
  auto nLanes = cfgc.options().get<int>("tpc-lanes");
  auto inputType = cfgc.options().get<std::string>("input-type");
  if (inputType == "digitizer") {
    // the digitizer is using a different lane setup so we have to force this for the moment
    laneConfiguration.resize(nLanes);
    std::iota(laneConfiguration.begin(), laneConfiguration.end(), 0);
  } else {
    laneConfiguration = tpcSectors;
  }

  return o2::event_visualisation::writer_workflow::getWorkflow(tpcSectors,                                    // sector configuration
                                                               laneConfiguration,                             // lane configuration
                                                               not cfgc.options().get<bool>("disable-mc"),    //
                                                               nLanes,                                        //
                                                               inputType,                                     //
                                                               cfgc.options().get<std::string>("output-type") //
  );
}
