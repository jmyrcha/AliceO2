//
// Created by maja on 15.07.19.
//
// From o2SyncReconstructionDummy.cxx - shortened
//

#include "EventVisualisationWorkflow/SimpleWorkflow.h"

#include "Framework/ControlService.h"

#include <chrono>

using namespace o2::framework;

// Add workflow options, note that customization needs to be declared before
// including Framework/runDataProcessing.
//void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
//{
//  std::vector<o2::framework::ConfigParamSpec> options {
//    // option name      // option value type  // default value    // description, available options
//    { "epn-roundrobin-delay", VariantType::Int, 28, {"Fake delay for waiting from the network for a new timeframe"} },
//    { "input-type", o2::framework::VariantType::String, "digits", { "digitizer, digits, raw, clusters" } },
//    { "output-type", o2::framework::VariantType::String, "tracks", { "digits, raw, clusters, tracks" } }
//  };
//  std::swap(workflowOptions, options);
//}

#include "Framework/runDataProcessing.h" // the main driver

// Main function - defining a workflow.
// This function cannot be inside any other namespace.
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Possibly read the configuration and set parameters.

  return o2::event_visualisation::simple_workflow::getWorkflow();
}

namespace o2
{
namespace event_visualisation
{
namespace simple_workflow
{

WorkflowSpec getWorkflow()
{
  WorkflowSpec specs;

  // Possibly checking input arguments here.
  // E.g. you can define different data processor specs depending on the input arguments.

  // Add specs to workflow.
  specs.emplace_back(getInputSpec());
  specs.emplace_back(getTPCTrackingSpec());
  specs.emplace_back(getWriterSpec());

  return std::move(specs);
}

// Example of how to define a spec with lambda functions
DataProcessorSpec getInputSpec()
{
  std::string processorName = "flp-dpl-source";

  // Defining inputs and output with lambda functions.
  auto createInputSpecs = []() {
    std::vector<InputSpec> inputSpecs{
      // No input here.
      // Fields similar as for OutputSpec.
    };
    return std::move(inputSpecs);
  };

  auto createOutputSpecs = []() {
    std::vector<OutputSpec> outputSpecs{
      // label     // data origin // object type
      OutputSpec{ { "tpc-cluster" }, "TPC", "CLUSTERS" },
      // label - used to refer to specific outputs,
      //         InputSpec has a single string instead of a vector of strings.
      // data origin - usually a detector??
      // object type - CLUSTERS, TRACKS etc.
      // optional parameters: generic subspecification (e.g. to address TPC sector),
      //                      lifetime of the message e.g. Lifetime::Timeframe.
      // O2 attaches to a message: its origin, a description and a generic subspecification.
    };
    return std::move(outputSpecs);
  };

  // AlgorithmSpec can utilise 3 callbacks:
  // using ProcessCallback = std::function<void(ProcessingContext &)>;
  // using InitCallback = std::function<ProcessCallback(InitContext &)>;
  // using ErrorCallback = std::function<void(ErrorContext &)>;

  // InitCallback - to be invoked at (re-)start of the job.
  // Useful for stateful processing - to init some state than is then passed to ProcessCallback.
  // Setup contains options set by the user - according to ConfigParamSpecs defined above.
  // ProcessCallback definition is contained in InitCallback definition, and this is the function returned.
  auto initFunction = [](InitContext& setup) {
    // Setting the states of some variables.
    // They need to be shared_ptr to be processed properly.
    // Getting a value from configuration options (e.g. command-line argument)
    auto delay = setup.options().get<int>("epn-roundrobin-delay");
    auto first = std::make_shared<bool>(true);
    // ProcessCallback - actual computation made in the workflow.
    return [delay, first](ProcessingContext& ctx) {
      // Example of how to create outputs.
      // Make 1 integer and route it to the output referred by OutputRef.
      // Here we refer to specific OutputSpecs by labels.
      ctx.outputs().make<int>(OutputRef{ "tpc-cluster" }, 1);
      //ctx.services().get<ControlService>().readyToQuit(false);
    };
  };

  // Defining input (e.g. command-line) arguments for this processor
  auto createConfig = []() {
    return ConfigParamSpec{ "epn-roundrobin-delay", VariantType::Int, 1, { "Fake delay for waiting from the network for a new timeframe" } };
  };

  return DataProcessorSpec{ processorName,
                            { createInputSpecs() },
                            { createOutputSpecs() },
                            AlgorithmSpec(initFunction),
                            { createConfig() } };
}

// Helper to create algorithm spec
AlgorithmSpec simplePipe(std::string const& what)
{
  return AlgorithmSpec{ [what](InitContext& ic) {
    // Set initial states
    auto messageSize = ic.options().get<int>("size");

    return [what, messageSize](ProcessingContext& ctx) {
      // Send to output <messageSize> chars
      auto msg = ctx.outputs().make<char>(OutputRef{ what }, messageSize);
    };
  } };
}

// Examples of how to declare a spec in a more declarative way
DataProcessorSpec getTPCTrackingSpec()
{
  return DataProcessorSpec{ "tpc-tracking",
                            { InputSpec{ "clusters", "TPC", "CLUSTERS" } },
                            { OutputSpec{ { "tracks" }, "TPC", "TRACKS" } },
                            simplePipe("tracks"),
                            { ConfigParamSpec{ "size", VariantType::Int, 1, { "Size of the output message" } } } };
}

DataProcessorSpec getWriterSpec()
{
  return DataProcessorSpec{
    "writer",
    { InputSpec{ "tracks-tpc", "TPC", "TRACKS" } },
    // No outputs
    {},
    AlgorithmSpec{ [](ProcessingContext& ctx) {
      ctx.services().get<ControlService>().readyToQuit(false);
    } }
  };
}

} // end namespace simple_workflow
} // end namespace event_visualisation
} // end namespace o2