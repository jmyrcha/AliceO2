//
// Created by maja on 15.07.19.
//

#ifndef O2EVE_SIMPLEWORKFLOW_H
#define O2EVE_SIMPLEWORKFLOW_H

#include "Framework/WorkflowSpec.h"
#include "Framework/ConfigContext.h"

using namespace o2::framework;

namespace o2
{
namespace event_visualisation
{
namespace simple_workflow
{

// Get spec for start of the job - creating TPC fake clusters
DataProcessorSpec getInputSpec();

// Helper to create algorithm spec
AlgorithmSpec simplePipe(std::string const& what);
// Get specs to analyse and transform input data
DataProcessorSpec getTPCTrackingSpec();

// Get spec to write the data (dummy here)
DataProcessorSpec getWriterSpec();

// Create simple workflow
WorkflowSpec getWorkflow();

} // end namespace simple_workflow
} // end namespace event_visualisation
} // end namespace o2

// Main function - defining a workflow.
// This function cannot be defined inside any other namespace
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc);

#endif //O2EVE_SIMPLEWORKFLOW_H
