//
// Created by maja on 16.07.19.
//

#include "Framework/WorkflowSpec.h"
#include "Framework/DataSpecUtils.h"
#include "DPLUtils/MakeRootTreeWriterSpec.h"
#include "DataFormatsTPC/TPCSectorHeader.h"
#include "Algorithm/RangeTokenizer.h"
#include "TPCBase/Digit.h"
#include "DataFormatsTPC/Constants.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"

#include "FairMQLogger.h"

#include "EventVisualisationWorkflow/ReaderWorkflow.h"
#include "EventVisualisationWorkflow/PublisherSpec.h"

#include <unordered_map>
#include <string>
#include <fstream>
#include <stdexcept>
#include <algorithm> // std::find
#include <tuple>     // make_tuple

namespace o2
{
namespace event_visualisation
{
namespace reader_workflow
{

using namespace framework;

template <typename T>
using BranchDefinition = MakeRootTreeWriterSpec::BranchDefinition<T>;

const std::unordered_map<std::string, InputType> InputMap{
  { "digitizer", InputType::Digitizer },
  { "digits", InputType::Digits },
  { "raw", InputType::Raw },
  { "clusters", InputType::Clusters },
};

framework::WorkflowSpec getWorkflow(std::vector<int> const& tpcSectors, std::vector<int> const& laneConfiguration, bool propagateMC, unsigned nLanes, std::string const& cfgInput, std::string const& cfgOutput)
{
  InputType inputType;

  try {
    inputType = InputMap.at(cfgInput);
  } catch (std::out_of_range&) {
    throw std::invalid_argument(std::string("invalid input type: ") + cfgInput);
  }

  WorkflowSpec specs;

  if (inputType == InputType::Digits) {
    specs.emplace_back(o2::event_visualisation::getPublisherSpec(PublisherConf{
                                                                   "tpc-digit-reader",
                                                                   "o2sim",
                                                                   { "digitbranch", "TPCDigit", "Digit branch" },
                                                                   { "mcbranch", "TPCDigitMCTruth", "MC label branch" },
                                                                   OutputSpec{ "TPC", "DIGITS" },
                                                                   OutputSpec{ "TPC", "DIGITSMCTR" },
                                                                   tpcSectors,
                                                                   laneConfiguration,
                                                                 },
                                                                 propagateMC));
  } else if (inputType == InputType::Raw) {
    specs.emplace_back(o2::event_visualisation::getPublisherSpec(PublisherConf{
                                                                   "tpc-raw-cluster-reader",
                                                                   "tpcraw",
                                                                   { "databranch", "TPCClusterHw", "Branch with TPC raw clusters" },
                                                                   { "mcbranch", "TPCClusterHwMCTruth", "MC label branch" },
                                                                   OutputSpec{ "TPC", "CLUSTERHW" },
                                                                   OutputSpec{ "TPC", "CLUSTERHWMCLBL" },
                                                                   tpcSectors,
                                                                   laneConfiguration,
                                                                 },
                                                                 propagateMC));
  } else if (inputType == InputType::Clusters) {
    specs.emplace_back(o2::event_visualisation::getPublisherSpec(PublisherConf{
                                                                   "tpc-native-cluster-reader",
                                                                   "tpcrec",
                                                                   { "clusterbranch", "TPCClusterNative", "Branch with TPC native clusters" },
                                                                   { "clustermcbranch", "TPCClusterNativeMCTruth", "MC label branch" },
                                                                   OutputSpec{ "TPC", "CLUSTERNATIVE" },
                                                                   OutputSpec{ "TPC", "CLNATIVEMCLBL" },
                                                                   tpcSectors,
                                                                   laneConfiguration,
                                                                 },
                                                                 propagateMC));
  }

  return std::move(specs);
}

} // end namespace reader_workflow
} // end namespace event_visualisation
} // end namespace o2