//
// Created by maja on 15.07.19.
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

#include "EventVisualisationWorkflow/TPCWorkflow.h"
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
namespace tpc_workflow
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

const std::unordered_map<std::string, OutputType> OutputMap{
  { "digits", OutputType::Digits },
  { "raw", OutputType::Raw },
  { "clusters", OutputType::Clusters },
  { "tracks", OutputType::Tracks },
};

framework::WorkflowSpec getWorkflow(std::vector<int> const& tpcSectors, std::vector<int> const& laneConfiguration, bool propagateMC, unsigned nLanes, std::string const& cfgInput, std::string const& cfgOutput)
{
  InputType inputType;

  try {
    inputType = InputMap.at(cfgInput);
  } catch (std::out_of_range&) {
    throw std::invalid_argument(std::string("invalid input type: ") + cfgInput);
  }
  std::vector<OutputType> outputTypes;
  try {
    outputTypes = RangeTokenizer::tokenize<OutputType>(cfgOutput, [](std::string const& token) { return OutputMap.at(token); });
  } catch (std::out_of_range&) {
    throw std::invalid_argument(std::string("invalid output type: ") + cfgOutput);
  }
  auto isEnabled = [&outputTypes](OutputType type) {
    return std::find(outputTypes.begin(), outputTypes.end(), type) != outputTypes.end();
  };

  if (inputType == InputType::Raw && isEnabled(OutputType::Digits)) {
    throw std::invalid_argument("input/output type mismatch, can not produce 'digits' from 'raw'");
  }
  if (inputType == InputType::Clusters && (isEnabled(OutputType::Digits) || isEnabled(OutputType::Raw))) {
    throw std::invalid_argument("input/output type mismatch, can not produce 'digits', nor 'raw' from 'clusters'");
  }

  WorkflowSpec specs;
  // Possibly checking input arguments

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

  //////////////////////////////////////////////////////////////////////////////////////////////
  //
  // generation of processor specs for various types of outputs
  // based on generic RootTreeWriter and MakeRootTreeWriterSpec generator
  //
  // -------------------------------------------------------------------------------------------
  // the callbacks for the RootTreeWriter
  //
  // The generic writer needs a way to associate incoming data with the individual branches for
  // the TPC sectors. The sector number is transmitted as part of the sector header, the callback
  // finds the corresponding index in the vector of configured sectors
  auto getIndex = [tpcSectors](o2::framework::DataRef const& ref) {
    auto const* tpcSectorHeader = o2::framework::DataRefUtils::getHeader<o2::tpc::TPCSectorHeader*>(ref);
    if (!tpcSectorHeader) {
      throw std::runtime_error("TPC sector header missing in header stack");
    }
    if (tpcSectorHeader->sector < 0) {
      // special data sets, don't write
      return ~(size_t)0;
    }
    size_t index = 0;
    for (auto const& sector : tpcSectors) {
      if (sector == tpcSectorHeader->sector) {
        return index;
      }
      ++index;
    }
    throw std::runtime_error("sector " + std::to_string(tpcSectorHeader->sector) + " not configured for writing");
  };
  auto getName = [tpcSectors](std::string base, size_t index) {
    return base + "_" + std::to_string(tpcSectors.at(index));
  };

  auto checkReady = [](o2::framework::DataRef const& ref) {
    auto const* tpcSectorHeader = o2::framework::DataRefUtils::getHeader<o2::tpc::TPCSectorHeader*>(ref);
    // sector number -1 indicates end-of-data
    if (tpcSectorHeader != nullptr) {
      if (tpcSectorHeader->sector == -1) {
        // indicate normal processing if not ready and skip if ready
        return std::make_tuple(MakeRootTreeWriterSpec::TerminationCondition::Action::SkipProcessing, true);
      }
    }
    return std::make_tuple(MakeRootTreeWriterSpec::TerminationCondition::Action::DoProcessing, false);
  };

  // -------------------------------------------------------------------------------------------
  // helper to create writer specs for different types of output
  auto makeWriterSpec = [tpcSectors, laneConfiguration, propagateMC, getIndex, getName, checkReady](const char* processName, const char* defaultFileName, const char* defaultTreeName, auto&& databranch, auto&& mcbranch) {
    if (tpcSectors.size() == 0) {
      throw std::invalid_argument(std::string("writer process configuration needs list of TPC sectors"));
    }

    auto amendInput = [tpcSectors, laneConfiguration](InputSpec& input, size_t index) {
      input.binding += std::to_string(laneConfiguration[index]);
      DataSpecUtils::updateMatchingSubspec(input, laneConfiguration[index]);
    };
    auto amendBranchDef = [laneConfiguration, propagateMC, amendInput, tpcSectors, getIndex, getName](auto&& def) {
      def.keys = mergeInputs(def.keys, laneConfiguration.size(), amendInput);
      def.nofBranches = tpcSectors.size();
      def.getIndex = getIndex;
      def.getName = getName;
      return std::move(def);
    };

    // depending on the MC propagation flag, the RootTreeWriter spec is created with two
    // or one branch definition
    if (propagateMC) {
      return std::move(MakeRootTreeWriterSpec(processName, defaultFileName, defaultTreeName,
                                              MakeRootTreeWriterSpec::TerminationCondition{ checkReady },
                                              std::move(amendBranchDef(databranch)),
                                              std::move(amendBranchDef(mcbranch)))());
    }
    return std::move(MakeRootTreeWriterSpec(processName, defaultFileName, defaultTreeName,
                                            MakeRootTreeWriterSpec::TerminationCondition{ checkReady },
                                            std::move(amendBranchDef(databranch)))());
  };

  //////////////////////////////////////////////////////////////////////////////////////////////
  //
  // a writer process for digits
  //
  // selected by output type 'digits'
  if (isEnabled(OutputType::Digits)) {
    using DigitOutputType = std::vector<o2::tpc::Digit>;
    using MCLabelContainer = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
    specs.push_back(makeWriterSpec("tpc-digits-writer",
                                   inputType == InputType::Digits ? "tpc-filtered-digits.root" : "tpcdigits.root",
                                   "o2sim",
                                   BranchDefinition<DigitOutputType>{ InputSpec{ "data", "TPC", "DIGITS", 0 },
                                                                      "TPCDigit", "digit-branch-name" },
                                   BranchDefinition<MCLabelContainer>{ InputSpec{ "mc", "TPC", "DIGITSMCTR", 0 },
                                                                       "TPCDigitMCTruth", "digitmc-branch-name" }));
  }

  //////////////////////////////////////////////////////////////////////////////////////////////
  //
  // a writer process for raw hardware clusters
  //
  // selected by output type 'raw'
  if (isEnabled(OutputType::Raw)) {
    using MCLabelContainer = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
    specs.push_back(makeWriterSpec("tpc-raw-cluster-writer",
                                   inputType == InputType::Raw ? "tpc-filtered-raw-clusters.root" : "tpc-raw-clusters.root",
                                   "tpcraw",
                                   BranchDefinition<const char*>{ InputSpec{ "data", "TPC", "CLUSTERHW", 0 },
                                                                  "TPCClusterHw", "databranch" },
                                   BranchDefinition<MCLabelContainer>{ InputSpec{ "mc", "TPC", "CLUSTERHWMCLBL", 0 },
                                                                       "TPCClusterHwMCTruth", "mcbranch" }));
  }

  //////////////////////////////////////////////////////////////////////////////////////////////
  //
  // a writer process for TPC native clusters
  //
  // selected by output type 'clusters'
  if (isEnabled(OutputType::Clusters)) {
    using MCLabelCollection = std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>>;
    specs.push_back(makeWriterSpec("tpc-native-cluster-writer",
                                   inputType == InputType::Clusters ? "tpc-filtered-native-clusters.root" : "tpc-native-clusters.root",
                                   "tpcrec",
                                   BranchDefinition<const char*>{ InputSpec{ "data", "TPC", "CLUSTERNATIVE", 0 },
                                                                  "TPCClusterNative", "databranch" },
                                   BranchDefinition<MCLabelCollection>{ InputSpec{ "mc", "TPC", "CLNATIVEMCLBL", 0 },
                                                                        "TPCClusterNativeMCTruth", "mcbranch" }));
  }

  return std::move(specs);
}

} // end namespace tpc_workflow
} // end namespace event_visualisation
} // end namespace o2
