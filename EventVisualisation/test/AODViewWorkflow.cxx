// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/WorkflowSpec.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/CompletionPolicy.h"
#include "Framework/DeviceSpec.h"
#include "AODReaderSpec.h"
#include "AODViewerSpec.h"

using namespace o2::framework;

// ------------------------------------------------------------------
#include "Framework/runDataProcessing.h"


/// This function is required to be implemented to define the workflow
/// specifications
WorkflowSpec defineDataProcessing(ConfigContext const& configcontext)
{
  // Reserve one entry which fill be filled with the SimReaderSpec
  // at the end. This places the processor at the beginning of the
  // workflow in the upper left corner of the GUI.
  WorkflowSpec specs;

  // put the reader
  specs.emplace_back(o2::viewer::getAODReaderSpec());

  // put the viewer
  specs.emplace_back(o2::viewer::getAODViewerSpec());

  return specs;
}