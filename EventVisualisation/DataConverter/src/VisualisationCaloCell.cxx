// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file    VisualisationCaloCell.cxx
/// \author  Maja Kabus
///

#include "EventVisualisationDataConverter/VisualisationCaloCell.h"

#include <TMath.h>

namespace o2
{
namespace event_visualisation
{

VisualisationCaloCell::VisualisationCaloCell() = default;

VisualisationCaloCell::VisualisationCaloCell(
  int absID,
  int module,
  double pos[],
  float amplitude,
  int type,
  float phi,
  float eta)
  : mAbsID(absID),
    mModule(module),
    mPos{pos[0], pos[1], pos[2]},
    mAmplitude(amplitude),
    mType(type),
    mPhi(phi),
    mEta(eta)
{
  // Clamp to [-pi, pi]
  if (mPhi > TMath::Pi()) {
    mPhi -= 2 * TMath::Pi();
  }
}

} // namespace event_visualisation
} // namespace o2
