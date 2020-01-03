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
/// \file    VisualisationCaloCell.h
/// \author  Maja Kabus
///

#ifndef ALICE_O2_EVENTVISUALISATION_DATACONVERTER_VISUALISATIONCALOCELL_H
#define ALICE_O2_EVENTVISUALISATION_DATACONVERTER_VISUALISATIONCALOCELL_H

#include "ConversionConstants.h"

#include <iosfwd>
#include <string>
#include <vector>
#include <array>
#include <cmath>

namespace o2
{
namespace event_visualisation
{

/// Minimalistic description of a calorimeter cell
///
/// This class is used mainly for visualisation purpose.
/// It keeps basic information about a calorimeter cell such as
/// its absolute ID, type (EMCAL vs PHOS), module, position, phi, eta and amplitude

class VisualisationCaloCell
{
 public:
  // Default constructor
  VisualisationCaloCell();

  // Constructor with properties initialisation
  VisualisationCaloCell(
    int absID,
    int module,
    double pos[],
    float amplitude,
    int type,
    float phi,
    float eta);

  // Cell absolute ID getter
  int getAbsID() { return mAbsID; }
  // Cell module getter
  int getModule() { return mModule; }
  // Cell amplitude getter
  float getAmplitude() { return mAmplitude; }
  // Cell type getter
  int getType() { return mType; }
  // Cell polar angle getter
  float getPhi() { return mPhi; }
  // Cell pseudorapidity getter
  float getEta() { return mEta; }
  // Cell position getters
  double X() { return mPos[0]; }
  double Y() { return mPos[1]; }
  double Z() { return mPos[2]; }

 private:
  int mAbsID;       /// Absolute cell ID
  int mModule;      /// Cell module
  double mPos[3];   /// Relative cell position in module
  float mAmplitude; /// Cell amplitude (= energy)
  int mType;        /// Cell type (-1 is undefined, 0 is PHOS, 1 is EMCAL)
  float mPhi;       /// Cell polar angle
  float mEta;       /// Cell pseudorapidity
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_DATACONVERTER_VISUALISATIONCALOCELL_H
