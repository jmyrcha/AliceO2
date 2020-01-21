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
/// \file    VisualisationConstants.h
/// \author  Jeremi Niedziela
/// \author  julian.myrcha@cern.ch
///

#ifndef ALICE_O2_EVENTVISUALISATION_BASE_VISUALISATIONCONSTANTS_H
#define ALICE_O2_EVENTVISUALISATION_BASE_VISUALISATIONCONSTANTS_H

#include <string>

namespace o2
{
namespace event_visualisation
{

enum EVisualisationGroup {
  ACO,
  EMC,
  HMP,
  MCH,
  PHS,
  RPH,
  TOF,
  SDD,
  SPD,
  SSD,
  ITS,
  TPC,
  TRD,
  AOD,
  NvisualisationGroups
};

const std::string gVisualisationGroupName[NvisualisationGroups] = {
  "ACO",
  "EMC",
  "HMP",
  "MCH",
  "PHS",
  "RPH",
  "TOF",
  "SDD",
  "SPD",
  "SSD",
  "ITS",
  "TPC",
  "TRD",
  "AOD"};

enum EVisualisationDataType {
  Raw,       ///< Raw data
  Hits,      ///< Hits
  Digits,    ///< Digits
  Clusters,  ///< Reconstructed clusters (RecPoints)
  Tracks,    ///< Reconstructed tracks
  Calo,      ///< Calorimeter cells (for AOD only)
  Muon,      ///< Muon tracks (for AOD only)
  NoData,    ///< no data was loaded
  NdataTypes ///< number of supported data types
};

const std::string gDataTypeNames[NdataTypes] = {
  "Raw",
  "Hits",
  "Digits",
  "Clusters",
  "Tracks",
  "Calo",
  "Muon",
  "NoData"};

enum EDataSource {
  SourceOnline,  ///< Online reconstruction is a source of events
  SourceOffline, ///< Local files are the source of events
  SourceHLT      ///< HLT reconstruction is a source of events
};

// FIXME: No reconstruction flag constants available yet??
// Remove the enum below once it is done
enum ERecoFlag {
  kITSin = 0x1,
  kITSout = 0x2,
  kITSrefit = 0x4,
  kITSpid = 0x8,
  kTPCin = 0x10,
  kTPCout = 0x20,
  kTPCrefit = 0x40,
  kTPCpid = 0x80,
  kTRDin = 0x100,
  kTRDout = 0x200,
  kTRDrefit = 0x400,
  kTRDpid = 0x800,
  kTOFin = 0x1000,
  kTOFout = 0x2000,
  kTOFrefit = 0x4000,
  kTOFpid = 0x8000,
  kHMPIDout = 0x10000,
  kHMPIDpid = 0x20000,
  kEMCALmatch = 0x40000,
  kTRDbackup = 0x80000,
  kTOFmismatch = 0x100000,
  kPHOSmatch = 0x200000,
  kITSupg = 0x400000,     // ITSupgrade reco
  kSkipFriend = 0x800000, // skip friend storage
  kGlobalMerge = 0x1000000,
  kMultInV0 = 0x2000000, // BIT(25): assumed to belong to V0 in multiplicity estimates
  kMultSec = 0x4000000,  // BIT(26): assumed to be secondary (due to the DCA) in multiplicity estimates
  kEmbedded = 0x8000000, // BIT(27), 1<<27: if it is a track that has been embedded into the event
  kITSpureSA = 0x10000000,
  kTRDStop = 0x20000000,
  kESDpid = 0x40000000,
  kTIME = 0x80000000
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_BASE_VISUALISATIONCONSTANTS_H
