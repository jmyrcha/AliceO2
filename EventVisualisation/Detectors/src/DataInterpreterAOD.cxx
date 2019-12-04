// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DataInterpreterAOD.cxx
/// \brief converting AOD data to Event Visualisation primitives
/// \author Maja Kabus

#include "EventVisualisationDetectors/DataInterpreterAOD.h"
#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationDataConverter/VisualisationEvent.h"
#include "EventVisualisationDetectors/CaloMatrix.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/PID.h"
#include "EMCALBase/Geometry.h"
#include "PHOSBase/Geometry.h"
#include "DataFormatsEMCAL/Constants.h"
#include "CCDB/Manager.h"
#include "CCDB/IdPath.h"
#include "CCDB/Condition.h"

#include "FairLogger.h"

#include <TEveManager.h>
#include <TEveTrackPropagator.h>
#include <TEveTrack.h>
#include <TGListTree.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector2.h>
#include <TEveQuadSet.h>
#include <TGeoNode.h>
#include <TGeoManager.h>
#include <TClonesArray.h>

#include <gsl/span>

using namespace o2::ccdb;

namespace o2
{
namespace event_visualisation
{

DataInterpreterAOD::DataInterpreterAOD()
{
  // TODO: This is a Run 2 geometry. Change to Run 3 once available
  mEMCALGeom = o2::emcal::Geometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM");
  mPHOSGeom = o2::phos::Geometry::GetInstance("Run2");

  // Version 1: EMCAL alignment matrices from O2CDB
  // Use o2-fill-o2cdb-emcal-phos-matrix to put ideal alignment matrices into O2CDB.
  //  LOG(INFO) << "Loading matrices from OCDB";
  //  IdPath path("EMCAL", "Align", "Data");
  //  Condition* cond = o2::ccdb::Manager::Instance()->getCondition(path, 1);
  //  if (!cond) {
  //    LOG(FATAL) << "Couldn't load EMCAL alignment matrices from CDB!";
  //  }
  //  cond->setOwner(0);
  //  TClonesArray* array = (TClonesArray*)cond->getObject();
  //  TClonesArray& alobj = *array;
  //
  //  for (int mod = 0; mod < mEMCALGeom->GetNumberOfSuperModules(); mod++) {
  //    if (!mEMCALGeom->GetMatrixForSuperModuleFromArray(mod)) {
  //      TGeoHMatrix* matrix = (TGeoHMatrix*)alobj[mod];
  //      mEMCALGeom->SetMisalMatrix(matrix, mod);
  //    }
  //  }

  // Version 2: EMCAL hardcoded matrices from sample ESD (for visual comparison)
  LOG(INFO) << "Loading hardcoded matrices";
  for (int mod = 0; mod < mEMCALGeom->GetNumberOfSuperModules(); mod++) {
    if (!mEMCALGeom->GetMatrixForSuperModuleFromArray(mod)) {
      TGeoHMatrix* matrix = CaloMatrix::getEMCALMatrix(mod);
      mEMCALGeom->SetMisalMatrix(matrix, mod);
    }
  }

  // TODO: Setting PHOS matrices once it will be possible
  // Currently no setter and getter for PHOS matrices in PHOS geometry
  //  path = IdPath("PHOS", "Align", "Data");
  //  cond = o2::ccdb::Manager::Instance()->getCondition(path, 1);
  //  if (!cond) {
  //    LOG(FATAL) << "Couldn't load PHOS alignment matrices from CDB!";
  //  }
  //  cond->setOwner(0);
  //  array->Clear();
  //  array = (TClonesArray*)cond->getObject();
  //  alobj = *array;
  //
  //  for (int mod = 0; mod < 5; mod++) {
  //    if (!mPHOSGeom->GetMatrixForSuperModuleFromArray(mod)) {
  //      TGeoHMatrix* matrix = (TGeoHMatrix*) alobj[mod];
  //      mPHOSGeom->SetMisalMatrix(matrix, mod);
  //    }
  //  }
}

DataInterpreterAOD::~DataInterpreterAOD() = default;

void DataInterpreterAOD::interpretDataForType(TObject* data, EVisualisationDataType type, VisualisationEvent& event)
{
  TList* list = (TList*)data;
  Int_t eventId = ((TVector2*)list->At(1))->X();
  TFile* AODFile = (TFile*)list->At(0);

  if (type == Tracks) {
    interpretAODTracks(AODFile, eventId, event);
  } else if (type == Calo) {
    interpretAODCaloCells(AODFile, eventId, event);
  } else if (type == Muon) {
    interpretMuonTracks(AODFile, eventId, event);
  }
}

void DataInterpreterAOD::interpretAODTracks(TFile* AODFile, Int_t eventId, VisualisationEvent& event)
{
  TTree* tracks = (TTree*)AODFile->Get("O2tracks");

  // Read all tracks parameters to buffers
  Int_t trkID;          // The index of the collision vertex, to which the track is attached
  Float_t trkX;         // X coordinate for the point of parametrisation
  Float_t trkAlpha;     // Local <--> global coor.system rotation angle
  Float_t trkY;         // fP[0] local Y-coordinate of a track (cm)
  Float_t trkZ;         // fP[1] local Z-coordinate of a track (cm)
  Float_t trkSnp;       // fP[2] local sine of the track momentum azimuthal angle
  Float_t trkTgl;       // fP[3] tangent of the track momentum dip angle
  Float_t trkSigned1Pt; // fP[4] 1/pt (1/(GeV/c))
  o2::track::PID trkPID;
  ULong64_t trkFlags; // Reconstruction status flags - track quality parameters
  tracks->SetBranchAddress("fID4Tracks", &trkID);
  tracks->SetBranchAddress("fX", &trkX);
  tracks->SetBranchAddress("fAlpha", &trkAlpha);
  tracks->SetBranchAddress("fY", &trkY);
  tracks->SetBranchAddress("fZ", &trkZ);
  tracks->SetBranchAddress("fSnp", &trkSnp);
  tracks->SetBranchAddress("fTgl", &trkTgl);
  tracks->SetBranchAddress("fSigned1Pt", &trkSigned1Pt);
  tracks->SetBranchAddress("fPID", &trkPID);
  tracks->SetBranchAddress("fFlags", &trkFlags);

  TEveTrackList* trackList = new TEveTrackList("tracks");
  trackList->IncDenyDestroy();
  auto prop = trackList->GetPropagator();
  prop->SetMagField(0.5);

  int tracksCount = tracks->GetEntriesFast();
  for (int i = 0; i < tracksCount; i++) {
    tracks->GetEntry(i);

    // TODO: This is provisional.
    //  AOD is supposed to contain references from primary vertex to elements from given collision.
    if (trkID > eventId + 10)
      break; // End event search
    if (trkID != eventId)
      continue;

    o2::track::TrackPar rec(trkX, trkAlpha, { trkY, trkZ, trkSnp, trkTgl, trkSigned1Pt });
    std::array<float, 3> p;
    rec.getPxPyPzGlo(p);

    TEveRecTrackD t;
    t.fP = { p[0], p[1], p[2] };
    t.fSign = (rec.getSign() < 0) ? -1 : 1;
    TEveTrack* eve_track = new TEveTrack(&t, prop);
    eve_track->MakeTrack();

    auto start = eve_track->GetLineStart();
    auto end = eve_track->GetLineEnd();
    double track_start[3] = { start.fX, start.fY, start.fZ };
    double track_end[3] = { end.fX, end.fY, end.fZ };
    double track_p[3] = { p[0], p[1], p[2] };

    VisualisationTrack track(rec.getSign(), 0.0, trkID, trkPID, 0.0, trkSigned1Pt, track_start, track_end, track_p, 0, 0.0, 0.0, 0.0, 0, trkFlags);

    for (Int_t i = 0; i < eve_track->GetN(); ++i) {
      Float_t x, y, z;
      eve_track->GetPoint(i, x, y, z);
      track.addPolyPoint(x, y, z);
    }
    delete eve_track;

    event.addTrack(track);
  }
  delete trackList;
}

void DataInterpreterAOD::interpretAODCaloCells(TFile* AODFile, Int_t eventId, VisualisationEvent& event)
{
  TTree* calo = (TTree*)AODFile->Get("O2calo");
  Int_t caloID;           // The index of the collision vertex, to which the cell is attached
  Short_t caloCellNumber; // Cell absolute Id. number
  Float_t caloAmplitude;  // Cell amplitude (= energy!)
  Char_t caloType;        // Cell type (-1 is undefined, 0 is PHOS, 1 is EMCAL)
  calo->SetBranchAddress("fID4Calo", &caloID);
  calo->SetBranchAddress("fCellNumber", &caloCellNumber);
  calo->SetBranchAddress("fAmplitude", &caloAmplitude);
  calo->SetBranchAddress("fType", &caloType);

  int caloCount = calo->GetEntriesFast();

  // TODO: In old alieve it was calculated for a given cluster before cluster cuts. What in O2 instead?
  float maxCellEnergy = 0.0;
  int maxCellEnergyAbsId = -1;
  // Storing values in vectors so as the file will not have to be traversed multiple times
  std::vector<Int_t> caloAbsIds = std::vector<Int_t>();
  std::vector<Float_t> caloAmplitudes = std::vector<Float_t>();
  std::vector<Char_t> caloTypes = std::vector<Char_t>();

  for (int i = 0; i < caloCount; i++) {
    calo->GetEntry(i);

    // TODO: This is provisional.
    //  AOD is supposed to contain references from primary vertex to elements from given collision.
    if (caloID > eventId + 10)
      break; // End event search
    if (caloID != eventId)
      continue;

    if (caloAmplitude > maxCellEnergy) {
      maxCellEnergy = caloAmplitude;
      maxCellEnergyAbsId = caloCellNumber;
    }

    caloAbsIds.emplace_back(caloCellNumber);
    caloAmplitudes.emplace_back(caloAmplitude);
    caloTypes.emplace_back(caloType);
  }

  caloCount = caloAbsIds.size();
  for (int i = 0; i < caloCount; i++) {
    // Cell rejected
    if (cutCell(caloAmplitudes, caloAbsIds, caloCount, caloAmplitudes[i], caloTypes[i],
                maxCellEnergy, maxCellEnergyAbsId))
      continue;

    if (caloTypes[i] == 0) {
      interpretPHOSCell(caloAbsIds[i], caloAmplitudes[i], event);
    } else if (caloTypes[i] == 1) {
      interpretEMCALCell(caloAbsIds[i], caloAmplitudes[i], event);
    } else {
      LOG(FATAL) << "Wrong calorimeter type";
    }
  }
}

void DataInterpreterAOD::interpretEMCALCell(Int_t absId, Float_t amplitude, VisualisationEvent& event)
{
  auto [iSupMod, iTower, iIphi, iIeta] = mEMCALGeom->GetCellIndex(absId);
  // It should not happen, but in case the OCDB file is not the correct one.
  if (iSupMod >= mEMCALGeom->GetNumberOfSuperModules()) {
    LOG(WARNING) << "iSupMod: " << iSupMod << " bigger than EMCAL modules number: "
                 << mEMCALGeom->GetNumberOfSuperModules();
    return;
  }

  // Gives label of cell in eta-phi position per each supermodule
  //auto [iphi, ieta] = mEMCALGeom->GetCellPhiEtaIndexInSModule(iSupMod, iTower, iIphi, iIeta);
  auto [eta, phi] = mEMCALGeom->EtaPhiFromIndex(absId);

  const auto& pos = mEMCALGeom->RelPosCellInSModule(absId);
  double posArr[3] = { pos.X(), pos.Y(), pos.Z() };

  VisualisationCaloCell cell(absId, iSupMod, posArr, amplitude, 1, phi, eta);

  event.addCaloCell(cell);
}

void DataInterpreterAOD::interpretPHOSCell(Int_t absId, Float_t amplitude, VisualisationEvent& event)
{
  TVector3 xyz;
  Int_t relId[3];
  mPHOSGeom->AbsToRelNumbering(absId, relId);
  double x, z;
  mPHOSGeom->AbsIdToRelPosInModule(absId, x, z);
  double posArr[3] = { x, 0.0, z };

  // TODO: No Local2Global in PHOS geometry yet. Uncomment it once available.
  // mPHOSGeom->Local2Global(relId[0], x, z, xyz);
  // eta = xyz.Eta(), phi = xyz.Phi()
  // fHistoPH->Fill(xyz.Eta(),GetPhi(xyz.Phi()),amp);

  VisualisationCaloCell cell(absId, relId[0] - 1, posArr, amplitude, 0, relId[1], 0);

  //  printf("\t PHOS cell: ID %d, energy %2.2f,eta %2.2f, phi %2.2f, phi not modified: %2.2f\n",
  //         absId, amplitude, eta, cell.getPhi()*TMath::RadToDeg(), phi);
  //  printf("\t SM %d iEta %d,  iPhi %d, x %3.3f, y %3.3f, z %3.3f \n",
  //         iSupMod, ieta, iphi, pos.X(), pos.Y(), pos.Z());

  event.addCaloCell(cell);
}

// TODO: Provisional cuts for all cells - in AliRoot cuts were applied for each cluster separately.
Bool_t DataInterpreterAOD::cutCell(std::vector<Float_t> caloAmplitudes, std::vector<Int_t> caloAbsIds, int caloCount, Float_t amplitude, Char_t caloType, Float_t maxCellEnergy, Int_t maxCellEnergyAbsId)
{
  //  if(caloCount < mNumMinCellsCut[caloType]) return kTRUE;
  //  if(caloCount > mNumMaxCellsCut[caloType]) return kTRUE;

  if (amplitude < mEnergyCut[caloType] / 2.0)
    return kTRUE;

  //  if(fCaloCluster->GetM02()    < m02lowCut[caloType]) return kTRUE ;
  //  if(fCaloCluster->GetM02()    > m02higCut[caloType]) return kTRUE ;
  //  if(fCaloCluster->GetM20()    < m20lowCut[caloType]) return kTRUE ;

  //  if (maxCellEnergy < mEnergyCut[caloType] / 2.0 || maxCellEnergyAbsId < 0) {
  //    return kTRUE;
  //  }

  //    auto [ism, irow, icol] = getModuleNumberColAndRow(amplitude, caloType);
  //    Float_t eCross = getECross(caloAmplitudes, caloAbsIds, caloType, ism, icol, irow);
  //    if (1 - eCross / maxCellEnergy > mExoCut) {
  //      return kTRUE;
  //    }

  return kFALSE;
}

/// Adapted from emcal_esdclustercells.C in AliRoot
/// Get for a given cell absolute ID the (super)module number, column and row
///
/// \param absId absolute cell identity number
/// \param type 1 for EMCAL, 0 for PHOS
/// \return tuple with reference cell module, row and column
std::tuple<Int_t, Int_t, Int_t> DataInterpreterAOD::getModuleNumberColAndRow(Int_t absId, Char_t caloType)
{
  Int_t iSM = -1, irow = -1, icol = -1;
  if (absId < 0)
    return std::make_tuple(-1, -1, -1);

  if (caloType == 1) {
    auto [iSupMod, iTower, iIphi, iIeta] = mEMCALGeom->GetCellIndex(absId);
    if (iSupMod >= mEMCALGeom->GetNumberOfSuperModules())
      return std::make_tuple(-1, -1, -1);
    std::tie(irow, icol) = mEMCALGeom->GetCellPhiEtaIndexInSModule(iSupMod, iTower, iIphi, iIeta);
  }    // EMCAL
  else // PHOS
  {
    Int_t relId[4];
    mPHOSGeom->AbsToRelNumbering(absId, relId);
    irow = relId[2];
    icol = relId[3];
    iSM = relId[0] - 1;
  } // PHOS
  return std::make_tuple(iSM, irow, icol);
}

/// Adapted from emcal_esdclustercells.C in AliRoot
/// \return sum of energy in neighbor cells (cross) to the reference cell
///
/// \param isEMCAL bool true for EMCAL, false for PHOS
/// \param imod reference cell module
/// \param icol reference cell column
/// \param irow reference cell row
///
Float_t DataInterpreterAOD::getECross(std::vector<Float_t>& caloAmplitudes, std::vector<Int_t>& caloAbsIds,
                                      Char_t caloType, Int_t imod, Int_t icol, Int_t irow)
{
  Float_t ecell1 = 0, ecell2 = 0, ecell3 = 0, ecell4 = 0;
  Int_t absId1 = -1, absId2 = -1, absId3 = -1, absId4 = -1;

  if (caloType == 1) {
    // Get close cells index, energy and time, not in corners
    Int_t rowMax = mEMCALGeom->GetNPhi() * 2; // AliEMCALGeoParams::fgkEMCALRows;
    Int_t colMax = mEMCALGeom->GetNEta() * 2; // AliEMCALGeoParams::fgkEMCALCols;

    if (imod == 11 || imod == 10 || imod == 18 || imod == 19)
      rowMax /= 3;
    // rowMax = AliEMCALGeoParams::fgkEMCALRows/3;

    if (imod > 12 && imod < 18)
      colMax = colMax * 2 / 3;
    // colMax = AliEMCALGeoParams::fgkEMCALCols*2/3;

    if (irow < rowMax)
      absId1 = mEMCALGeom->GetAbsCellIdFromCellIndexes(imod, irow + 1, icol);
    if (irow > 0)
      absId2 = mEMCALGeom->GetAbsCellIdFromCellIndexes(imod, irow - 1, icol);

    if (icol < colMax - 1)
      absId3 = mEMCALGeom->GetAbsCellIdFromCellIndexes(imod, irow, icol + 1);
    if (icol > 0)
      absId4 = mEMCALGeom->GetAbsCellIdFromCellIndexes(imod, irow, icol - 1);

    // In case of cell in eta = 0 border, depending on SM shift the cross cell index
    if (imod > 11 && imod < 18) {
      if (icol == colMax - 1 && !(imod % 2)) {
        absId3 = mEMCALGeom->GetAbsCellIdFromCellIndexes(imod + 1, irow, 0);
        absId4 = mEMCALGeom->GetAbsCellIdFromCellIndexes(imod, irow, icol - 1);
      } else if (icol == 0 && imod % 2) {
        absId3 = mEMCALGeom->GetAbsCellIdFromCellIndexes(imod, irow, icol + 1);
        absId4 = mEMCALGeom->GetAbsCellIdFromCellIndexes(imod - 1, irow, colMax - 1);
      }
    }
  }
  //  else // PHOS
  //  {
  //    Int_t relId1[] = { imod+1, 0, irow+1, icol   };
  //    Int_t relId2[] = { imod+1, 0, irow-1, icol   };
  //    Int_t relId3[] = { imod+1, 0, irow  , icol+1 };
  //    Int_t relId4[] = { imod+1, 0, irow  , icol-1 };
  //
  //    // TODO: RelToAbsNumbering() not implemented yet in PHOS!
  //    mPHOSGeom->RelToAbsNumbering(relId1, absId1);
  //    mPHOSGeom->RelToAbsNumbering(relId2, absId2);
  //    mPHOSGeom->RelToAbsNumbering(relId3, absId3);
  //    mPHOSGeom->RelToAbsNumbering(relId4, absId4);
  //  }

  if (absId1 > 0)
    ecell1 = getCaloCellAmplitude(caloAmplitudes, caloAbsIds, absId1);
  if (absId2 > 0)
    ecell2 = getCaloCellAmplitude(caloAmplitudes, caloAbsIds, absId2);
  if (absId3 > 0)
    ecell3 = getCaloCellAmplitude(caloAmplitudes, caloAbsIds, absId3);
  if (absId4 > 0)
    ecell4 = getCaloCellAmplitude(caloAmplitudes, caloAbsIds, absId4);

  return ecell1 + ecell2 + ecell3 + ecell4;
}

Float_t DataInterpreterAOD::getCaloCellAmplitude(std::vector<Float_t> caloAmplitudes,
                                                 std::vector<Int_t> caloAbsIds, Int_t absId)
{
  auto it = std::find_if(caloAbsIds.begin(), caloAbsIds.end(), [absId](int id) { return id == absId; });
  if (it == caloAbsIds.end()) {
    LOG(WARNING) << "getCaloCellAmplitude: No cell of given abs ID!" << std::endl;
    return 0.0;
  }
  int ind = it - caloAbsIds.begin();
  return caloAmplitudes[ind];
}

void DataInterpreterAOD::interpretMuonTracks(TFile* AODFile, Int_t eventId, VisualisationEvent& event)
{
  TTree* muon = (TTree*)AODFile->Get("O2mu");

  Int_t muonID;                       // The index of the collision vertex, to which the muon is attached
  Float_t muonInverseBendingMomentum; // Inverse bending momentum (GeV/c ** -1) times the charge
  Float_t muonThetaX;                 // Angle of track at vertex in X direction (rad)
  Float_t muonThetaY;                 // Angle of track at vertex in Y direction (rad)
  Float_t muonZ;                      // Z coordinate (cm)
  Float_t muonBendingCoor;            // bending coordinate (cm)
  Float_t muonNonBendingCoor;         // non bending coordinate (cm)
  Float_t muonCovariances[15];        // reduced covariance matrix of UNCORRECTED track parameters AT FIRST CHAMBER
  Float_t muonChi2;                   // chi2 in the MUON track fit
  Float_t muonChi2MatchTrigger;       // chi2 of trigger/track matching
  muon->SetBranchAddress("fID4mu", &muonID);
  muon->SetBranchAddress("fInverseBendingMomentum", &muonInverseBendingMomentum);
  muon->SetBranchAddress("fThetaX", &muonThetaX);
  muon->SetBranchAddress("fThetaY", &muonThetaY);
  muon->SetBranchAddress("fZ", &muonZ);
  muon->SetBranchAddress("fBendingCoor", &muonBendingCoor);
  muon->SetBranchAddress("fNonBendingCoor", &muonNonBendingCoor);
  muon->SetBranchAddress("fCovariances", &muonCovariances);
  muon->SetBranchAddress("fChi2", &muonChi2);
  muon->SetBranchAddress("fChi2MatchTrigger", &muonChi2MatchTrigger);

  // load trigger circuit
  //    static AliMUONTriggerCircuit* gTriggerCircuit = 0x0;
  //    if (!gTriggerCircuit)
  //    {
  //      AliEveEventManager::Instance()->AssertGeometry();
  //      AliMUONGeometryTransformer* fMUONGeometryTransformer = new AliMUONGeometryTransformer();
  //      fMUONGeometryTransformer->LoadGeometryData();
  //      gTriggerCircuit = new AliMUONTriggerCircuit(fMUONGeometryTransformer);
  //    }

  TEveRecTrack recTrack;
  TEveTrack* track;
  //TEveTrackList* matchedList;
  //TEveTrackList* unmatchedList;
  TEveTrackList* ghostList = new TEveTrackList("ghost");
  ghostList->IncDenyDestroy();
  auto prop = ghostList->GetPropagator();
  prop->SetMagField(0.5);

  int muonCount = muon->GetEntriesFast();
  for (int i = 0; i < muonCount; i++) {
    muon->GetEntry(i);

    // TODO: This is provisional.
    //  AOD is supposed to contain references from primary vertex to elements from given collision.
    if (muonID > eventId + 10)
      break; // End event search
    if (muonID != eventId)
      continue;

    //recTrack.fIndex = muonID;
    //recTrack.fLabel = emt->GetLabel();
    //recTrack.fV.Set(muonNonBendingCoor, muonBendingCoor, muonZ);
    //      recTrack.fStatus = emt->GetMatchTrigger();
    //      recTrack.fSign = emt->Charge();
    //      recTrack.fP.Set(emt->PxAtDCA(),emt->PyAtDCA(),emt->PzAtDCA());
    //      recTrack.fBeta = ( emt->E() > 0 ) ? emt->P()/emt->E() : 0;

    //track = new TEveTrack(&recTrack, matched->GetPropagator());
    //track->SetName(Form("%cmu", emt->Charge()>0 ? '+':'-'));
    //track->SetStdTitle();
    //track->SetSourceObject(emt); // WARNING: Change the UniqueID of the object!!

    // TODO: Recognize different kinds of muons - matched / unmatched / ghost
    // TODO: Chi test for tracker??

    // Lines from ghost muon tracks
    recTrack.fStatus = 0;
    //recTrack.fSign = emt->Charge();
    //Double_t z11 = (muonZ < -1.0) ? muonZ : (Double_t)AliMUONConstants::DefaultChamberZ(10);
    recTrack.fV.Set(muonNonBendingCoor, muonBendingCoor, muonZ);
    recTrack.fP.Set(-TMath::Tan(muonThetaX), -TMath::Tan(muonThetaY), -1.0);

    // produce eve track
    track = new TEveTrack(&recTrack, prop);
    track->SetName("mu");
    track->SetTitle("Trigger only");
    //track->SetSourceObject(emt);

    // add the track to proper list
    //      track->SetAttLineAttMarker(ghostList);
    //      ghost->AddElement(track);

    // TODO: How to get proper track coordinates??
    track->MakeTrack();
    auto start = track->GetLineStart();
    auto end = track->GetLineEnd();
    double track_start[3] = { start.fX, start.fY, start.fZ };
    double track_end[3] = { end.fX, end.fY, end.fZ };
    double track_p[3] = { 0.0, 0.0, 0.0 }; // { p[0], p[1], p[2] };

    VisualisationTrack visTrack(recTrack.fSign, 0.0, 0, o2::track::PID::Muon, 0.0, 0.0, track_start, track_end, track_p, 0, 0.0, 0.0, 0.0, 0, 0);

    for (Int_t i = 0; i < track->GetN(); ++i) {
      Float_t x, y, z;
      track->GetPoint(i, x, y, z);
      visTrack.addPolyPoint(x, y, z);
    }

    delete track;
    event.addMuonTrack(visTrack);
  }
  delete ghostList;
}

} // namespace event_visualisation
} // namespace o2
