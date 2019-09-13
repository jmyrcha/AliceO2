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

#include <iostream>
#include <gsl/span>

using namespace std;
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
  IdPath path("EMCAL", "Align", "Data");
  Condition* cond = o2::ccdb::Manager::Instance()->getCondition(path, 1);
  if (!cond) {
    LOG(FATAL) << "Couldn't load EMCAL alignment matrices from CDB!";
  }
  cond->setOwner(0);
  TClonesArray* array = (TClonesArray*)cond->getObject();
  TClonesArray& alobj = *array;

  for (int mod = 0; mod < mEMCALGeom->GetNumberOfSuperModules(); mod++) {
    if (!mEMCALGeom->GetMatrixForSuperModuleFromArray(mod)) {
      TGeoHMatrix* matrix = (TGeoHMatrix*)alobj[mod];
      mEMCALGeom->SetMisalMatrix(matrix, mod);
    }
  }

  // Version 2: EMCAL hardcoded matrices from sample ESD (for visual comparison)
  //for (int mod = 0; mod < mEMCALGeom->GetNumberOfSuperModules(); mod++) {
  //  if (!mEMCALGeom->GetMatrixForSuperModuleFromArray(mod)) {
  //    TGeoHMatrix* matrix = CaloMatrix::getEMCALMatrix(mod);
  //    mEMCALGeom->SetMisalMatrix(matrix, mod);
  //  }
  //}

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

std::unique_ptr<VisualisationEvent> DataInterpreterAOD::interpretDataForType(TObject* data, EVisualisationDataType type)
{
  TList* list = (TList*)data;
  Int_t event = ((TVector2*)list->At(1))->X();
  TFile* AODFile = (TFile*)list->At(0);

  if (type == Tracks) {
    return this->interpretAODTracks(AODFile, event);
  } else if (type == Calo) {
    return this->interpretAODCaloCells(AODFile, event);
  } else if (type == Muon) {
    return this->interpretMuonTracks(AODFile, event);
  }
  return std::make_unique<VisualisationEvent>(0, 0, 0, 0, "", 0);
}

std::unique_ptr<VisualisationEvent> DataInterpreterAOD::interpretAODTracks(TFile* AODFile, Int_t event)
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
  tracks->SetBranchAddress("fID4Tracks", &trkID);
  tracks->SetBranchAddress("fX", &trkX);
  tracks->SetBranchAddress("fAlpha", &trkAlpha);
  tracks->SetBranchAddress("fY", &trkY);
  tracks->SetBranchAddress("fZ", &trkZ);
  tracks->SetBranchAddress("fSnp", &trkSnp);
  tracks->SetBranchAddress("fTgl", &trkTgl);
  tracks->SetBranchAddress("fSigned1Pt", &trkSigned1Pt);

  auto ret_event = std::make_unique<VisualisationEvent>(0, 0, 0, 0, "", 0);
  TEveTrackList* trackList = new TEveTrackList("tracks");
  trackList->IncDenyDestroy();
  auto prop = trackList->GetPropagator();
  prop->SetMagField(0.5);

  int tracksCount = tracks->GetEntriesFast();

  for (int i = 0; i < tracksCount; i++) {
    tracks->GetEntry(i);

    // TODO: This is provisional.
    //  AOD is supposed to contain references from primary vertex to elements from given collision.
    if (trkID > event + 10)
      break; // End event search
    if (trkID != event)
      continue;

    o2::track::TrackPar rec(trkX, trkAlpha, { trkY, trkZ, trkSnp, trkTgl, trkSigned1Pt });
    std::array<float, 3> p;
    // TODO: Convert coordinates local --> global - correct?
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

    VisualisationTrack track(rec.getSign(), 0.0, trkID, 0, 0.0, trkSigned1Pt, track_start, track_end, track_p, 0, 0.0, 0.0, 0.0, 0);

    for (Int_t i = 0; i < eve_track->GetN(); ++i) {
      Float_t x, y, z;
      eve_track->GetPoint(i, x, y, z);
      track.addPolyPoint(x, y, z);
    }
    delete eve_track;

    ret_event->addTrack(track);
  }
  delete trackList;
  return ret_event;
}

std::unique_ptr<VisualisationEvent> DataInterpreterAOD::interpretAODCaloCells(TFile* AODFile, Int_t eventID)
{
  auto ret_event = std::make_unique<VisualisationEvent>(0, 0, 0, 0, "", 0);

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
  for (int i = 0; i < caloCount; i++) {
    calo->GetEntry(i);

    // TODO: This is provisional.
    //  AOD is supposed to contain references from primary vertex to elements from given collision.
    if (caloID > eventID + 10)
      break; // End event search
    if (caloID != eventID)
      continue;

    // Cell rejected
    if (cutCell(caloCellNumber, caloAmplitude, caloType))
      continue;

    if (caloType == 0) {
      ret_event = interpretPHOSCell(caloCellNumber, caloAmplitude, std::move(ret_event));
    } else if (caloType == 1) {
      ret_event = interpretEMCALCell(caloCellNumber, caloAmplitude, std::move(ret_event));
    } else {
      Error("DataInterpreterAOD::interpretAODCaloCells", "Wrong calorimeter type");
    }
  }
  return ret_event;
}

std::unique_ptr<VisualisationEvent> DataInterpreterAOD::interpretEMCALCell(Int_t absID, Float_t amplitude, std::unique_ptr<VisualisationEvent> event)
{
  auto [iSupMod, iTower, iIphi, iIeta] = mEMCALGeom->GetCellIndex(absID);
  // It should not happen, but in case the OCDB file is not the correct one.
  if (iSupMod >= mEMCALGeom->GetNumberOfSuperModules())
    return event;

  // Gives label of cell in eta-phi position per each supermodule
  //auto [iphi, ieta] = mEMCALGeom->GetCellPhiEtaIndexInSModule(iSupMod, iTower, iIphi, iIeta);
  auto [eta, phi] = mEMCALGeom->EtaPhiFromIndex(absID);

  const auto& pos = mEMCALGeom->RelPosCellInSModule(absID);
  double posArr[3] = { pos.X(), pos.Y(), pos.Z() };

  VisualisationCaloCell cell(absID, iSupMod, posArr, amplitude, 1, phi, eta);

  event->addCaloCell(cell);
  return event;
}

std::unique_ptr<VisualisationEvent> DataInterpreterAOD::interpretPHOSCell(Int_t absID, Float_t amplitude, std::unique_ptr<VisualisationEvent> event)
{
  TVector3 xyz;
  Int_t relId[3];
  mPHOSGeom->AbsToRelNumbering(absID, relId);
  double x, z;
  mPHOSGeom->AbsIdToRelPosInModule(absID, x, z);
  double posArr[3] = { x, 0.0, z };

  // TODO: No Local2Global in PHOS geometry yet. Uncomment it once available.
  // mPHOSGeom->Local2Global(relId[0], x, z, xyz);
  // fHistoPH->Fill(xyz.Eta(),GetPhi(xyz.Phi()),amp);

  VisualisationCaloCell cell(absID, relId[0], posArr, amplitude, 0, relId[1], 0);

  //  printf("\t PHOS cell: ID %d, energy %2.2f,eta %2.2f, phi %2.2f, phi not modified: %2.2f\n",
  //         absID, amplitude, eta, cell.getPhi()*TMath::RadToDeg(), phi);
  //  printf("\t SM %d iEta %d,  iPhi %d, x %3.3f, y %3.3f, z %3.3f \n",
  //         iSupMod, ieta, iphi, pos.X(), pos.Y(), pos.Z());

  event->addCaloCell(cell);
  return event;
}

Bool_t DataInterpreterAOD::cutCell(Int_t absID, Float_t amplitude, Int_t type)
{
  //  if(fCaloCluster->GetNCells() < nMinCellsCut[type]) return kTRUE ;
  //  if(fCaloCluster->GetNCells() > nMaxCellsCut[type]) return kTRUE ;

  //  if(amplitude < mEnergyCut[type]) return kTRUE;

  //  if(fCaloCluster->GetM02()    < m02lowCut[type]) return kTRUE ;
  //  if(fCaloCluster->GetM02()    > m02higCut[type]) return kTRUE ;

  //  if(fCaloCluster->GetM20()    < m20lowCut[type]) return kTRUE ;

  if (amplitude < mEnergyCut[type] / 2.0 || absID < 0) {
    return kTRUE;
  }

  //  Int_t ism  =  -1, icol = -1, irow = -1;
  //  GetModuleNumberColAndRow(absID, type, ism,icol,irow)  ;
  //  Float_t eCross = GetECross(type, ism, icol, irow);
  //  if(1-eCross/eMax > mExoCut)
  //  {
  //    return kTRUE;
  //  }

  return kFALSE;
}

std::unique_ptr<VisualisationEvent> DataInterpreterAOD::interpretMuonTracks(TFile* AODFile, Int_t eventID)
{
  auto ret_event = std::make_unique<VisualisationEvent>(0, 0, 0, 0, "", 0);
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
    if (muonID > eventID + 10)
      break; // End event search
    if (muonID != eventID)
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

    VisualisationTrack visTrack(recTrack.fSign, 0.0, 0, 0, 0.0, 0.0, track_start, track_end, track_p, 0, 0.0, 0.0, 0.0, 0);

    for (Int_t i = 0; i < track->GetN(); ++i) {
      Float_t x, y, z;
      track->GetPoint(i, x, y, z);
      visTrack.addPolyPoint(x, y, z);
    }

    delete track;
    ret_event->addMuonTrack(visTrack);
  }
  delete ghostList;
  return ret_event;
}

} // namespace event_visualisation
} // namespace o2
