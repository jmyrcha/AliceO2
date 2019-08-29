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
/// \file    EventManager.cxx
/// \author  Jeremi Niedziela
/// \author julian.myrcha@cern.ch
/// \author p.nowakowski@cern.ch
/// \author Maja Kabus

#include "EventVisualisationView/EventManager.h"
#include "EventVisualisationBase/GeometryManager.h"
#include "EventVisualisationDataConverter/VisualisationEvent.h"
#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationView/MultiView.h"
#include "EventVisualisationBase/DataSource.h"
#include "EventVisualisationBase/DataInterpreter.h"
#include "EventVisualisationBase/DataSourceOffline.h"
#include "EventVisualisationDetectors/DataReaderVSD.h"
#include "EventVisualisationDetectors/DataReaderITS.h"
#include "EventVisualisationDetectors/CaloMatrix.h"
#include "EMCALBase/Geometry.h"
#include "PHOSBase/Geometry.h"

#include <FairLogger.h>

#include <TEveManager.h>
#include <TEveProjectionManager.h>
#include <TEveTrackPropagator.h>
#include <TSystem.h>
#include <TEnv.h>
#include <TEveElement.h>
#include <TGListTree.h>
#include <TEveQuadSet.h>
#include <TEveTrans.h>
#include <TGeoNode.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TEveTrans.h>
#include <TStyle.h>
#include <TEveCalo.h>
#include <TEveCaloData.h>
#include <TH2.h>

#include <iostream>

using namespace std;

namespace o2
{
namespace event_visualisation
{

EventManager* EventManager::mInstance = nullptr;

EventManager& EventManager::getInstance()
{
  if (mInstance == nullptr) {
    mInstance = new EventManager();
  }
  return *mInstance;
}

EventManager::EventManager() : TEveEventManager("Event", "")
{
  for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; i++) {
    mDataInterpreters[i] = nullptr;
    mDataReaders[i] = nullptr;
  }
}

void EventManager::Open()
{
  switch (mCurrentDataSourceType) {
    case SourceOnline:
      break;
    case SourceOffline: {
      DataSourceOffline* source = new DataSourceOffline();
      for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; i++) {
        if (mDataInterpreters[i] != nullptr) {
          mDataReaders[i]->open();
          source->registerReader(mDataReaders[i], static_cast<EVisualisationGroup>(i));
        }
      }
      setDataSource(source);
    } break;
    case SourceHLT:
      break;
  }
}

void EventManager::GotoEvent(Int_t no)
{
  //-1 means last event
  if (no == -1) {
    no = getDataSource()->GetEventCount() - 1;
  }

  this->mCurrentEvent = no;

  MultiView::getInstance()->destroyAllEvents();

  for (int i = 0; i < EVisualisationDataType::NdataTypes; ++i) {
    mDataTypeLists[i] = new TEveElementList(gDataTypeNames[i].c_str());
  }

  for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; ++i) {
    DataInterpreter* interpreter = mDataInterpreters[i];
    if (interpreter) {
      for (int dataType = 0; dataType < EVisualisationDataType::NdataTypes; ++dataType) {
        TObject* data = getDataSource()->getEventData(no, (EVisualisationGroup)i);
        std::unique_ptr<VisualisationEvent> event = interpreter->interpretDataForType(data, (EVisualisationDataType)dataType);
        displayVisualisationEvent(*event, gVisualisationGroupName[i]);
      }
    }
  }

  for (int i = 0; i < EVisualisationDataType::NdataTypes; ++i) {
    MultiView::getInstance()->registerElement(mDataTypeLists[i]);
  }

  MultiView::getInstance()->redraw3D();
}

void EventManager::NextEvent()
{
  Int_t event = (this->mCurrentEvent + 1) % getDataSource()->GetEventCount();
  GotoEvent(event);
}

void EventManager::PrevEvent()
{
  GotoEvent(this->mCurrentEvent - 1);
}

void EventManager::Close()
{
  delete this->mDataSource;
  this->mDataSource = nullptr;
}

void EventManager::AfterNewEventLoaded()
{
  TEveEventManager::AfterNewEventLoaded();
}

void EventManager::AddNewEventCommand(const TString& cmd)
{
  TEveEventManager::AddNewEventCommand(cmd);
}

void EventManager::RemoveNewEventCommand(const TString& cmd)
{
  TEveEventManager::RemoveNewEventCommand(cmd);
}

void EventManager::ClearNewEventCommands()
{
  TEveEventManager::ClearNewEventCommands();
}

EventManager::~EventManager()
{
  for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; i++) {
    if (mDataInterpreters[i] != nullptr) {
      delete mDataInterpreters[i];
      delete mDataReaders[i];

      mDataInterpreters[i] = nullptr;
      mDataReaders[i] = nullptr;
    }
  }
  mInstance = nullptr;
}

void EventManager::displayVisualisationEvent(VisualisationEvent& event, const std::string& detectorName)
{
  displayTracks(event, detectorName);
  displayClusters(event, detectorName);
  if (detectorName == "AOD") {
    displayCalo(event);
    displayMuonTracks(event);
  }
}

void EventManager::displayTracks(VisualisationEvent& event, const std::string& detectorName)
{
  size_t trackCount = event.getTrackCount();
  auto* list = new TEveTrackList(detectorName.c_str());
  list->IncDenyDestroy();

  for (size_t i = 0; i < trackCount; ++i) {
    VisualisationTrack track = event.getTrack(i);
    TEveRecTrackD t;
    double* p = track.getMomentum();
    t.fP = { p[0], p[1], p[2] };
    t.fSign = track.getCharge() > 0 ? 1 : -1;
    auto* vistrack = new TEveTrack(&t, &TEveTrackPropagator::fgDefault);
    vistrack->SetLineColor(kMagenta);
    size_t pointCount = track.getPointCount();
    vistrack->Reset(pointCount);

    for (size_t j = 0; j < pointCount; ++j) {
      auto point = track.getPoint(j);
      vistrack->SetNextPoint(point[0], point[1], point[2]);
    }
    list->AddElement(vistrack);
  }

  if (trackCount != 0) {
    mDataTypeLists[EVisualisationDataType::Tracks]->AddElement(list);
  }
}

void EventManager::displayMuonTracks(VisualisationEvent& event)
{
  size_t muonCount = event.getMuonTrackCount();
  std::cout << "Number of muon tracks: " << muonCount << std::endl;

  auto* muonList = new TEveElementList("AOD Muon tracks");
  muonList->SetTitle(Form("N=%d", muonCount));
  muonList->IncDenyDestroy();

  auto* match = new TEveTrackList("Matched");
  match->IncDenyDestroy();
  match->SetRnrPoints(kFALSE);
  match->SetRnrLine(kTRUE);
  match->SetLineColor(kGreen);
  setupMuonTrackPropagator(match->GetPropagator(), kTRUE, kTRUE);
  muonList->AddElement(match);

  auto* nomatch = new TEveTrackList("Not matched");
  nomatch->IncDenyDestroy();
  nomatch->SetRnrPoints(kFALSE);
  nomatch->SetRnrLine(kTRUE);
  nomatch->SetLineColor(kGreen);
  setupMuonTrackPropagator(nomatch->GetPropagator(), kTRUE, kFALSE);
  muonList->AddElement(nomatch);

  auto* ghost = new TEveTrackList("Ghost");
  ghost->IncDenyDestroy();
  ghost->SetRnrPoints(kFALSE);
  ghost->SetRnrLine(kTRUE);
  ghost->SetLineColor(kGreen);
  setupMuonTrackPropagator(ghost->GetPropagator(), kFALSE, kTRUE);
  muonList->AddElement(ghost);

  //int muonCount = event.getMuonTrackCount();
  for (size_t i = 0; i < muonCount; ++i) {
    VisualisationTrack track = event.getMuonTrack(i);
    TEveRecTrackD t;
    double* p = track.getMomentum();
    t.fP = { p[0], p[1], p[2] };
    t.fSign = track.getCharge() > 0 ? 1 : -1;
    auto* vistrack = new TEveTrack(&t, ghost->GetPropagator()); // &TEveTrackPropagator::fgDefault);
    size_t pointCount = track.getPointCount();
    vistrack->Reset(pointCount);

    for (size_t j = 0; j < pointCount; ++j) {
      auto point = track.getPoint(j);
      vistrack->SetNextPoint(point[0], point[1], point[2]);
    }
    vistrack->SetAttLineAttMarker(ghost);
    ghost->AddElement(vistrack);
  }

  if (muonCount != 0 || muonCount != 0) {
    mDataTypeLists[EVisualisationDataType::Muon]->AddElement(muonList);
  }
}

void EventManager::setupMuonTrackPropagator(TEveTrackPropagator* prop, Bool_t tracker, Bool_t trigger)
{
  // TODO: Set magnetic field properly
  //  if (AliMUONTrackExtrap::IsFieldON())
  //  {
  //    prop->SetMagFieldObj(new AliEveMagField);
  //  }
  //  else
  //  {
  //    prop->SetMagField(0.0);
  //  }
  prop->SetMagField(0.5);
  prop->SetStepper(TEveTrackPropagator::kRungeKutta);

  prop->SetMaxR(1000);
  // TODO: Find corresponding constants in O2
  //if (trigger) prop->SetMaxZ(-AliMUONConstants::DefaultChamberZ(13) + 10.);
  //else prop->SetMaxZ(-AliMUONConstants::MuonFilterZBeg());

  // Go through pathmarks
  prop->SetFitDaughters(kFALSE);
  prop->SetFitReferences(kTRUE);
  prop->SetFitDecay(kFALSE);
  prop->SetFitCluster2Ds(kFALSE);

  // Render the ref pathmarks
  prop->SetRnrReferences(kTRUE);
  prop->RefPMAtt().SetMarkerSize(0.5);
  if (trigger)
    prop->RefPMAtt().SetMarkerColor(kGreen);
  else
    prop->RefPMAtt().SetMarkerColor(kAzure);

  // Render first vertex
  if (tracker) {
    prop->SetRnrFV(kTRUE);
    if (trigger)
      prop->RefFVAtt().SetMarkerColor(kGreen);
    else
      prop->RefFVAtt().SetMarkerColor(kAzure);
  }
}

void EventManager::displayClusters(VisualisationEvent& event, const std::string& detectorName)
{
  size_t clusterCount = event.getClusterCount();
  auto* point_list = new TEvePointSet(detectorName.c_str());
  point_list->IncDenyDestroy();
  point_list->SetMarkerColor(kBlue);

  for (size_t i = 0; i < clusterCount; ++i) {
    VisualisationCluster cluster = event.getCluster(i);
    point_list->SetNextPoint(cluster.X(), cluster.Y(), cluster.Z());
  }

  if (clusterCount != 0) {
    mDataTypeLists[EVisualisationDataType::Clusters]->AddElement(point_list);
  }
}

void EventManager::displayCalo(VisualisationEvent& event)
{
  size_t caloCount = event.getCaloCellsCount();
  if (caloCount == 0)
    return;

  auto* caloList = new TEveElementList("3D Histogram");
  caloList->IncDenyDestroy();
  auto* emcalList = new TEveElementList("EMCAL");
  emcalList->IncDenyDestroy();
  auto* phosList = new TEveElementList("PHOS");
  phosList->IncDenyDestroy();

  // 3D calorimeter histogram
  double pi = TMath::Pi();
  TEveCaloDataHist* data = new TEveCaloDataHist();
  TH2F* histoEM = new TH2F("histoEMcell", "EMCal Cell #eta vs #phi vs E",
                           100, -1.5, 1.5, 80, -pi, pi);
  TH2F* histoPH = new TH2F("histoPHcell", "PHOS Cell #eta vs #phi vs E",
                           100, -1.5, 1.5, 80, -pi, pi);
  data->AddHistogram(histoEM);
  data->RefSliceInfo(0).Setup("EMCell:", 0, kOrange + 7);
  data->AddHistogram(histoPH);
  data->RefSliceInfo(1).Setup("PHCell:", 0, kYellow);

  data->GetEtaBins()->SetTitleFont(120);
  data->GetEtaBins()->SetTitle("h");
  data->GetPhiBins()->SetTitleFont(120);
  data->GetPhiBins()->SetTitle("f");
  data->IncDenyDestroy();

  TEveCalo3D* calo3d = new TEveCalo3D(data);
  calo3d->SetBarrelRadius(600);
  calo3d->SetEndCapPos(550);
  calo3d->SetMaxTowerH(300);
  calo3d->SetFrameTransparency(100);
  gEve->AddElement(calo3d, caloList);

  // Version with TEveCalo3D as separate top node on the list (outside 'Event' and 'EMCAL').
  // Warning: the elements are not properly removed when moving to next / prev event

  //TEveScene* g_histo2d_s2 = gEve->SpawnNewScene("3D Histogram", "3D Histogram");
  //gEve->GetDefaultViewer()->AddScene(g_histo2d_s2);
  //MultiView::getInstance()->getView(MultiView::EViews::View3d)->AddScene(g_histo2d_s2);
  //g_histo2d_s2->SetElementName("3D Histogram Scene");
  //g_histo2d_s2->AddElement(calo3d);
  //MultiView::getInstance()->registerElement(calo3d);

  // Warning: Geometries need to be initialised before
  // We assume that they were initialised in AOD interpreter
  const auto& emcalGeom = o2::emcal::Geometry::GetInstance();
  const auto& phosGeom = o2::phos::Geometry::GetInstance();

  int numberOfSuperModules = emcalGeom->GetNumberOfSuperModules();
  TEveQuadSet* emcalQuads[numberOfSuperModules];
  memset(emcalQuads, 0, numberOfSuperModules * sizeof(TEveQuadSet*));

  // TODO: Any way not to hardcode number of PHOS modules?
  TEveQuadSet* phosQuads[4];
  memset(phosQuads, 0, 4 * sizeof(TEveQuadSet*));

  // Quad size
  Float_t quadSizeEMCAL = 6; // cm, tower side size
  Float_t quadSizePHOS = 2.2;

  for (Int_t sm = 0; sm < numberOfSuperModules; ++sm) {
    emcalQuads[sm] = new TEveQuadSet(Form("SM %d", sm + 1));

    // Warning: It will crash if there is no matrix.
    // We assume all matrices are already set by AOD interpreter.
    setCaloQuadSet(quadSizeEMCAL, emcalGeom->GetMatrixForSuperModule(sm), emcalQuads[sm]);
    gEve->AddElement(emcalQuads[sm], emcalList);
  }

  for (Int_t mod = 0; mod < 4; ++mod) {
    phosQuads[mod] = new TEveQuadSet(Form("Mod %d", mod + 1)); // Why just mod in AliRoot?

    // TODO: Setting PHOS matrices once it will be possible

    setCaloQuadSet(quadSizeEMCAL, new TGeoHMatrix, phosQuads[mod]);
    gEve->AddElement(phosQuads[mod], phosList);
  }

  for (size_t i = 0; i < caloCount; ++i) {
    VisualisationCaloCell caloCell = event.getCaloCell(i);

    // Cells = blue quads
    int module = caloCell.getModule();
    if (emcalQuads[module]) {
      emcalQuads[module]->AddQuad(caloCell.Y(), caloCell.Z());
      emcalQuads[module]->QuadValue(caloCell.getAmplitude() * 1000);
    }

    // Histogram = orange boxes
    float eta = caloCell.getEta();
    if (TMath::Abs(eta) < 0.7) {
      // FIXME: Scale - looks almost flat, but in AliRoot there is the same amplitude value in the histogram??
      histoEM->Fill(eta, caloCell.getPhi(), caloCell.getAmplitude());
      //      std::cout << "Histo for cell: " << caloCell.getAbsID() << " eta: " << eta << " phi: " << caloCell.getPhi() << " energy: " << caloCell.getAmplitude() << std::endl;
    }
  }

  mDataTypeLists[EVisualisationDataType::Calo]->AddElement(caloList);
  mDataTypeLists[EVisualisationDataType::Calo]->AddElement(emcalList);
  mDataTypeLists[EVisualisationDataType::Calo]->AddElement(phosList);
}

void EventManager::setCaloQuadSet(const Float_t quadSize, const TGeoHMatrix* matrix, TEveQuadSet* quadSet)
{
  quadSet->SetOwnIds(kTRUE);
  // Type of object to be displayed, rectangle with cell size
  quadSet->Reset(TEveQuadSet::kQT_RectangleYZFixedDimX, kFALSE, 32);
  quadSet->SetDefWidth(quadSize);
  quadSet->SetDefHeight(quadSize);
  quadSet->RefMainTrans().SetFrom(*matrix);

  // Define energy range for the color palette
  Int_t maxEMCalE = 2000; // MeV
  Int_t minEMCalE = 100;  // MeV

  gStyle->SetPalette(1, 0);
  TEveRGBAPalette* pal = new TEveRGBAPalette(minEMCalE, maxEMCalE);
  quadSet->SetPalette(pal);

  quadSet->RefitPlex();
}

void EventManager::registerDetector(DataReader* reader, DataInterpreter* interpreter, EVisualisationGroup type)
{
  mDataReaders[type] = reader;
  mDataInterpreters[type] = interpreter;
}

} // namespace event_visualisation
} // namespace o2
