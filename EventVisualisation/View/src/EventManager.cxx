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
/// \author  julian.myrcha@cern.ch
/// \author  p.nowakowski@cern.ch
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#include "EventVisualisationView/EventManager.h"
#include "EventVisualisationDataConverter/VisualisationEvent.h"
#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationView/MultiView.h"
#include "EventVisualisationBase/DataSource.h"
#include "EventVisualisationBase/DataInterpreter.h"
#include "EventVisualisationBase/DataSourceOffline.h"
#include "ReconstructionDataFormats/PID.h"
#include "EMCALBase/Geometry.h"
#include "PHOSBase/Geometry.h"

#include "FairLogger.h"

#include <TEveManager.h>
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
#include <TStyle.h>
#include <TEveCalo.h>
#include <TEveCaloData.h>

#include <string>
#include <cctype>
#include <vector>

namespace o2
{
namespace event_visualisation
{

EventManager* EventManager::sInstance = nullptr;

EventManager& EventManager::getInstance()
{
  if (sInstance == nullptr) {
    sInstance = new EventManager();
  }
  return *sInstance;
}

EventManager::EventManager() : TEveEventManager("Event", "")
{
  for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; i++) {
    mDataInterpreters[i] = nullptr;
    mDataReaders[i] = nullptr;
  }

  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);
  mWidth = settings.GetValue("tracks.width", 1);

  // 3D calorimeter histogram
  double pi = TMath::Pi();
  mEmcalHistogram = new TH2F("histoEMcell", "EMCal Cell #eta vs #phi vs E",
                             100, -1.5, 1.5, 80, -pi, pi);
  mPhosHistogram = new TH2F("histoPHcell", "PHOS Cell #eta vs #phi vs E",
                            100, -1.5, 1.5, 80, -pi, pi);
}

void EventManager::Open()
{
  switch (mCurrentDataSourceType) {
    case SourceOnline:
      break;
    case SourceOffline: {
      DataSourceOffline* source = new DataSourceOffline();
      for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; i++) {
        if (mDataReaders[i] != nullptr) {
          source->registerReader(mDataReaders[i], static_cast<EVisualisationGroup>(i));
          mDataReaders[i]->open();
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
    no = getDataSource()->getEventCount() - 1;
  }

  mCurrentEventNumber = no;

  MultiView::getInstance()->destroyAllEvents();

  for (int i = 0; i < EVisualisationDataType::NdataTypes; ++i) {
    mDataTypeLists[i] = new TEveElementList(gDataTypeNames[i].c_str());
  }

  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);

  for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; ++i) {
    DataInterpreter* interpreter = mDataInterpreters[i];
    if (interpreter) {
      auto dataSource = getDataSource();
      if (dataSource->hasEventData(no, (EVisualisationGroup)i)) {
        mCurrentEvent = std::make_unique<VisualisationEvent>(0, 0, 0, 0, "", 0);
        for (int dataType = 0; dataType < EVisualisationDataType::NdataTypes; ++dataType) {
          if (settings.GetValue(Form("%c%s.show", std::tolower(gDataTypeNames[dataType][0]), gDataTypeNames[dataType].substr(1).data()), false)) {
            TObject* data = dataSource->getEventData(no, (EVisualisationGroup)i, (EVisualisationDataType)dataType);
            if (data) {
              interpreter->interpretDataForType(data, (EVisualisationDataType)dataType, *mCurrentEvent);
            }
          }
        }
        displayVisualisationEvent(*mCurrentEvent, gVisualisationGroupName[i]);
      }
    }
  }

  for (int i = 0; i < EVisualisationDataType::NdataTypes; ++i) {
    MultiView::getInstance()->registerElement(mDataTypeLists[i]);
  }

  MultiView::getInstance()->redraw3D();

  if (settings.GetValue("tracks.animate", false))
    animateTracks();
}

void EventManager::NextEvent()
{
  Int_t event = (mCurrentEventNumber + 1) % getDataSource()->getEventCount();
  GotoEvent(event);
}

void EventManager::PrevEvent()
{
  GotoEvent(mCurrentEventNumber - 1);
}

void EventManager::Close()
{
  delete mDataSource;
  mDataSource = nullptr;
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
  MultiView::getInstance()->destroyAllEvents();
  for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; i++) {
    delete mDataInterpreters[i];
    mDataInterpreters[i] = nullptr;
    delete mDataReaders[i];
    mDataReaders[i] = nullptr;
  }
  delete mDataSource;
  mDataSource = nullptr;
  delete mEmcalHistogram;
  mEmcalHistogram = nullptr;
  delete mPhosHistogram;
  mPhosHistogram = nullptr;
  sInstance = nullptr;
}

void EventManager::displayVisualisationEvent(VisualisationEvent& event, const std::string& detectorName)
{
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);
  if (settings.GetValue("tracks.show", false)) {
    displayTracks(event, detectorName);
  }
  if (settings.GetValue("clusters.show", false)) {
    displayClusters(event, detectorName);
  }

  if (detectorName == "AOD") {
    if (settings.GetValue("calo.show", false)) {
      displayCalo(event);
    }
    if (settings.GetValue("muon.show", false)) {
      displayMuonTracks(event);
    }
  }
}

void EventManager::displayTracks(VisualisationEvent& event, const std::string& detectorName)
{
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);

  if (settings.GetValue("tracks.byPt.show", false)) {
    displayTracksByPt(event, detectorName);
    return;
  }
  if (settings.GetValue("tracks.byType.show", false) && detectorName == "AOD") {
    displayTracksByType(event, detectorName);
    return;
  }

  size_t trackCount = event.getTrackCount();
  if (trackCount == 0)
    return;

  auto* list = new TEveTrackList(detectorName.c_str());
  list->IncDenyDestroy();
  list->SetLineWidth(mWidth);

  const float magF = 0.1 * 5; // FIXME: Get it from OCDB / event
  const float maxR = settings.GetValue("tracks.animate", false) ? 0 : 520;
  auto prop = list->GetPropagator();
  prop->SetMagField(magF);
  prop->SetMaxR(maxR);

  for (size_t i = 0; i < trackCount; ++i) {
    VisualisationTrack track = event.getTrack(i);
    TEveRecTrackD t;
    double* p = track.getMomentum();
    t.fP = {p[0], p[1], p[2]};
    t.fSign = track.getCharge() > 0 ? 1 : -1;
    auto* vistrack = new TEveTrack(&t, prop);
    vistrack->SetMainColor(kMagenta);
    size_t pointCount = track.getPointCount();
    vistrack->Reset(pointCount);

    for (size_t j = 0; j < pointCount; ++j) {
      auto point = track.getPoint(j);
      vistrack->SetNextPoint(point[0], point[1], point[2]);
    }
    list->AddElement(vistrack);
  }

  mDataTypeLists[EVisualisationDataType::Tracks]->AddElement(list);
}

void EventManager::displayTracksByPt(VisualisationEvent& event, const std::string& detectorName)
{
  size_t trackCount = event.getTrackCount();
  if (trackCount == 0)
    return;

  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);

  auto* trackList = new TEveElementList(Form("%s tracks by Pt", detectorName.c_str()));
  trackList->IncDenyDestroy();

  const Int_t nCont = 6;
  const float magF = 0.1 * 5; // FIXME: Get it from OCDB / event
  const float maxR = settings.GetValue("tracks.animate", false) ? 0 : 520;

  TEveTrackList* tl[nCont];
  Int_t tc[nCont];
  Int_t count = 0;

  Color_t colors[nCont];
  // default color scheme by type:
  colors[0] = kGreen;
  colors[1] = kSpring + 10;
  colors[2] = kYellow + 1;
  colors[3] = kOrange;
  colors[4] = kOrange - 3;
  colors[5] = kRed;

  tl[0] = new TEveTrackList("pt < 0.2");
  tl[1] = new TEveTrackList("0.2 < pt < 0.4");
  tl[2] = new TEveTrackList("0.4 < pt < 0.7");
  tl[3] = new TEveTrackList("0.7 < pt < 1.1");
  tl[4] = new TEveTrackList("1.1 < pt < 1.6");
  tl[5] = new TEveTrackList("pt > 1.6");

  for (int i = 0; i < nCont; i++) {
    tc[i] = 0;
    auto prop = tl[i]->GetPropagator();
    prop->SetMagField(magF);
    prop->SetMaxR(maxR);
    tl[i]->SetMainColor(colors[i]);
    tl[i]->SetLineWidth(mWidth);
    trackList->AddElement(tl[i]);
  }

  for (size_t i = 0; i < trackCount; ++i) {
    VisualisationTrack track = event.getTrack(i);

    int index;
    double pt = abs(track.getSignedPt());

    if (pt <= 0.2)
      index = 0;
    else if (pt > 0.2 && pt <= 0.4)
      index = 1;
    else if (pt > 0.4 && pt <= 0.7)
      index = 2;
    else if (pt > 0.7 && pt <= 1.1)
      index = 3;
    else if (pt > 1.1 && pt <= 1.6)
      index = 4;
    else
      index = 5;

    TEveTrackList* tlist = tl[index];
    ++tc[index];
    ++count;

    TEveRecTrackD t;
    double* p = track.getMomentum();
    t.fP = {p[0], p[1], p[2]};
    t.fSign = track.getCharge() > 0 ? 1 : -1;
    auto* vistrack = new TEveTrack(&t, tlist->GetPropagator());
    size_t pointCount = track.getPointCount();
    vistrack->Reset(pointCount);

    for (size_t j = 0; j < pointCount; ++j) {
      auto point = track.getPoint(j);
      vistrack->SetNextPoint(point[0], point[1], point[2]);
    }
    vistrack->SetName(Form("Track no = %lu, pt = %f", i, pt));
    vistrack->SetTitle(Form("Track no = %lu, pt = %f", i, pt));
    vistrack->SetAttLineAttMarker(tlist);
    tlist->AddElement(vistrack);
  }

  for (Int_t ti = 0; ti < nCont; ++ti) {
    TEveTrackList* tlist = tl[ti];
    tlist->SetName(Form("%s [%d]", tlist->GetName(), tlist->NumChildren()));
    tlist->SetTitle(Form("N tracks = %d", tc[ti]));
    tlist->MakeTracks();
  }
  trackList->SetTitle(Form("N all tracks = %d", count));

  mDataTypeLists[EVisualisationDataType::Tracks]->AddElement(trackList);
}

void EventManager::displayTracksByType(VisualisationEvent& event, const std::string& detectorName)
{
  size_t trackCount = event.getTrackCount();
  if (trackCount == 0)
    return;

  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);

  auto* trackList = new TEveElementList(Form("%s tracks by type", detectorName.c_str()));
  trackList->IncDenyDestroy();

  const Int_t nCont = 15;
  const float magF = 0.1 * 5; // FIXME: Get it from OCDB / event
  const float maxR = settings.GetValue("tracks.animate", false) ? 0 : 520;

  TEveTrackList* tl[nCont];
  Int_t tc[nCont];
  Int_t count = 0;

  Color_t colors[15];
  // default color scheme by type:
  colors[0] = settings.GetValue("tracks.byType.electron", 600);
  colors[1] = settings.GetValue("tracks.byType.muon", 416);
  colors[2] = settings.GetValue("tracks.byType.pion", 632);
  colors[3] = settings.GetValue("tracks.byType.kaon", 400);
  colors[4] = settings.GetValue("tracks.byType.proton", 797);
  colors[5] = settings.GetValue("tracks.byType.deuteron", 797);
  colors[6] = settings.GetValue("tracks.byType.triton", 797);
  colors[7] = settings.GetValue("tracks.byType.he3", 797);
  colors[8] = settings.GetValue("tracks.byType.alpha", 403);
  colors[9] = settings.GetValue("tracks.byType.photon", 0);
  colors[10] = settings.GetValue("tracks.byType.pi0", 616);
  colors[11] = settings.GetValue("tracks.byType.neutron", 900);
  colors[12] = settings.GetValue("tracks.byType.kaon0", 801);
  colors[13] = settings.GetValue("tracks.byType.elecon", 920);
  colors[14] = settings.GetValue("tracks.byType.unknown", 920);

  tl[0] = new TEveTrackList("Electrons");
  tl[1] = new TEveTrackList("Muons");
  tl[2] = new TEveTrackList("Pions");
  tl[3] = new TEveTrackList("Kaons");
  tl[4] = new TEveTrackList("Protons");
  tl[5] = new TEveTrackList("Deuterons");
  tl[6] = new TEveTrackList("Tritons");
  tl[7] = new TEveTrackList("He3");
  tl[8] = new TEveTrackList("Alpha");
  tl[9] = new TEveTrackList("Photons");
  tl[10] = new TEveTrackList("Pi0");
  tl[11] = new TEveTrackList("Neutrons");
  tl[12] = new TEveTrackList("Kaon0");
  tl[13] = new TEveTrackList("EleCon");
  tl[14] = new TEveTrackList("Unknown");

  for (int i = 0; i < nCont; i++) {
    tc[i] = 0;
    auto prop = tl[i]->GetPropagator();
    prop->SetMagField(magF);
    prop->SetMaxR(maxR);
    tl[i]->SetMainColor(colors[i]);
    tl[i]->SetLineWidth(mWidth);
    trackList->AddElement(tl[i]);
  }

  bool shading = true;
  int shade = -3;
  VisualisationTrack firstTrack = event.getTrack(0);
  int firstPid = firstTrack.getPID().getID();

  // Check if all tracks have the same PID. If no, turn off shading
  for (size_t i = 1; i < trackCount; ++i) {
    const VisualisationTrack& track = event.getTrack(i);
    if (track.getPID().getID() != firstPid) {
      shading = false;
      break;
    }
  }

  for (size_t i = 0; i < trackCount; ++i) {
    VisualisationTrack track = event.getTrack(i);

    if (!trackSelected(track))
      continue;

    int pid = track.getPID().getID();

    TEveTrackList* tlist = tl[pid];
    ++tc[pid];
    ++count;

    TEveRecTrackD t;
    double* p = track.getMomentum();
    t.fP = {p[0], p[1], p[2]};
    t.fSign = track.getCharge() > 0 ? 1 : -1;
    auto* vistrack = new TEveTrack(&t, tlist->GetPropagator());
    size_t pointCount = track.getPointCount();
    vistrack->Reset(pointCount);

    for (size_t j = 0; j < pointCount; ++j) {
      auto point = track.getPoint(j);
      vistrack->SetNextPoint(point[0], point[1], point[2]);
    }

    if (shading) {
      if ((kGreen + shade) < 0) {
        shade = 0;
      }
      vistrack->SetMainColor(kGreen + shade);
      shade++;
      if (shade > 3)
        shade = -3;
    } else {
      vistrack->SetAttLineAttMarker(tlist);
    }

    vistrack->SetName(Form("Track no = %lu, PID = %d", i, pid));
    vistrack->SetTitle(getTrackTitle(track));
    tlist->AddElement(vistrack);
  }

  for (Int_t ti = 0; ti < nCont; ++ti) {
    TEveTrackList* tlist = tl[ti];
    tlist->SetName(Form("%s [%d]", tlist->GetName(), tlist->NumChildren()));
    tlist->SetTitle(Form("N tracks = %d", tc[ti]));
    tlist->MakeTracks();
  }
  trackList->SetTitle(Form("N all tracks = %d", count));

  mDataTypeLists[EVisualisationDataType::Tracks]->AddElement(trackList);
}

// From AliRoot/EVE/EveBase/AliEveESDTracks.cxx
// Used for drawing ESD tracks by PID, TPC, ITS, ITS standalone tracks
TString EventManager::getTrackTitle(VisualisationTrack& track)
{
  TString s;

  Int_t index = track.getID();
  TString idx(index == kMinInt ? "<undef>" : Form("%d", index));

  double* v = track.getVertex();
  double* p = track.getMomentum();
  double pt = TMath::Abs(track.getSignedPt());
  double ptsig = TMath::Sqrt(track.getC1Pt21Pt2());
  double ptsq = pt * pt;
  double ptm = pt / (1.0 + pt * ptsig);
  double ptM = pt / (1.0 - pt * ptsig);

  // FIXME: Label and PDG to be completed once available
  s = Form(
    "Index=%s, Label=%s\nChg=%d, Pdg=%d\n"
    "pT = %.3f + %.3f - %.3f [%.3f]\n"
    "P  = (%.3f, %.3f, %.3f)\n"
    "V  = (%.3f, %.3f, %.3f)\n",
    idx.Data(), "0", track.getCharge(), 0,
    pt, ptM - pt, pt - ptm, ptsig * ptsq,
    p[0], p[1], p[2],
    v[0], v[1], v[2]);

  int o;
  s += "Det (in, out, refit, pid):\n";
  o = ERecoFlag::kITSin;
  s += Form("ITS (%d, %d, %d, %d)  ", track.isRecoFlagSet(o), track.isRecoFlagSet(o << 1), track.isRecoFlagSet(o << 2), track.isRecoFlagSet(o << 3));
  o = ERecoFlag::kTPCin;
  s += Form("TPC(%d, %d, %d, %d)\n", track.isRecoFlagSet(o), track.isRecoFlagSet(o << 1), track.isRecoFlagSet(o << 2), track.isRecoFlagSet(o << 3));
  o = ERecoFlag::kTRDin;
  s += Form("TRD(%d, %d, %d, %d) ", track.isRecoFlagSet(o), track.isRecoFlagSet(o << 1), track.isRecoFlagSet(o << 2), track.isRecoFlagSet(o << 3));
  o = ERecoFlag::kTOFin;
  s += Form("TOF(%d, %d, %d, %d)\n", track.isRecoFlagSet(o), track.isRecoFlagSet(o << 1), track.isRecoFlagSet(o << 2), track.isRecoFlagSet(o << 3));
  o = ERecoFlag::kHMPIDout;
  s += Form("HMPID(out = %d, pid = %d)\n", track.isRecoFlagSet(o), track.isRecoFlagSet(o << 1));
  s += Form("ESD pid = %d", track.isRecoFlagSet(ERecoFlag::kESDpid));

  return s;
}

bool EventManager::trackSelected(const VisualisationTrack& track)
{
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);
  std::string trackSelection = settings.GetValue("tracks.selection", "");

  if (trackSelection == "ITSin_noTPCin") {
    return track.isRecoFlagSet(ERecoFlag::kITSin) && !track.isRecoFlagSet(ERecoFlag::kTPCin);
  }
  if (trackSelection == "noTISpureSA") {
    return !track.isRecoFlagSet(ERecoFlag::kITSpureSA);
  }
  if (trackSelection == "TPCrefit") {
    return track.isRecoFlagSet(ERecoFlag::kTPCrefit);
  }
  if (trackSelection == "TPCrefit_ITSrefit") {
    return track.isRecoFlagSet(ERecoFlag::kTPCrefit) && track.isRecoFlagSet(ERecoFlag::kITSrefit);
  }
  return true;
}

void EventManager::animateTracks()
{
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);

  std::string tracksString = "";
  if (settings.GetValue("tracks.byPt.show", false)) {
    tracksString = " tracks by Pt";
  } else if (settings.GetValue("tracks.byType.show", false)) {
    tracksString = " tracks by type";
  }

  for (int i = 0; i < EVisualisationGroup::NvisualisationGroups; ++i) {
    if (mDataInterpreters[i]) {
      TEveElementList* tracks = (TEveElementList*)mDataTypeLists[EVisualisationDataType::Tracks]->FindChild(Form("%s%s", gVisualisationGroupName[i].c_str(), tracksString.c_str()));
      if (tracks) {
        std::vector<TEveTrackPropagator*> propagators;

        for (TEveElement::List_i i = tracks->BeginChildren(); i != tracks->EndChildren(); i++) {
          TEveTrackList* trkList = ((TEveTrackList*)*i);
          propagators.push_back(trkList->GetPropagator());
        }
        int animationSpeed = settings.GetValue("tracks.animation.speed", 10);
        for (int R = 0; R <= 520; R += animationSpeed) {
          animationSpeed -= 10;
          if (animationSpeed < 10)
            animationSpeed = 10;
          for (int propIter = 0; propIter < propagators.size(); propIter++) {
            propagators[propIter]->SetMaxR(R);
          }
          // Unlike in alieve, full redraw is needed to update the visualisation (slows down the animation)
          gSystem->ProcessEvents();
          gEve->FullRedraw3D();
        }
      }
    }
  }
}

void EventManager::displayMuonTracks(VisualisationEvent& event)
{
  size_t muonCount = event.getMuonTrackCount();

  auto* match = new TEveTrackList("Matched");
  match->IncDenyDestroy();
  match->SetRnrPoints(kFALSE);
  match->SetRnrLine(kTRUE);
  match->SetMainColor(kGreen);
  setupMuonTrackPropagator(match->GetPropagator(), kTRUE, kTRUE);

  auto* nomatch = new TEveTrackList("Not matched");
  nomatch->IncDenyDestroy();
  nomatch->SetRnrPoints(kFALSE);
  nomatch->SetRnrLine(kTRUE);
  nomatch->SetMainColor(kGreen);
  setupMuonTrackPropagator(nomatch->GetPropagator(), kTRUE, kFALSE);

  auto* ghost = new TEveTrackList("Ghost");
  ghost->IncDenyDestroy();
  ghost->SetRnrPoints(kFALSE);
  ghost->SetRnrLine(kTRUE);
  ghost->SetMainColor(kGreen);
  setupMuonTrackPropagator(ghost->GetPropagator(), kFALSE, kTRUE);

  for (size_t i = 0; i < muonCount; ++i) {
    VisualisationTrack track = event.getMuonTrack(i);
    TEveRecTrackD t;
    double* p = track.getMomentum();
    t.fP = {p[0], p[1], p[2]};
    t.fV = track.getVertex();
    t.fSign = track.getCharge() > 0 ? 1 : -1;
    t.fIndex = track.getID();
    auto* vistrack = new TEveTrack(&t, ghost->GetPropagator());
    size_t pointCount = track.getPointCount();
    vistrack->Reset(pointCount);

    for (size_t j = 0; j < pointCount; ++j) {
      auto point = track.getPoint(j);
      vistrack->SetNextPoint(point[0], point[1], point[2]);
    }
    vistrack->SetAttLineAttMarker(ghost);

    // For ghost muon tracks
    vistrack->SetName("mu");
    vistrack->SetTitle("Trigger only");
    // For other muon tracks
    //vistrack->SetName(Form("%cmu", track.getCharge() > 0 ? '+' : '-'));
    //vistrack->SetStdTitle();

    ghost->AddElement(vistrack);
  }

  ghost->SetTitle(Form("N all tracks = %lu", muonCount));

  if (muonCount != 0) {
    mDataTypeLists[EVisualisationDataType::Muon]->AddElement(match);
    mDataTypeLists[EVisualisationDataType::Muon]->AddElement(nomatch);
    mDataTypeLists[EVisualisationDataType::Muon]->AddElement(ghost);
  }
}

void EventManager::setupMuonTrackPropagator(TEveTrackPropagator* prop, Bool_t tracker, Bool_t trigger)
{
  // TODO: Set magnetic field properly
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

  // Warning: Geometries need to be initialised before
  // We assume that they were initialised in AOD interpreter
  const auto& emcalGeom = o2::emcal::Geometry::GetInstance();
  const auto& phosGeom = o2::phos::Geometry::GetInstance();

  mEmcalHistogram->Reset();
  mPhosHistogram->Reset();

  int numberOfSuperModules = emcalGeom->GetNumberOfSuperModules();
  TEveQuadSet* emcalQuads[numberOfSuperModules];
  memset(emcalQuads, 0, numberOfSuperModules * sizeof(TEveQuadSet*));

  // TODO: Any way not to hardcode number of PHOS modules?
  TEveQuadSet* phosQuads[4];
  memset(phosQuads, 0, 4 * sizeof(TEveQuadSet*));

  // Quad size
  float quadSizeEMCAL = 6; // cm, tower side size
  float quadSizePHOS = 2.2;

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

    setCaloQuadSet(quadSizePHOS, new TGeoHMatrix, phosQuads[mod]);
    gEve->AddElement(phosQuads[mod], phosList);
  }

  for (size_t i = 0; i < caloCount; ++i) {
    VisualisationCaloCell caloCell = event.getCaloCell(i);

    int module = caloCell.getModule();
    float eta = caloCell.getEta();
    if (caloCell.getType() == 1) {
      // Cells = blue quads
      if (emcalQuads[module]) {
        emcalQuads[module]->AddQuad(caloCell.Y(), caloCell.Z());
        emcalQuads[module]->QuadValue(caloCell.getAmplitude() * 1000);
      }

      // Histogram = orange boxes
      if (TMath::Abs(eta) < 0.7) {
        mEmcalHistogram->Fill(eta, caloCell.getPhi(), caloCell.getAmplitude());
        //      printf("\t CaloCell %d, energy %2.2f,eta %2.2f, phi %2.2f\n",
        //             caloCell.getAbsID(), caloCell.getAmplitude(), caloCell.getEta(), caloCell.getPhi()*TMath::RadToDeg());
      } else {
        LOG(WARNING) << "Wrong eta value for EMCAL calorimeter cell, active workaround!!!";
      }
    }
    // TODO: PHOS Geometry not set yet, so not possible to use
    //    else {
    //      if(phosQuads[module]) {
    //        phosQuads[module]->AddQuad(caloCell.X(), caloCell.Z());
    //        phosQuads[module]->QuadValue(caloCell.getAmplitude() * 1000);
    //      }
    //    }
  }

  TEveCaloDataHist* data = new TEveCaloDataHist();
  data->AddHistogram(mEmcalHistogram);
  data->RefSliceInfo(0).Setup("EMCell:", 0, kOrange + 7);
  data->AddHistogram(mPhosHistogram);
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

  mDataTypeLists[EVisualisationDataType::Calo]->AddElement(caloList);
  mDataTypeLists[EVisualisationDataType::Calo]->AddElement(emcalList);
  mDataTypeLists[EVisualisationDataType::Calo]->AddElement(phosList);
}

void EventManager::setCaloQuadSet(const float quadSize, const TGeoHMatrix* matrix, TEveQuadSet* quadSet)
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
