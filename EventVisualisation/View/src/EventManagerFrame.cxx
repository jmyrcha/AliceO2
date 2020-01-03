// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file EventManagerFrame.cxx
/// \brief GUI (bottom buttons) for visualisation
/// \author julian.myrcha@cern.ch
/// \author p.nowakowski@cern.ch

#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationBase/GeometryManager.h"
#include "EventVisualisationView/EventManagerFrame.h"
#include "EventVisualisationView/MultiView.h"
#include "EventVisualisationBase/DataSourceOffline.h"

#include <TEveManager.h>
#include <TSystem.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>

#include <Rtypes.h>

ClassImp(o2::event_visualisation::EventManagerFrame);

namespace o2
{
namespace event_visualisation
{

EventManagerFrame::EventManagerFrame(o2::event_visualisation::EventManager& eventManager)
  : TGMainFrame(gClient->GetRoot(), 400, 100, kVerticalFrame)
{
  mEventManager = &eventManager;

  const TString cls("o2::event_visualisation::EventManagerFrame");
  TGTextButton* b = nullptr;
  TGHorizontalFrame* f = new TGHorizontalFrame(this);
  {
    Int_t width = 50;
    this->AddFrame(f, new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 2));

    b = EventManagerFrame::makeButton(f, "First", width);
    b->Connect("Clicked()", cls, this, "DoFirstEvent()");
    b = EventManagerFrame::makeButton(f, "Prev", width);
    b->Connect("Clicked()", cls, this, "DoPrevEvent()");

    mEventId = new TGNumberEntry(f, 0, 5, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 0, 10000);
    f->AddFrame(mEventId, new TGLayoutHints(kLHintsNormal, 10, 5, 0, 0));
    mEventId->Connect("ValueSet(Long_t)", cls, this, "DoSetEvent()");
    TGLabel* infoLabel = new TGLabel(f);
    f->AddFrame(infoLabel, new TGLayoutHints(kLHintsNormal, 5, 10, 4, 0));

    b = EventManagerFrame::makeButton(f, "Next", width);
    b->Connect("Clicked()", cls, this, "DoNextEvent()");
    b = EventManagerFrame::makeButton(f, "Last", width);
    b->Connect("Clicked()", cls, this, "DoLastEvent()");
    b = EventManagerFrame::makeButton(f, "Screenshot", 2 * width);
    b->Connect("Clicked()", cls, this, "DoScreenshot()");
    b = EventManagerFrame::makeButton(f, "Refresh", 2 * width);
    b->Connect("Clicked()", cls, this, "DoRefresh()");

    b = EventManagerFrame::makeButton(f, "R2 Geometry", 2 * width);
    b->Connect("Clicked()", cls, this, "DoR2Geometry()");
    b = EventManagerFrame::makeButton(f, "R3 Geometry", 2 * width);
    b->Connect("Clicked()", cls, this, "DoR3Geometry()");
  }
  SetCleanup(kDeepCleanup);
  Layout();
  MapSubwindows();
  MapWindow();
}

TGTextButton* EventManagerFrame::makeButton(TGCompositeFrame* p, const char* txt,
                                            Int_t width, Int_t lo, Int_t ro, Int_t to, Int_t bo)
{
  TGTextButton* b = new TGTextButton(p, txt);

  if (width > 0) {
    b->SetWidth(width);
    b->ChangeOptions(b->GetOptions() | kFixedWidth);
  }
  p->AddFrame(b, new TGLayoutHints(kLHintsNormal, lo, ro, to, bo));
  return b;
}

void EventManagerFrame::DoFirstEvent()
{
  mEventManager->GotoEvent(0);
  mEventId->SetIntNumber(mEventManager->getCurrentEvent());
}

void EventManagerFrame::DoPrevEvent()
{
  mEventManager->PrevEvent();
  mEventId->SetIntNumber(mEventManager->getCurrentEvent());
}

void EventManagerFrame::DoNextEvent()
{
  mEventManager->NextEvent();
  mEventId->SetIntNumber(mEventManager->getCurrentEvent());
}

void EventManagerFrame::DoLastEvent()
{
  mEventManager->GotoEvent(-1); /// -1 means last available
  mEventId->SetIntNumber(mEventManager->getCurrentEvent());
}

void EventManagerFrame::DoSetEvent()
{
  mEventManager->GotoEvent(mEventId->GetNumberEntry()->GetIntNumber());
}

void EventManagerFrame::DoScreenshot()
{
}

void EventManagerFrame::DoRefresh()
{
  refresh(false);
}

void EventManagerFrame::DoR2Geometry()
{
}

void EventManagerFrame::DoR3Geometry()
{
}

void EventManagerFrame::setupGeometry()
{
  // read path to geometry files from config file
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);

  // get geometry from Geometry Manager and register in multiview
  auto multiView = MultiView::getInstance();

  for (auto detName : gVisualisationGroupName) {
    if (settings.GetValue((detName + ".draw").c_str(), false)) {
      if (detName == "TPC" || detName == "MCH" || detName == "MTR" || detName == "MID" || detName == "MFT" || detName == "AD0" || detName == "FMD") { // don't load MUON+MFT and AD and standard TPC to R-Phi view
        multiView->drawGeometryForDetector(detName, true, false);
      } else if (detName == "RPH") { // special TPC geom from R-Phi view
        multiView->drawGeometryForDetector(detName, false, true, false);
      } else { // default
        multiView->drawGeometryForDetector(detName);
      }
    }
  }
}

void EventManagerFrame::setupCamera()
{
  // move and rotate sub-views
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);

  // read settings from config file
  const double angleHorizontal = settings.GetValue("camera.3D.rotation.horizontal", -0.4);
  const double angleVertical = settings.GetValue("camera.3D.rotation.vertical", 1.0);

  double zoom[MultiView::NumberOfViews];
  zoom[MultiView::View3d] = settings.GetValue("camera.3D.zoom", 1.0);
  zoom[MultiView::ViewRphi] = settings.GetValue("camera.R-Phi.zoom", 1.0);
  zoom[MultiView::ViewZrho] = settings.GetValue("camera.Rho-Z.zoom", 1.0);

  // get necessary elements of the multiview and set camera position
  auto multiView = MultiView::getInstance();

  for (int viewIter = 0; viewIter < MultiView::NumberOfViews; ++viewIter) {
    TGLViewer* glView = multiView->getView(static_cast<MultiView::EViews>(viewIter))->GetGLViewer();
    glView->CurrentCamera().Reset();

    if (viewIter == 0) {
      glView->CurrentCamera().RotateRad(angleHorizontal, angleVertical);
    }
    glView->CurrentCamera().Dolly(zoom[viewIter], kFALSE, kTRUE);
  }
}

void EventManagerFrame::setupBackground()
{
  // get viewers of multiview and change color to the value from config file
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);
  Color_t col = settings.GetValue("background.color", 1);

  for (int viewIter = 0; viewIter < MultiView::NumberOfViews; ++viewIter) {
    TEveViewer* view = MultiView::getInstance()->getView(static_cast<MultiView::EViews>(viewIter));

    view->GetGLViewer()->SetClearColor(col);
  }
}

void EventManagerFrame::refresh(bool firstTime)
{
  MultiView::getInstance()->destroyAllGeometries();
  MultiView::getInstance()->destroyAllEvents();

  setupGeometry();

  // MUST have Redraw3D(true) at application startup
  // so as R-Phi and Rho-Z scenes geometries are displayed properly
  // MUST have FullRedraw3D(true) for correct camera angle as well
  if (firstTime) {
    gEve->FullRedraw3D(true);
  }

  setupCamera();

  setupBackground();

  mEventManager->GotoEvent(mEventManager->getCurrentEvent());
  gSystem->ProcessEvents();
  gEve->Redraw3D();
}

} // namespace event_visualisation
} // namespace o2
