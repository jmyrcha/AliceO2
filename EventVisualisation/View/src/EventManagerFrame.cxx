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
/// \file    EventManagerFrame.cxx
/// \brief   GUI (bottom buttons) for visualisation
/// \author  julian.myrcha@cern.ch
/// \author  p.nowakowski@cern.ch
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationBase/GeometryManager.h"
#include "EventVisualisationView/EventManagerFrame.h"
#include "EventVisualisationView/EventManager.h"
#include "EventVisualisationView/MultiView.h"

#include "FairLogger.h"

#include <TEveManager.h>
#include <TSystem.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>
#include <TEveBrowser.h>
#include <TASImage.h>
#include <TGFileDialog.h>
#include <TTimeStamp.h>

#include <Rtypes.h>

#include <cstring>

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
    int width = 50;
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
                                            int width, int lo, int ro, int to, int bo)
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
  mEventId->SetIntNumber(mEventManager->getCurrentEventNumber());
}

void EventManagerFrame::DoPrevEvent()
{
  mEventManager->PrevEvent();
  mEventId->SetIntNumber(mEventManager->getCurrentEventNumber());
}

void EventManagerFrame::DoNextEvent()
{
  mEventManager->NextEvent();
  mEventId->SetIntNumber(mEventManager->getCurrentEventNumber());
}

void EventManagerFrame::DoLastEvent()
{
  mEventManager->GotoEvent(-1); /// -1 means last available
  mEventId->SetIntNumber(mEventManager->getCurrentEventNumber());
}

void EventManagerFrame::DoSetEvent()
{
  mEventManager->GotoEvent(mEventId->GetNumberEntry()->GetIntNumber());
}

void EventManagerFrame::DoScreenshot()
{
  saveScreenshot();
}

void EventManagerFrame::DoRefresh()
{
  refresh(false);
}

void EventManagerFrame::DoR2Geometry()
{
  auto& geomMan = GeometryManager::getInstance();
  geomMan.setR2Geometry(true);
  refresh(false);
}

void EventManagerFrame::DoR3Geometry()
{
  auto& geomMan = GeometryManager::getInstance();
  geomMan.setR2Geometry(false);
  refresh(false);
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

  mEventManager->GotoEvent(mEventManager->getCurrentEventNumber());
  gSystem->ProcessEvents();
  gEve->Redraw3D();
}

std::tuple<int, int, bool, bool, bool, const char*, const char*> EventManagerFrame::getScreenshotOptions()
{
  // 1,440 x 0,900 (2K)
  // 3,840 × 2,160 (4K)
  // 7,680 × 4,320 (8K)
  // 15,360 × 8,640 (16K)
  int width = 3840;
  int height = 2160;

  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);

  bool logo = settings.GetValue("screenshot.logo.draw", true);
  bool info = settings.GetValue("screenshot.info.draw", true);
  bool projections = settings.GetValue("screenshot.projections.draw", true);
  const char* energyLabel = settings.GetValue("screenshot.force.energy", "Energy: unknown");
  const char* systemLabel = settings.GetValue("screenshot.force.system", "Colliding system: unknown");
  if (std::strcmp(energyLabel, "") == 0)
    energyLabel = "Energy: unknown";
  if (std::strcmp(systemLabel, "") == 0)
    systemLabel = "Colliding system: unknown";

  return std::make_tuple(width, height, logo, info, projections, energyLabel, systemLabel);
}

void EventManagerFrame::drawScreenshotLogo(TASImage* compositeImg, int width)
{
  TASImage* aliceLogo = new TASImage(
    Form("%s/share/EventVisualisation/resources/alice_logo_big.png", gSystem->Getenv("O2_ROOT")));
  if (aliceLogo) {
    double ratio = 1434.0 / 1939.0;
    aliceLogo->Scale(0.08 * width, 0.08 * width / ratio);
    compositeImg->Merge(aliceLogo, "alphablend", 20, 20);
    delete aliceLogo;
  }
}

void EventManagerFrame::drawScreenshotInfo(TASImage* compositeImg, const char* energyLabel, const char* systemLabel, int height)
{
  auto& currentEvent = mEventManager->getCurrentEvent();
  TTimeStamp ts(currentEvent.getTimeStamp());
  const char* runNumber = Form("Run: %d", currentEvent.getRunNumber());
  const char* timeStamp = Form("Timestamp: %s (UTC)", ts.AsString("s"));

  const char* system;
  if (currentEvent.getCollidingSystem() == "") {
    system = systemLabel;
  } else {
    system = Form("Colliding system: %s", currentEvent.getCollidingSystem().data());
  }

  const char* energy;
  if (currentEvent.getEnergy() >= 0.0000001) {
    energy = Form("Energy: %.0f TeV", currentEvent.getEnergy() / 1000);
  } else {
    energy = energyLabel;
  }

  int fontSize = 0.015 * height;
  compositeImg->BeginPaint();
  compositeImg->DrawText(10, height - 25 - 4 * fontSize, runNumber, fontSize, "#BBBBBB", "FreeSansBold.otf");
  compositeImg->DrawText(10, height - 20 - 3 * fontSize, timeStamp, fontSize, "#BBBBBB", "FreeSansBold.otf");
  compositeImg->DrawText(10, height - 15 - 2 * fontSize, system, fontSize, "#BBBBBB", "FreeSansBold.otf");
  compositeImg->DrawText(10, height - 10 - 1 * fontSize, energy, fontSize, "#BBBBBB", "FreeSansBold.otf");
  compositeImg->EndPaint();
}

void EventManagerFrame::saveScreenshotToFile(TASImage* compositeImg)
{
  TGFileInfo fileinfo;
  const char* filetypes[] = {"All types", "*", 0, 0};
  fileinfo.fFileTypes = filetypes;
  fileinfo.fIniDir = StrDup(".");
  new TGFileDialog(gClient->GetDefaultRoot(), gClient->GetDefaultRoot(), kFDOpen, &fileinfo);
  if (!fileinfo.fFilename) {
    LOG(WARNING) << "SaveScreenshot: Couldn't get image path from dialog window!!!";
    return;
  }
  compositeImg->WriteImage(fileinfo.fFilename);
  LOG(INFO) << "Image saved to: " << fileinfo.fFilename;
}

void EventManagerFrame::saveScreenshot()
{
  auto [width, height, logo, info, projections, energyLabel, systemLabel] = getScreenshotOptions();

  gEve->GetBrowser()->RaiseWindow();
  gEve->DoRedraw3D();

  TEveViewerList* viewers = gEve->GetViewers();
  int viewersCount = viewers->NumChildren() - 1;

  TASImage* compositeImg = new TASImage(width, height);

  // 3D view size
  int width3DView = projections ? TMath::FloorNint(2.0 * width / 3.0) : width; // width of the 3D view
  int height3DView = height;                                                   // height of the 3D view
  float aspectRatio = (float)width3DView / (float)height3DView;                // 3D View aspect ratio

  // Children view size
  int widthChildView = TMath::FloorNint((float)width / 3.0);
  int heightChildView = TMath::FloorNint((float)height / viewersCount);
  float childAspectRatio = (float)widthChildView / (float)heightChildView;

  int index = 0;       // iteration counter
  int x = width3DView; // x position of the child view
  int y = 0;           // y position of the child view

  for (TEveElement::List_i i = viewers->BeginChildren(); i != viewers->EndChildren(); i++) {
    TEveViewer* view = ((TEveViewer*)*i);

    // Save OpenGL view in file and read it back using BB (missing method in Root returning TASImage)
    TString viewFilename = Form("view-%d.png", index);
    view->GetGLViewer()->SavePictureUsingBB(viewFilename);
    TASImage* viewImg = new TASImage(viewFilename);

    // Second option is to use FBO instead of BB
    // This improves the quality of pictures in some specific cases
    // but is causes a bug (moving mouse over views makes them disappear
    // on new event being loaded)

    int currentHeight, currentWidth;
    float currentAspectRatio;
    int currentX, currentY;

    if (index == 0) {
      currentHeight = height3DView;
      currentWidth = width3DView;
      currentAspectRatio = aspectRatio;
      currentX = 0;
      currentY = 0;
    } else {
      currentHeight = heightChildView;
      currentWidth = widthChildView;
      currentAspectRatio = childAspectRatio;
      currentX = x;
      currentY = y;
    }

    // For FBO version
    //TASImage* viewImg = (TASImage*)view->GetGLViewer()->GetPictureUsingFBO(currentWidth, currentHeight);

    if (viewImg) {
      if (viewImg->GetWidth() < currentAspectRatio * viewImg->GetHeight()) {
        viewImg->Crop(0, (viewImg->GetHeight() - viewImg->GetWidth() / currentAspectRatio), viewImg->GetWidth(),
                      viewImg->GetWidth() / currentAspectRatio);
      } else {
        viewImg->Crop((viewImg->GetWidth() - viewImg->GetHeight() * currentAspectRatio), 0, viewImg->GetHeight() * currentAspectRatio,
                      viewImg->GetHeight());
      }
      viewImg->Scale(currentWidth, currentHeight);
      viewImg->CopyArea(compositeImg, 0, 0, viewImg->GetWidth(), viewImg->GetHeight(), currentX, currentY);

      if (index != 0) {
        compositeImg->DrawRectangle(currentX, currentY, currentWidth, currentHeight, "#C0C0C0"); // draw a border around child views
      }
      delete viewImg;
    }

    if (index != 0) { // skip 3D View
      y += heightChildView;
    }

    if (!projections) {
      break;
    }

    index++;
  }

  // Draw ALICE Logo
  if (logo) {
    drawScreenshotLogo(compositeImg, width);
  }

  // Draw info
  if (info) {
    drawScreenshotInfo(compositeImg, energyLabel, systemLabel, height);
  }

  // Save screenshot to file
  saveScreenshotToFile(compositeImg);

  delete compositeImg;
}

} // namespace event_visualisation
} // namespace o2
