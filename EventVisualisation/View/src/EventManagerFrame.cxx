//
// Created by jmy on 23.02.19.
//


#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>
#include <EventVisualisationView/EventManagerFrame.h>
#include <EventVisualisationView/MultiView.h>
#include <EventVisualisationBase/DataSourceOffline.h>
#include <EventVisualisationBase/DataReaderVSD.h>
#include "EventVisualisationBase/ConfigurationManager.h"
#include <Rtypes.h>
#include <iostream>

ClassImp(o2::event_visualisation::EventManagerFrame)

namespace o2 {
namespace event_visualisation {

EventManagerFrame::EventManagerFrame(o2::event_visualisation::EventManager& eventManager)
:TGMainFrame(gClient->GetRoot(), 400, 100, kVerticalFrame) {
    std::cout << "EventManagerFrame constructor" << std::endl;

    fM = &eventManager;

    const TString cls("o2::event_visualisation::EventManagerFrame");
    TGTextButton *b = 0;
    TGHorizontalFrame *f = new TGHorizontalFrame(this);
    {
        Int_t width = 50;
        this->AddFrame(f, new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 2));

        fFirstEvent = b = EventManagerFrame::makeButton(f, "First", width);
        b->Connect("Clicked()", cls, this, "DoFirstEvent()");
        fPrevEvent = b = EventManagerFrame::makeButton(f, "Prev", width);
        b->Connect("Clicked()", cls, this, "DoPrevEvent()");

        fEventId = new TGNumberEntry(f, 0, 5, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
                                     TGNumberFormat::kNELLimitMinMax, 0, 10000);
        f->AddFrame(fEventId, new TGLayoutHints(kLHintsNormal, 10, 5, 0, 0));
        fEventId->Connect("ValueSet(Long_t)", cls, this, "DoSetEvent()");
        fInfoLabel = new TGLabel(f);
        f->AddFrame(fInfoLabel, new TGLayoutHints(kLHintsNormal, 5, 10, 4, 0));

        fNextEvent = b = EventManagerFrame::makeButton(f, "Next", width);
        b->Connect("Clicked()", cls, this, "DoNextEvent()");
        fLastEvent = b = EventManagerFrame::makeButton(f, "Last", width);
        b->Connect("Clicked()", cls, this, "DoLastEvent()");
        fScreenshot = b = EventManagerFrame::makeButton(f, "Screenshot", 2 * width);
        b->Connect("Clicked()", cls, this, "DoScreenshot()");
        fOldGeom = b = EventManagerFrame::makeButton(f, "Run 2 geometry", 2 * width);
        b->Connect("Clicked()", cls, this, "DoOldGeometry()");
        fNewGeom = b = EventManagerFrame::makeButton(f, "O2 geometry", 2 * width);
        b->Connect("Clicked()", cls, this, "DoNewGeometry()");
    }

    SetCleanup(kDeepCleanup);

    Layout();
    MapSubwindows();
    MapWindow();
}

EventManagerFrame::~EventManagerFrame() {
    std::cout << "EventManagerFrame::~EventManagerFrame()" << std::endl;
}

TGTextButton* EventManagerFrame::makeButton(TGCompositeFrame *p, const char *txt,
        Int_t width, Int_t lo, Int_t ro, Int_t to, Int_t bo) {
    TGTextButton* b = new TGTextButton(p, txt);

    //b->SetFont("-adobe-helvetica-bold-r-*-*-48-*-*-*-*-*-iso8859-1");

    if (width > 0) {
        b->SetWidth(width);
        b->ChangeOptions(b->GetOptions() | kFixedWidth);
    }
    p->AddFrame(b, new TGLayoutHints(kLHintsNormal, lo,ro,to,bo));
    return b;
}

void EventManagerFrame::DoFirstEvent() {
    std::cout << "DoFirstEvent: Getting instance of Multiview" << std::endl;
    auto multiView = MultiView::getInstance();
    std::cout << "DoFirstEvent: Getting EventRegistration instance" << std::endl;
    auto eventReg = EventRegistration::getInstance();
    fM->GotoEvent(0);
    std::cout << "DoFirstEvent: After fM Getting instance of Multiview" << std::endl;
    multiView = MultiView::getInstance();
    std::cout << "DoFirstEvent: After fM Getting EventRegistration instance" << std::endl;
    eventReg = EventRegistration::getInstance();
    fEventId->SetIntNumber(fM->getCurrentEvent());
//    DataReaderVSD* offline = dynamic_cast<DataReaderVSD*>(fM->getDataSource());
//    std::cout << "DoFirstEvent" << std::endl;
//    fEventId->SetIntNumber(0);
//    offline->GotoEvent(0);
}

void EventManagerFrame::DoPrevEvent() {
    fM->PrevEvent();
    fEventId->SetIntNumber(fM->getCurrentEvent());
//    std::cout << "DoPrevEvent" << std::endl;
//    fEventId->SetIntNumber(fM->gotoEvent(fEventId->GetNumber()-1));
//    DataReaderVSD* offline = dynamic_cast<DataReaderVSD*>(fM->getDataSource());
//    offline->GotoEvent(fEventId->GetNumber());
}

void EventManagerFrame::DoNextEvent() {
  std::cout << "DoNextEvent: Getting instance of Multiview" << std::endl;
  auto multiView = MultiView::getInstance();
  std::cout << "DoNextEvent: Getting EventRegistration instance" << std::endl;
  auto eventReg = EventRegistration::getInstance();
    fM->NextEvent();
    fEventId->SetIntNumber(fM->getCurrentEvent());
//    std::cout << "DoNextEvent" << std::endl;
//    fEventId->SetIntNumber(fM->gotoEvent(fEventId->GetNumber()+1));
//    DataReaderVSD* offline = dynamic_cast<DataReaderVSD*>(fM->getDataSource());
//    offline->GotoEvent(fEventId->GetNumber());
}

void EventManagerFrame::DoLastEvent() {
    fM->GotoEvent(-1);  /// -1 means last available
    fEventId->SetIntNumber(fM->getCurrentEvent());
//    std::cout << "DoLastEvent" << std::endl;
//    fEventId->SetIntNumber(fM->gotoEvent(-1));      // -1 means last available
//    DataReaderVSD* offline = dynamic_cast<DataReaderVSD*>(fM->getDataSource());
//    offline->GotoEvent(fEventId->GetNumber());
}

void EventManagerFrame::DoSetEvent() {
}

void EventManagerFrame::DoScreenshot() {

}

void EventManagerFrame::DoOldGeometry() {
  auto multiView = MultiView::getInstance();
  std::cout << "DoOldGeometry: Getting EventRegistration instance" << std::endl;
  auto eventReg = EventRegistration::getInstance();
  multiView->destroyAllGeometries();
  setupGeometry(true);
}

void EventManagerFrame::DoNewGeometry() {
  auto multiView = MultiView::getInstance();
  std::cout << "DoNewGeometry: Getting EventRegistration instance" << std::endl;
  auto eventReg = EventRegistration::getInstance();
  multiView->destroyAllGeometries();
  setupGeometry(false);
}

void EventManagerFrame::setupGeometry(bool oldGeom)
{
  // read path to geometry files from config file
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);

  // get geometry from Geometry Manager and register in multiview
  auto multiView = MultiView::getInstance();
  std::cout << "EventManagerFrame::setupGeometry: Getting EventRegistration instance" << std::endl;
  auto eventReg = EventRegistration::getInstance();

  for(int iDet=0;iDet<NvisualisationGroups;++iDet){
    EVisualisationGroup det = static_cast<EVisualisationGroup>(iDet);
    std::string detName = gVisualisationGroupName[det];
    if(settings.GetValue((detName+".draw").c_str(), false))
    {
      if(   detName=="TPC" || detName=="MCH" || detName=="MTR"
            || detName=="MID" || detName=="MFT" || detName=="AD0"
            || detName=="FMD"){// don't load MUON+MFT and AD and standard TPC to R-Phi view

        multiView->drawGeometryForDetector(detName, oldGeom, true,false);
      }
      else if(detName=="RPH"){// special TPC geom from R-Phi view

        multiView->drawGeometryForDetector(detName, oldGeom, false,true,false);
      }
      else{// default
        multiView->drawGeometryForDetector(detName, oldGeom);
      }
    }
  }
}

}
}
