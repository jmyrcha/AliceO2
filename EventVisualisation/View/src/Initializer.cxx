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
/// \file    Initializer.cxx
/// \author  Jeremi Niedziela
///

#include "EventVisualisationView/Initializer.h"

#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationBase/EventManager.h"
#include "EventVisualisationBase/GeometryManager.h"
#include "EventVisualisationView/MultiView.h"
#include "EventVisualisationBase/VisualisationConstants.h"
#include "EventVisualisationView/EventManagerFrame.h"
#include "EventVisualisationBase/DataSourceOffline.h"
#include "EventVisualisationDetectors/DataSourceOfflineITS.h"
#include "EventVisualisationDetectors/DataInterpreterITS.h"

#include <TGTab.h>
#include <TEnv.h>
#include <TEveBrowser.h>
#include <TEveManager.h>
#include <TRegexp.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TROOT.h>
#include <TEveWindowManager.h>
#include <iostream>
#include <TFile.h>
using namespace std;


namespace o2  {
namespace event_visualisation {

Initializer::Initializer(EventManager::EDataSource defaultDataSource)
{
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);
  
  const bool fullscreen      = settings.GetValue("fullscreen.mode",false);       // hide left and bottom tabs
  const string ocdbStorage   = settings.GetValue("OCDB.default.path","local://$ALICE_ROOT/OCDB");// default path to OCDB
  cout << "Initializer -- OCDB path:" << ocdbStorage << endl;
  
  auto &eventManager = EventManager::getInstance();
  //eventManager.setDataSourceType(defaultDataSource);
  eventManager.setCdbPath(ocdbStorage);

//  DataSourceOffline *ds = new DataSourceOfflineVSD();
//DataSourceOffline *ds = new DataSourceOffline();
//  ds->Open("AliVSD.root");
  //DataSourceOffline *ds = new DataSourceOfflineITS();
  //ds->Open("o2trac_its.root");
  //eventManager.setDataSource(ds);

  eventManager.setDataSourceType(EventManager::EDataSource::SourceOffline);
  eventManager.setDataSourcePath("o2trac_its.root");
  eventManager.Open();

  //gEve->AddEvent(&eventManager);
  
  setupGeometry();
  gSystem->ProcessEvents();
  gEve->Redraw3D(kTRUE);
  
  setupBackground();
  
  // Setup windows size, fullscreen and focus
  TEveBrowser *browser = gEve->GetBrowser();
  browser->GetTabRight()->SetTab(1);
  browser->MoveResize(0, 0, gClient->GetDisplayWidth(),gClient->GetDisplayHeight() - 32);

  browser->StartEmbedding(TRootBrowser::kBottom);
  EventManagerFrame *frame =new EventManagerFrame(eventManager);
  browser->StopEmbedding("EventCtrl Balbinka");

  if(fullscreen){
    ((TGWindow*)gEve->GetBrowser()->GetTabLeft()->GetParent())->Resize(1,0);
    ((TGWindow*)gEve->GetBrowser()->GetTabBottom()->GetParent())->Resize(0,1);
  }
  gEve->GetBrowser()->Layout();
  gSystem->ProcessEvents();
  
  setupCamera();
  
  // Temporary:
  // For the time being we draw single random event on startup.
  // Later this will be triggered by button, and finally moved to configuration.
  gEve->AddEvent(&EventManager::getInstance());
  //MultiView::getInstance()->drawRandomEvent();
  //gEve->AddEvent(&eventManager);

  // AOD is dummy here, as there is no separate type for tracks
  MultiView::getInstance()->drawITSEvent(EDataType::AOD);
  frame->DoFirstEvent();
}

Initializer::~Initializer() = default;

void Initializer::setupGeometry()
{
  o2::base::GeometryManager::loadGeometry("O2geometry.root", "FAIRGeom");
  auto gman = o2::its::GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::utils::bit2Mask(o2::TransformType::T2L, o2::TransformType::T2GRot,
                                            o2::TransformType::L2G));

//  // read path to geometry files from config file
//  TEnv settings;
//  ConfigurationManager::getInstance().getConfig(settings);
//
//  // get geometry from Geometry Manager and register in multiview
//  auto multiView = MultiView::getInstance();
//
//  for(int iDet=0;iDet<NvisualisationGroups;++iDet){
//    EVisualisationGroup det = static_cast<EVisualisationGroup>(iDet);
//    string detName = gVisualisationGroupName[det];
//    if(settings.GetValue((detName+".draw").c_str(), false))
//    {
//      if(   detName=="TPC" || detName=="MCH" || detName=="MTR"
//         || detName=="MID" || detName=="MFT" || detName=="AD0"
//         || detName=="FMD"){// don't load MUON+MFT and AD and standard TPC to R-Phi view
//
//        multiView->drawGeometryForDetector(detName,true,false);
//      }
//      else if(detName=="RPH"){// special TPC geom from R-Phi view
//
//        multiView->drawGeometryForDetector(detName,false,true,false);
//      }
//      else{// default
//        multiView->drawGeometryForDetector(detName);
//      }
//    }
//  }
}
 
void Initializer::setupCamera()
{
  // move and rotate sub-views
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);
  
  // read settings from config file
  const double angleHorizontal = settings.GetValue("camera.3D.rotation.horizontal",-0.4);
  const double angleVertical   = settings.GetValue("camera.3D.rotation.vertical",1.0);
  
  double zoom[MultiView::NumberOfViews];
  zoom[MultiView::View3d]   = settings.GetValue("camera.3D.zoom",1.0);
  zoom[MultiView::ViewRphi] = settings.GetValue("camera.R-Phi.zoom",1.0);
  zoom[MultiView::ViewZrho] = settings.GetValue("camera.Rho-Z.zoom",1.0);
  
  // get necessary elements of the multiview and set camera position
  auto multiView = MultiView::getInstance();
  
  for(int viewIter=0;viewIter<MultiView::NumberOfViews;++viewIter){
    TGLViewer *glView = multiView->getView(static_cast<MultiView::EViews>(viewIter))->GetGLViewer();
    glView->CurrentCamera().Reset();
    
    if(viewIter==0){
        glView->CurrentCamera().RotateRad(angleHorizontal, angleVertical);
    }
    glView->CurrentCamera().Dolly(zoom[viewIter], kFALSE, kTRUE);
  }
}

void Initializer::setupBackground()
{
  // get viewers of multiview and change color to the value from config file
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);
  Color_t col = settings.GetValue("background.color", 1);

  for(int viewIter=0; viewIter<MultiView::NumberOfViews;++viewIter){
    TEveViewer *view = MultiView::getInstance()->getView(static_cast<MultiView::EViews>(viewIter));
    view->GetGLViewer()->SetClearColor(col);
  }
}

}
}
