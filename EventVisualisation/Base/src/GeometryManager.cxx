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
/// \file    GeometryManager.cxx
/// \author  Jeremi Niedziela

#include "EventVisualisationBase/GeometryManager.h"

#include "EventVisualisationBase/ConfigurationManager.h"

#include <TFile.h>
#include <TGLViewer.h>
#include <TEnv.h>
#include <TEveGeoShapeExtract.h>
#include <TEveManager.h>
#include <TEveProjectionManager.h>
#include <TSystem.h>

#include <iostream>

using namespace std;

namespace o2  {
namespace event_visualisation {
  
GeometryManager& GeometryManager::getInstance()
{
  static GeometryManager instance;
  return instance;
}

TEveGeoShape* GeometryManager::getGeometryForDetector(string detectorName, bool run2)
{
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings, run2);

  // read geometry path from config file
  string geomPath = settings.GetValue("simple.geom.path","");

  // TODO:
  // we need a way to set O2 installation path here
  //
  const string o2basePathSpecifier = "${ALICE_ROOT}";
  const string o2basePath = "";//= gSystem->Getenv("ALICE_ROOT");
  const size_t o2pos = geomPath.find(o2basePathSpecifier);

  if(o2pos != string::npos){
    geomPath.replace(o2pos,o2pos+o2basePathSpecifier.size(),o2basePath);
  }

  const string runPath = run2 ? "/run2" : "/O2";

  // load ROOT file with geometry
  TFile *f = TFile::Open(Form("%s%s/simple_geom_%s.root", geomPath.c_str(), runPath.c_str(),detectorName.c_str()));
  //cout << "GeometryManager::GetSimpleGeom opening geometry for: " << detectorName << endl;
  if(!f){
    cout<<"GeometryManager::GetSimpleGeom -- no file with geometry found for: "<<detectorName<<"!"<<endl;
    return nullptr;
  }
  
  TEveGeoShapeExtract *geomShapeExtract = static_cast<TEveGeoShapeExtract*>(f->Get(detectorName.c_str()));
  TEveGeoShape *geomShape = TEveGeoShape::ImportShapeExtract(geomShapeExtract);
  f->Close();
  geomShape->SetName(detectorName.c_str());
  
  // tricks for different R-Phi geom of TPC:
  if(detectorName=="RPH"){  // use all other parameters of regular TPC geom
    detectorName = "TPC";
  }

  // prepare geometry to be drawn including all children
  drawDeep(geomShape, run2,
           settings.GetValue(Form("%s.color",detectorName.c_str()),-1),
           settings.GetValue(Form("%s.trans",detectorName.c_str()),-1),
           settings.GetValue(Form("%s.line.color",detectorName.c_str()),-1));

  gEve->GetDefaultGLViewer()->UpdateScene();

  return geomShape;
}

void GeometryManager::drawDeep(TEveGeoShape *geomShape, bool run2, Color_t color, Char_t transparency, Color_t lineColor)
{
  if(geomShape->HasChildren()){
    geomShape->SetRnrSelf(false);
    
    if(run2 && strcmp(geomShape->GetElementName(),"TPC_Drift_1")==0){// hack for TPC drift chamber
      geomShape->SetRnrSelf(kTRUE);
      if(color>=0) geomShape->SetMainColor(color);
      if(lineColor>=0){
        geomShape->SetLineColor(lineColor);
        geomShape->SetLineWidth(0.1);
        geomShape->SetDrawFrame(true);
      }
      else{
        geomShape->SetDrawFrame(false);
      }
      if(transparency>=0){
        geomShape->SetMainTransparency(transparency);
      }
    }
    
    for(TEveElement::List_i i = geomShape->BeginChildren(); i != geomShape->EndChildren(); ++i){
      drawDeep(static_cast<TEveGeoShape*>(*i), run2, color, transparency, lineColor);
    }
  }
  else
  {
    geomShape->SetRnrSelf(true);
    if(color>=0) geomShape->SetMainColor(color);
    if(lineColor>=0){
      geomShape->SetLineColor(lineColor);
      geomShape->SetLineWidth(0.1);
      geomShape->SetDrawFrame(true);
    }
    else{
      geomShape->SetDrawFrame(false);
    }
    if(transparency>=0){
      geomShape->SetMainTransparency(transparency);
    }
      
    if(run2 && strcmp(geomShape->GetElementName(),"PHOS_5")==0){// hack for PHOS module which is not installed
      geomShape->SetRnrSelf(false);
    }
  }
}
  
}
}
