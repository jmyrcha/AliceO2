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
/// \author  Julian Myrcha
///

#include "EventVisualisationBase/GeometryManager.h"
#include "EventVisualisationBase/ConfigurationManager.h"

#include "FairLogger.h"

#include <TFile.h>
#include <TGLViewer.h>
#include <TEnv.h>
#include <TEveGeoShapeExtract.h>
#include <TEveManager.h>
#include <TEveProjectionManager.h>
#include <TSystem.h>

namespace o2
{
namespace event_visualisation
{

GeometryManager* GeometryManager::sInstance = nullptr;

GeometryManager::GeometryManager()
{
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);
  mR2Geometry = settings.GetValue("simple.geom.default", "R2") == "R2";
}

GeometryManager& GeometryManager::getInstance()
{
  if (sInstance == nullptr) {
    sInstance = new GeometryManager();
  }
  return *sInstance;
}

TEveGeoShape* GeometryManager::getGeometryForDetector(std::string detectorName)
{
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);

  // read geometry path from config file
  std::string geomPath = settings.GetValue(this->mR2Geometry ? "simple.geom.R2.path" : "simple.geom.R3.path", "");

  // load ROOT file with geometry
  TFile* f = TFile::Open(Form("%s/simple_geom_%s.root", geomPath.c_str(), detectorName.c_str()));
  if (!f) {
    LOG(ERROR) << "GeometryManager::GetSimpleGeom -- no file with geometry found for: " << detectorName << "!";
    return nullptr;
  }
  LOG(INFO) << "GeometryManager::GetSimpleGeom for: " << detectorName << " from " << Form("%s/simple_geom_%s.root", geomPath.c_str(), detectorName.c_str());

  TEveGeoShapeExtract* geomShapeExtract = static_cast<TEveGeoShapeExtract*>(f->Get(detectorName.c_str()));
  TEveGeoShape* geomShape = TEveGeoShape::ImportShapeExtract(geomShapeExtract);
  f->Close();

  geomShape->SetName(detectorName.c_str());
  // tricks for different R-Phi geom of TPC:
  if (detectorName == "RPH") { // use all other parameters of regular TPC geom
    detectorName = "TPC";
  }

  // prepare geometry to be drawn including all children
  drawDeep(geomShape,
           settings.GetValue(Form("%s.color", detectorName.c_str()), -1),
           settings.GetValue(Form("%s.trans", detectorName.c_str()), -1),
           settings.GetValue(Form("%s.line.color", detectorName.c_str()), -1));

  gEve->GetDefaultGLViewer()->UpdateScene();

  return geomShape;
}

void GeometryManager::drawDeep(TEveGeoShape* geomShape, Color_t color, char transparency, Color_t lineColor)
{
  if (geomShape->HasChildren()) {
    geomShape->SetRnrSelf(false);

    if (mR2Geometry && strcmp(geomShape->GetElementName(), "TPC_Drift_1") == 0) { // hack for TPC drift chamber
      geomShape->SetRnrSelf(kTRUE);
      if (color >= 0)
        geomShape->SetMainColor(color);
      if (lineColor >= 0) {
        geomShape->SetLineColor(lineColor);
        geomShape->SetLineWidth(0.1);
        geomShape->SetDrawFrame(true);
      } else {
        geomShape->SetDrawFrame(false);
      }
      if (transparency >= 0) {
        geomShape->SetMainTransparency(transparency);
      }
    }

    for (TEveElement::List_i i = geomShape->BeginChildren(); i != geomShape->EndChildren(); ++i) {
      drawDeep(static_cast<TEveGeoShape*>(*i), color, transparency, lineColor);
    }
  } else {
    geomShape->SetRnrSelf(true);
    if (color >= 0)
      geomShape->SetMainColor(color);
    if (lineColor >= 0) {
      geomShape->SetLineColor(lineColor);
      geomShape->SetLineWidth(0.1);
      geomShape->SetDrawFrame(true);
    } else {
      geomShape->SetDrawFrame(false);
    }
    if (transparency >= 0) {
      geomShape->SetMainTransparency(transparency);
    }

    if (mR2Geometry && strcmp(geomShape->GetElementName(), "PHOS_5") == 0) { // hack for PHOS module which is not installed
      geomShape->SetRnrSelf(false);
    }
  }
}

} // namespace event_visualisation
} // namespace o2
