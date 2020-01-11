// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "Framework/InputSpec.h"
#include "Framework/CallbackService.h"
#include "Framework/ControlService.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/ParallelContext.h"
#include "Framework/runDataProcessing.h"
#include "Framework/InputRecord.h"
#include "Framework/Logger.h"
#include "Framework/AnalysisDataModel.h"

#include "EventVisualisationView/EventManager.h"
#include "EventVisualisationView/Initializer.h"
#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationBase/DataSourceOnline.h"
#include "EventVisualisationDetectors/DataReaderAOD.h"

#include <TApplication.h>
#include <TEveBrowser.h>
#include <TEveManager.h>
#include <TEnv.h>
#include <TList.h>

#include <chrono>
#include <iostream>
#include <thread>

using namespace o2::framework;
using namespace o2::event_visualisation;
using DataHeader = o2::header::DataHeader;

WorkflowSpec defineDataProcessing(ConfigContext const&) {
  return WorkflowSpec
  {
    {
      "trackDisplay",
      Inputs{
        {"tracks", "AOD", "TRACKPAR"}},
      Outputs{},
      AlgorithmSpec{
        adaptStateful(
          [](CallbackService& callbacks) {
            LOGF(info, "Welcome in O2 event visualisation tool");
            // No ITS, no TPC, only AOD
            o2::event_visualisation::Options options { false, false, true };
            TEnv settings;
            ConfigurationManager::getInstance().getConfig(settings);

            std::array<const char*, 7> keys = {
              "Gui.DefaultFont", "Gui.MenuFont", "Gui.MenuHiFont",
              "Gui.DocFixedFont", "Gui.DocPropFont", "Gui.IconFont", "Gui.StatusFont"
            };
            for (const auto& key : keys) {
              if (settings.Defined(key)) {
                gEnv->SetValue(key, settings.GetValue(key, ""));
              }
            }

            // Create ROOT application environment
            TApplication* app = new TApplication("o2eve", nullptr, nullptr);
            app->Connect("TEveBrowser", "CloseWindow()", "TApplication", app, "Terminate()");

            LOGF(info, "Initializing TEveManager");
            if (!TEveManager::Create(kTRUE, "FI")) {
              LOGF(fatal, "Could not create TEveManager!");
            }

            // Initialize o2 Event Visualisation internals
            Initializer::setup(options, EDataSource::SourceOnline);

            // Launch ROOT application in a separate thread
            // FIXME: Can it be made non-blocking?
            auto rootTAppThread = std::make_shared<std::thread>([] {
              LOGF(info, "Starting the application...");
              // Start the application
              gApplication->Run(kTRUE);

              // The thread continue here only after user closes the event display window

              LOGF(info, "Closing the application...");
              // Terminate application
              TEveManager::Terminate();
              gApplication->Terminate(0);
              LOGF(info, "ROOT app thread finishing...");
            });

            return adaptStateless([rootTAppThread](InputRecord& inputs, ControlService& control) {
              auto input = inputs.get<TableConsumer>("tracks");

              o2::aod::Tracks myTracks{input->asArrowTable()};

              int maxEventID = 0;
              TList* onlineTracks = new TList();
              for (auto& track : myTracks) {
                if (track.collisionId() > maxEventID) {
                  maxEventID = track.collisionId();
                }
                onlineTracks->Add(new AODTrack(track.collisionId(), track.x(), track.alpha(), track.y(), track.z(), track.snp(), track.tgl(), track.signed1Pt(), 0, 0));
              }

              LOGF(info, "Number of tracks in the workflow: %d", onlineTracks->GetEntries());

              auto& eventManager = EventManager::getInstance();
              DataSourceOnline* dataSource = (DataSourceOnline*)eventManager.getDataSource();
              dataSource->setEventCount(maxEventID + 1, EVisualisationGroup::AOD);
              dataSource->setOnlineEventData(onlineTracks, EVisualisationGroup::AOD, EVisualisationDataType::Tracks);

              // Wait till event display terminates properly
              rootTAppThread->join();

              control.readyToQuit(QuitRequest::All);
            });
          })
      }
    }
  };
}
