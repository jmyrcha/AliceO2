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
/// \file    EventManagerFrame.h
/// \brief   GUI (bottom buttons) for visualisation
/// \author  julian.myrcha@cern.ch
/// \author  p.nowakowski@cern.ch
/// \author  Maja Kabus <maja.kabus@cern.ch>
///

#ifndef ALICE_O2_EVENTVISUALISATION_VIEW_EVENTMANAGERFRAME_H
#define ALICE_O2_EVENTVISUALISATION_VIEW_EVENTMANAGERFRAME_H

#include "EventVisualisationView/EventManager.h"

#include <TGMdiMainFrame.h>

#include <tuple>

class TGTextButton;
class TGCompositeFrame;
class TGNumberEntry;
class TGLabel;

namespace o2
{
namespace event_visualisation
{

class EventManagerFrame : public TGMainFrame
{
 private:
  static TGTextButton* makeButton(TGCompositeFrame* p, const char* txt, Int_t width = 0,
                                  Int_t lo = 0, Int_t ro = 0, Int_t to = 0, Int_t bo = 0);

 protected:
  o2::event_visualisation::EventManager* mEventManager; // Model object.
  TGNumberEntry* mEventId;                              // Display/edit current event id

 public:
  EventManagerFrame(o2::event_visualisation::EventManager& eventManager);
  ~EventManagerFrame() override = default;
  ClassDefOverride(EventManagerFrame, 0); // GUI window for EventManager.

  /// Loads geometry for all detectors
  void setupGeometry();
  /// Sets up background color
  void setupBackground();
  /// Sets up camera position
  void setupCamera();
  /// Clears and draws everything anew, used also by Initializer on startup
  void refresh(bool firstTime);
  /// Saves a screenshot
  void saveScreenshot();
  std::tuple<int, int, bool, bool, bool, const char*, const char*> getScreenshotOptions();

 public: // slots
  void DoFirstEvent();
  void DoPrevEvent();
  void DoNextEvent();
  void DoLastEvent();
  void DoSetEvent();
  void DoScreenshot();
  void DoRefresh();
  void DoR2Geometry();
  void DoR3Geometry();
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_VIEW_EVENTMANAGERFRAME_H
