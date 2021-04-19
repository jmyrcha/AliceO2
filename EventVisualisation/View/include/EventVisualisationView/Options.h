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
/// \file    Options.cxx
/// \author  Julian Myrcha

#ifndef ALICE_O2_EVENTVISUALISATION_VIEW_OPTIONS_H
#define ALICE_O2_EVENTVISUALISATION_VIEW_OPTIONS_H

#include <string>

namespace o2
{
namespace event_visualisation
{

class Options
{
 private:
  static Options instance;
  bool mRandomTracks;    // -r
  bool mVsd;             // -v
  bool mItc;             // -i
  bool mJSON;            // -j
  std::string mFileName; // -f 'data.root'
  bool saveTojson(std::string filename);
  bool readFromjson(std::string filename);

 public:
  static Options* Instance() { return &instance; }
  std::string printOptions();
  std::string usage();
  bool processCommandLine(int argc, char* argv[]);

  bool randomTracks() { return this->mRandomTracks; }
  std::string fileName() { return this->mFileName; }
  bool vsd() { return this->mVsd; }
  bool itc() { return this->mItc; }
  bool json() { return this->mJSON; }
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_VIEW_OPTIONS_H
