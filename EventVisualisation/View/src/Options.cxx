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

#include "EventVisualisationView/Options.h"
#include "FairLogger.h"
#include <unistd.h>
#include <sstream>

using namespace std;

namespace o2
{
namespace event_visualisation
{

Options Options::instance;

std::string Options::printOptions()
{
  static const char* str[2] = {"false", "true"};
  stringstream ss;
  ss << "fileName    : " << this->fileName() << std::endl;
  ss << "randomTracks: " << str[this->randomTracks()] << std::endl;
  ss << "itc         : " << str[this->itc()] << std::endl;
  ss << "json        : " << str[this->json()] << std::endl;
  ss << "vsd         : " << str[this->vsd()] << std::endl;
  return ss.str();
  ;
}

std::string Options::usage()
{
  stringstream ss;
  ss << "usage:" << std::endl;
  ss << "\t"
     << "o2eve <options>" << std::endl;
  ss << "\t\t"
     << "where <options> are any from the following:" << std::endl;
  ss << "\t\t"
     << "-f name        name of the data file" << std::endl;
  ss << "\t\t"
     << "-i             use itc reading from files as a source" << std::endl;
  ss << "\t\t"
     << "-j             use json files as a source" << std::endl;
  ss << "\t\t"
     << "-o name        name of the options file" << std::endl;
  ss << "\t\t"
     << "-r             use random tracks" << std::endl;
  ss << "\t\t"
     << "-s             save options to o2eve.json in current folder" << std::endl;
  ss << "\t\t"
     << "-v             use vsd files as a source" << std::endl;
  ss << "\tdefault values are always taken from o2eve.json in current folder if present" << std::endl;
  return ss.str();
}
bool Options::processCommandLine(int argc, char* argv[])
{
  int opt;
  bool save = false;
  std::string optionsFileName = "o2eve.json"; // name with options to use

  // put ':' in the starting of the
  // string so that program can
  //distinguish between '?' and ':'
  while ((opt = getopt(argc, argv, ":f:ijo:rsv")) != -1) {
    switch (opt) {
      case 'f':
        this->mFileName = optarg;
        break;
      case 'i':
        this->mItc = true;
        break;
      case 'j':
        this->mItc = true;
        break;
      case 'o':
        optionsFileName = optarg;
        break;
      case 'r':
        this->mRandomTracks = true;
        break;
      case 's':
        save = true;
        break;
      case 'v':
        this->mVsd = true;
        break;
      case ':':
        LOG(ERROR) << "option needs a value: " << char(optopt);
        LOG(INFO) << usage();
        return false;
      case '?':
        LOG(ERROR) << "unknown option: " << char(optopt);
        LOG(INFO) << usage();
        return false;
    }
  }

  // optind is for the extra arguments
  // which are not parsed
  for (; optind < argc; optind++) {
    LOG(ERROR) << "extra arguments: " << argv[optind];
    LOG(INFO) << usage();
    return false;
  }

  return true;
}

} // namespace event_visualisation
} // namespace o2