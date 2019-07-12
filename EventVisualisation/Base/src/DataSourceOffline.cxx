//
//  Created by jmy on 26.02.19.
//

#include "EventVisualisationBase/DataSourceOffline.h"

//#include <TSystem.h>
//#include <TEveTreeTools.h>
//#include <TEveTrack.h>
//#include <TEveManager.h>
//#include <TFile.h>
//#include <TPRegexp.h>
//#include <TEveTrackPropagator.h>
//#include <TEveViewer.h>
//#include <TEveEventManager.h>
//#include <TEveVSD.h>
//#include <TVector3.h>

namespace o2 {
namespace event_visualisation {

DataSourceOffline::DataSourceOffline():
  fgRawFileName("raw.root"),
  fgClustersFileName("clusters.root"),
  fgTracksFileName("tracks.root")
{

}

void DataSourceOffline::Open(TString ESDFileName)
{
  if(fIsOpen)
  {
    std::cout << "Files already opened!" << std::endl;
    return;
  }

  this->fgESDFileName = ESDFileName;
  OpenRawFile();
  OpenTracksFile();
  fIsOpen = kTRUE;
}


}
}
