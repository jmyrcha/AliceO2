// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CaloMatrix.h
/// \brief temporary, to provide alignment matrix for event display
/// \author Maja Kabus

#ifndef ALICE_O2_EVENTVISUALISATION_DETECTORS_CALOMATRIX_H
#define ALICE_O2_EVENTVISUALISATION_DETECTORS_CALOMATRIX_H

/// This file contains EMCAL and PHOS alignment matrices from sample Run 2 AliESD.root.
/// It is used to properly visualise calorimeter cells in event display.
/// This file will be removed once it will be possible to obtain matrices
/// from OCDB / AOD / other valid sources

#include <TGeoMatrix.h>
#include <TError.h>

namespace o2
{
namespace event_visualisation
{

class CaloMatrix
{
 public:
  static TGeoHMatrix* getEMCALMatrix(int superModule)
  {
    if (superModule < 0 || superModule > 19) {
      Error("CaloMatrix::getEMCALMatrix", "Super module number out of bounds");
    }
    auto matrix = new TGeoHMatrix;
    matrix->SetRotation(mEMCALRotation[superModule]);
    matrix->SetTranslation(mEMCALTranslation[superModule]);
    return matrix;
  }

  static TGeoHMatrix* getPHOSMatrix(int module)
  {
    if (module < 0 || module > 4) {
      Error("CaloMatrix::getPHOSMatrix", "Module number out of bounds");
    }
    auto matrix = new TGeoHMatrix;
    matrix->SetRotation(mPHOSRotation[module]);
    matrix->SetTranslation(mPHOSTranslation[module]);
    return matrix;
  }

 private:
  static const double constexpr mEMCALRotation[20][9] = {
    { -0.006944, -0.999971, -0.003101, 0.999974, -0.006950, 0.001959, -0.001980, -0.003087, 0.999993 },
    { -0.006949, 0.999971, 0.003076, 0.999974, 0.006955, -0.001781, -0.001802, 0.003064, -0.999994 },
    { -0.341741, -0.939786, -0.003909, 0.939790, -0.341749, 0.001629, -0.002867, -0.003117, 0.999991 },
    { -0.342967, 0.939340, 0.003726, 0.939343, 0.342974, -0.001620, -0.002800, 0.002945, -0.999992 },
    { -0.640317, -0.768103, -0.003573, 0.768111, -0.640313, -0.002273, -0.000541, -0.004200, 0.999991 },
    { -0.637305, 0.770597, 0.004713, 0.770580, 0.637322, -0.005201, -0.007011, 0.000317, -0.999975 },
    { -0.865058, -0.501664, -0.002602, 0.501663, -0.865062, 0.001139, -0.002822, -0.000321, 0.999996 },
    { -0.856010, 0.516941, 0.004225, 0.516930, 0.856020, -0.003575, -0.005465, -0.000876, -0.999985 },
    { -0.983552, -0.180595, -0.003415, 0.180578, -0.983549, 0.004878, -0.004239, 0.004181, 0.999982 },
    { -0.982680, 0.185294, 0.002600, 0.185293, 0.982683, -0.000325, -0.002615, 0.000162, -0.999997 },
    { -0.987763, 0.155962, 0.000000, -0.155962, -0.987763, 0.000000, 0.000000, 0.000000, 1.000000 },
    { -0.987762, -0.155962, 0.001624, -0.155963, 0.987763, -0.000325, -0.001554, -0.000574, -0.999999 },
    { -0.004792, 0.999989, 0.000071, -0.999861, -0.004792, 0.015979, 0.015979, 0.000005, 0.999872 },
    { -0.000903, -1.000000, -0.000237, -0.999916, 0.000900, 0.012896, -0.012896, 0.000249, -0.999917 },
    { 0.335634, 0.941971, -0.006353, -0.941818, 0.335694, 0.016974, 0.018122, 0.000287, 0.999836 },
    { 0.340000, -0.940417, -0.004078, -0.940275, -0.340020, 0.016429, -0.016836, -0.001752, -0.999857 },
    { 0.641400, 0.767153, -0.009103, -0.767022, 0.641461, 0.014343, 0.016842, -0.002217, 0.999856 },
    { 0.639954, -0.768362, -0.008914, -0.768196, -0.640006, 0.016358, -0.018273, -0.003621, -0.999826 },
    { 0.875383, 0.483430, -0.000324, -0.483411, 0.875354, 0.008257, 0.004275, -0.007071, 0.999966 },
    { 0.871620, -0.490182, 0.000487, -0.490167, -0.871583, 0.008905, -0.003941, -0.008001, -0.999960 }
  };

  static const double constexpr mEMCALTranslation[20][3] = {
    { 1.103595, 446.179186, 176.156236 },
    { 2.184570, 445.524752, -173.845817 },
    { -153.438665, 418.122535, 176.985028 },
    { -152.120705, 417.547276, -173.025378 },
    { -289.189499, 340.175281, 176.609885 },
    { -287.694534, 339.699917, -173.890416 },
    { -389.876653, 220.892504, 177.093304 },
    { -388.546732, 220.295691, -173.268943 },
    { -443.438688, 73.950905, 177.319832 },
    { -442.373124, 73.110819, -173.013460 },
    { -452.004041, -35.679580, 176.963780 },
    { -451.719786, -35.736435, -173.072967 },
    { 1.297238, -451.169769, 221.666956 },
    { 0.527137, -451.522544, -221.016420 },
    { 154.636036, -424.769450, 220.994434 },
    { 154.843257, -424.798970, -221.042832 },
    { 291.562450, -346.984624, 220.697797 },
    { 290.694360, -347.212468, -221.155547 },
    { 366.183681, -270.685708, 175.669604 },
    { 366.269398, -270.673096, -174.801908 }
  };

  static const double constexpr mPHOSRotation[4][9] = {
    { 0.767338, -0.004312, 0.641229, 0.641188, -0.007942, -0.767343, 0.008401, 0.999959, -0.003329 },
    { 0.940798, 0.002136, 0.338961, 0.338959, 0.001227, -0.940800, -0.002425, 0.999997, 0.000430 },
    { 0.999976, -0.000756, -0.006876, -0.006877, -0.001172, -0.999976, 0.000748, 0.999999, -0.001177 },
    { 0.934365, -0.002092, -0.356311, -0.356314, -0.000801, -0.934366, 0.001669, 0.999997, -0.001494 }
  };

  static const double constexpr mPHOSTranslation[4][3] = {
    { 313.381746, -372.024764, -2.354855 },
    { 165.078518, -454.148853, -2.801675 },
    { -0.083626, -482.172358, -1.906440 },
    { -167.233337, -460.901723, -1.821245 }
  };
};

} // namespace event_visualisation
} // namespace o2

#endif //ALICE_O2_EVENTVISUALISATION_DETECTORS_CALOMATRIX_H
