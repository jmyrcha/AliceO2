# Copyright CERN and copyright holders of ALICE O2. This software is distributed
# under the terms of the GNU General Public License v3 (GPL Version 3), copied
# verbatim in the file "COPYING".
#
# See http://alice-o2.web.cern.ch/license for full licensing information.
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

o2_add_library(EventVisualisationDetectors
               SOURCES src/DataReaderAOD.cxx
                       src/DataReaderITS.cxx
                       src/DataReaderTPC.cxx
                       src/DataInterpreterAOD.cxx
                       src/DataInterpreterITS.cxx
                       src/DataInterpreterTPC.cxx
               PUBLIC_LINK_LIBRARIES O2::EventVisualisationBase
                                     O2::DataFormatsITSMFT
                                     O2::DataFormatsITS
                                     O2::ITSBase
                                     O2::TPCBase
                                     O2::EMCALBase
                                     O2::PHOSBase
                                     O2::DataFormatsTPC
                                     O2::ReconstructionDataFormats
                                     O2::CCDB
                                     FairLogger::FairLogger
               )

o2_add_executable(fill-o2cdb-emcal-phos-matrix
        SOURCES src/FillO2CDBEMCALPHOSMatrix.cxx
        PUBLIC_LINK_LIBRARIES O2::CCDB
                              O2::DetectorsBase
                              O2::DetectorsCommonDataFormats
                              O2::EventVisualisationBase
        )

o2_add_test(Detectors
        COMPONENT_NAME eve
        PUBLIC_LINK_LIBRARIES O2::EventVisualisationDetectors
        SOURCES test/testEventVisualisationDetectors.cxx
        LABELS eve detectors)
