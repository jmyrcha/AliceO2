# Copyright CERN and copyright holders of ALICE O2. This software is distributed
# under the terms of the GNU General Public License v3 (GPL Version 3), copied
# verbatim in the file "COPYING".
#
# See http://alice-o2.web.cern.ch/license for full licensing information.
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

o2_add_library(EventVisualisationWorkflow
               SOURCES src/SimpleWorkflow.cxx
                       src/ReaderWorkflow.cxx
                       src/TPCWorkflow.cxx
                       src/WriterWorkflow.cxx
                       src/PublisherSpec.cxx
               PUBLIC_LINK_LIBRARIES O2::DetectorsBase
                                     O2::Framework
                                     O2::DPLUtils
                                     O2::Algorithm
                                     O2::DataFormatsTPC
                                     O2::TPCBase
)


o2_add_executable(o2-simple-workflow
                  SOURCES src/SimpleWorkflow.cxx
                  PUBLIC_LINK_LIBRARIES O2::EventVisualisationWorkflow
                  )

o2_add_executable(o2-eve-workflow
                  SOURCES src/eve-workflow.cxx
                  PUBLIC_LINK_LIBRARIES O2::EventVisualisationWorkflow
                  )

o2_add_executable(o2-eve-reader-workflow
                  SOURCES src/eve-reader-workflow.cxx
                  PUBLIC_LINK_LIBRARIES O2::EventVisualisationWorkflow
                  )

o2_add_executable(o2-eve-writer-workflow
                  SOURCES src/eve-writer-workflow.cxx
                  PUBLIC_LINK_LIBRARIES O2::EventVisualisationWorkflow
                  )

