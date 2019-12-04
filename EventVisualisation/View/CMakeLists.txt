# Copyright CERN and copyright holders of ALICE O2. This software is distributed
# under the terms of the GNU General Public License v3 (GPL Version 3), copied
# verbatim in the file "COPYING".
#
# See http://alice-o2.web.cern.ch/license for full licensing information.
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

o2_add_library(EventVisualisationView
               SOURCES src/MultiView.cxx src/Initializer.cxx src/EventManagerFrame.cxx src/EventManager.cxx
               PUBLIC_LINK_LIBRARIES O2::EventVisualisationBase
                                     O2::EventVisualisationDetectors
                                     FairLogger::FairLogger)

o2_target_root_dictionary(EventVisualisationView
                          HEADERS include/EventVisualisationView/MultiView.h
                                  include/EventVisualisationView/Initializer.h
                                  include/EventVisualisationView/EventManager.h
                                  include/EventVisualisationView/EventManagerFrame.h
                          LINKDEF src/EventVisualisationViewLinkDef.h)

o2_add_executable(eve
                  SOURCES src/main.cxx
                  PUBLIC_LINK_LIBRARIES O2::EventVisualisationView
                  )

o2_add_test(View
        COMPONENT_NAME eve
        PUBLIC_LINK_LIBRARIES O2::EventVisualisationView
        SOURCES test/testEventVisualisationView.cxx
        LABELS eve)

# file(GLOB EVE_TESTS "/home/maja/CERN/alice/sw/BUILD/O2-latest-maja-unit-tests/O2/stage/tests/o2-test-eve-*")
#set(EVE_TESTS "/home/maja/CERN/alice/sw/BUILD/O2-latest-maja-unit-tests/O2/stage/tests/o2-test-eve-View")
#file(COPY "/home/maja/CERN/alice/sw/BUILD/O2-latest-maja-unit-tests/O2/stage/tests/o2-test-eve-View" DESTINATION "/home/maja/CERN/new_alieve/o2eve/root-files/tests")
