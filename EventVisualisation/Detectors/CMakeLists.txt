# @author  Jeremi Niedziela

set(MODULE_NAME "EventVisualisationDetectors")

O2_SETUP(NAME ${MODULE_NAME})

# Define the source and header files
set(SRCS
    src/DataInterpreterRND.cxx
    src/DataInterpreterVSD.cxx
    src/DataInterpreterITS.cxx
    src/DataSourceOfflineVSD.cxx
    src/DataSourceOfflineITS.cxx
    )

# HEADERS is not needed if we don't generate a dictionary.
set(HEADERS
    include/${MODULE_NAME}/DataInterpreterRND.h
    include/${MODULE_NAME}/DataInterpreterVSD.h
    include/${MODULE_NAME}/DataInterpreterITS.h
    include/${MODULE_NAME}/DataSourceOfflineVSD.h
    include/${MODULE_NAME}/DataSourceOfflineITS.h
    )

set(LIBRARY_NAME ${MODULE_NAME})
set(BUCKET_NAME event_visualisation_detectors_bucket)

O2_GENERATE_LIBRARY()
