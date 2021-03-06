# DO NOT EDIT THIS FILE! THIS FILE IS AUTOMATICALLY GENERATED BY xml2cmake
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8 FATAL_ERROR)
SET(VISIT_INCLUDE_DIR "/usr/local/visit/2.10.2/2.10.2/linux-x86_64/include")
SET(VISIT_LIBRARY_DIR "/usr/local/visit/2.10.2/2.10.2/linux-x86_64/lib")
SET(VISIT_BINARY_DIR "/usr/local/visit/2.10.2/2.10.2/linux-x86_64/bin")
SET(VISIT_ARCHIVE_DIR "/usr/local/visit/2.10.2/2.10.2/linux-x86_64/archives")
SET(VISIT_PLUGIN_DIR "/home/jiming/.visit/2.10.2/linux-x86_64/plugins")
INCLUDE("/usr/local/visit/2.10.2/2.10.2/linux-x86_64/include/PluginVsInstall.cmake")
INCLUDE("/usr/local/visit/2.10.2/2.10.2/linux-x86_64/include/VisItLibraryDependencies.cmake")

PROJECT(ATHENA_database)

SET(COMMON_SOURCES
ATHENAPluginInfo.C
ATHENACommonPluginInfo.C
)

SET(LIBI_SOURCES 
ATHENAPluginInfo.C
)

SET(LIBM_SOURCES
ATHENAMDServerPluginInfo.C
${COMMON_SOURCES}
avtATHENAFileFormat.C
)

SET(LIBE_SOURCES
ATHENAEnginePluginInfo.C
${COMMON_SOURCES}
avtATHENAFileFormat.C
)

INCLUDE_DIRECTORIES(
${CMAKE_CURRENT_SOURCE_DIR}
${HDF5_INCLUDE_DIR}
${VISIT_COMMON_INCLUDES}
${VISIT_INCLUDE_DIR}/visit/avt/DBAtts/MetaData
${VISIT_INCLUDE_DIR}/visit/avt/DBAtts/SIL
${VISIT_INCLUDE_DIR}/visit/avt/Database/Database
${VISIT_INCLUDE_DIR}/visit/avt/Database/Formats
${VISIT_INCLUDE_DIR}/visit/avt/Database/Ghost
${VISIT_INCLUDE_DIR}/visit/avt/FileWriter
${VISIT_INCLUDE_DIR}/visit/avt/Filters
${VISIT_INCLUDE_DIR}/visit/avt/MIR/Base
${VISIT_INCLUDE_DIR}/visit/avt/MIR/Tet
${VISIT_INCLUDE_DIR}/visit/avt/MIR/Zoo
${VISIT_INCLUDE_DIR}/visit/avt/Math
${VISIT_INCLUDE_DIR}/visit/avt/Pipeline/AbstractFilters
${VISIT_INCLUDE_DIR}/visit/avt/Pipeline/Data
${VISIT_INCLUDE_DIR}/visit/avt/Pipeline/Pipeline
${VISIT_INCLUDE_DIR}/visit/avt/Pipeline/Sinks
${VISIT_INCLUDE_DIR}/visit/avt/Pipeline/Sources
${VISIT_INCLUDE_DIR}/visit/avt/VisWindow/VisWindow
${VISIT_INCLUDE_DIR}/visit/visit_vtk/full
${VISIT_INCLUDE_DIR}/visit/visit_vtk/lightweight
${EAVL_INCLUDE_DIR} 
${VTK_INCLUDE_DIRS} 
)


LINK_DIRECTORIES(${VISIT_LIBRARY_DIR} ${EAVL_LIBRARY_DIR} ${VTK_LIBRARY_DIRS} ${HDF5_LIBRARY_DIR} )

ADD_LIBRARY(IATHENADatabase ${LIBI_SOURCES})
TARGET_LINK_LIBRARIES(IATHENADatabase visitcommon)
SET(INSTALLTARGETS IATHENADatabase)

IF(NOT VISIT_ENGINE_ONLY AND NOT VISIT_DBIO_ONLY)
    ADD_LIBRARY(MATHENADatabase ${LIBM_SOURCES}    )
    TARGET_LINK_LIBRARIES(MATHENADatabase visitcommon avtdbatts avtdatabase_ser ${HDF5_LIB} )
    ADD_TARGET_DEFINITIONS(MATHENADatabase MDSERVER)
    SET(INSTALLTARGETS ${INSTALLTARGETS} MATHENADatabase)
ENDIF(NOT VISIT_ENGINE_ONLY AND NOT VISIT_DBIO_ONLY)

ADD_LIBRARY(EATHENADatabase_ser ${LIBE_SOURCES})
TARGET_LINK_LIBRARIES(EATHENADatabase_ser visitcommon avtdatabase_ser avtpipeline_ser ${HDF5_LIB} )
ADD_TARGET_DEFINITIONS(EATHENADatabase_ser ENGINE)
SET(INSTALLTARGETS ${INSTALLTARGETS} EATHENADatabase_ser)

IF(VISIT_PARALLEL)
    ADD_PARALLEL_LIBRARY(EATHENADatabase_par ${LIBE_SOURCES})
    TARGET_LINK_LIBRARIES(EATHENADatabase_par visitcommon avtdatabase_par avtpipeline_par ${HDF5_LIB} )
    ADD_TARGET_DEFINITIONS(EATHENADatabase_par ENGINE)
    SET(INSTALLTARGETS ${INSTALLTARGETS} EATHENADatabase_par)
ENDIF(VISIT_PARALLEL)

VISIT_INSTALL_DATABASE_PLUGINS(${INSTALLTARGETS})
VISIT_PLUGIN_TARGET_RTOD(databases ${INSTALLTARGETS})

