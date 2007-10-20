INCLUDE_DIRECTORIES(
  ${PAULY_SOURCE_DIR}/Applications/QuickITKApps/BrainAtlas
  ${PAULY_SOURCE_DIR}/CoreLibraries/ITKExtensions/WellComposed
)

OPTION(USE_ITK_FFTW "Set this option if ITK has FFTW support" OFF)

ADD_EXECUTABLE(binary2mesh BinaryToMesh.cxx)
ADD_EXECUTABLE(viewmesh MeshViewer.cxx)
ADD_EXECUTABLE(maskimage MaskImage.cxx)
ADD_EXECUTABLE(threshimage ThresholdImage.cxx)
ADD_EXECUTABLE(dumpmeshpoints DumpMeshPoints.cxx)

# Use FFTW from ITK if requested
IF(USE_ITK_FFTW)
  SET(USE_FFTW 1)
  INCLUDE(${ITK_SOURCE_DIR}/CMake/FindFFTW.cmake)
  
  IF(FFTW_FOUND)
	INCLUDE_DIRECTORIES(${FFTW_INCLUDE_PATH})
	ADD_EXECUTABLE(maskair MaskAir.cxx)
	TARGET_LINK_LIBRARIES(maskair ${FFTW_LIB})
  ENDIF(FFTW_FOUND)
ENDIF(USE_ITK_FFTW)

SET(COMMON_LIBS
  ITKCommon 
  ITKBasicFilters 
  ITKIO 
  ITKNumerics
  ITKStatistics
  ITKVoxBoIO
  ${PAULY_CORE_LIBS}
  ${VTK_COMMON_LIBS}
  ${FLTK_LIBRARY}
)

TARGET_LINK_LIBRARIES(binary2mesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(viewmesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(maskimage ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(threshimage ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(dumpmeshpoints ${COMMON_LIBS})

# Add the option to use OpenMesh
OPTION(USE_OPENMESH "Use OpenMesh meshing library" OFF)
IF(USE_OPENMESH)
  INCLUDE(${PAULY_CMAKE_DIR}/FindOpenMesh.cmake)
ENDIF(USE_OPENMESH)

# Add the openmesh-dependent code
IF(OPENMESH_FOUND)
  INCLUDE_DIRECTORIES(${OPENMESH_INCLUDE_DIR})
  ADD_EXECUTABLE(sulcusfilter FilterSulci.cxx implicit.c)
  TARGET_LINK_LIBRARIES(sulcusfilter ${OPENMESH_LIBRARY})
ENDIF(OPENMESH_FOUND)

# Configure scripts
FIND_PATH(BATLAS_DATA_ROOT "atlas02" "~/data2004" "c:/data2004")

SET(BATLAS_SCRIPTS
  makemesh.sh
  makewhitemesh.sh
  metispartition.sh
)
  
FOREACH(SCRIPT ${BATLAS_SCRIPTS})
  CONFIGURE_FILE(
    ${PAULY_SOURCE_DIR}/Applications/QuickITKApps/BrainAtlas/Scripts/${SCRIPT}.in
    ${PAULY_BINARY_DIR}/Applications/QuickITKApps/BrainAtlas/${SCRIPT}
  )
ENDFOREACH(SCRIPT)

INSTALL_TARGETS(
  "/bin"
  maskimage threshimage binary2mesh viewmesh)
