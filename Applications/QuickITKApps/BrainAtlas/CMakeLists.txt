INCLUDE_DIRECTORIES(
  ${PAULY_SOURCE_DIR}/Applications/QuickITKApps/BrainAtlas
)

OPTION(USE_ITK_FFTW "Set this option if ITK has FFTW support" OFF)

# Tracer application GUI
SET(TRACER_SRCS
  BinaryHeap.h
  ShortestPath.cxx
  ShortestPath.h
  TracerCurves.cxx
  TracerCurves.h
  TracerData.cxx
  TracerData.h
  TracerMain.cxx
  TracerMainWindow.h
  TracerMainWindow.cxx
  TracerUserInterfaceBase.h
  TracerUserInterfaceLogic.cxx
  TracerUserInterfaceLogic.h
  Trackball.cxx
  Trackball.h
  VTKMeshShortestDistance.cxx
  VTKMeshShortestDistance.h
  VTKMeshVoronoiDiagram.cxx
  VTKMeshVoronoiDiagram.h
)

ADD_EXECUTABLE("binary2mesh" BinaryToMesh.cxx)
ADD_EXECUTABLE("viewmesh" MeshViewer.cxx)
ADD_EXECUTABLE("binary2metis" ImageToMetis.cxx)
ADD_EXECUTABLE("maskimage" MaskImage.cxx)
ADD_EXECUTABLE(braintracer ${TRACER_SRCS})

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

LINK_LIBRARIES(
  ITKCommon 
  ITKBasicFilters 
  ITKIO 
  ITKNumerics
  ITKStatistics
  ${PAULY_CORE_LIBS}
  vtkCommon 
  vtkFiltering
  vtkImaging
  vtkIO 
  vtkGraphics 
  vtkHybrid
  vtkRendering
  vtkPatented
  ${FLTK_LIBRARY}
)

SET(TRACER_FLUIDS TracerUserInterface.fl)

FLTK_WRAP_UI(braintracer ${TRACER_FLUIDS})

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
  
CONFIGURE_FILE(
  ${PAULY_SOURCE_DIR}/Applications/QuickITKApps/BrainAtlas/Utilities/InnoSetup/InstallBrainTracer.iss.in
  ${PAULY_BINARY_DIR}/Applications/QuickITKApps/BrainAtlas/Utilities/InnoSetup/InstallBrainTracer.iss)

FOREACH(SCRIPT ${BATLAS_SCRIPTS})
  CONFIGURE_FILE(
    ${PAULY_SOURCE_DIR}/Applications/QuickITKApps/BrainAtlas/Scripts/${SCRIPT}.in
    ${PAULY_BINARY_DIR}/Applications/QuickITKApps/BrainAtlas/${SCRIPT}
  )
ENDFOREACH(SCRIPT)
