INCLUDE_DIRECTORIES(
  ${PAULY_SOURCE_DIR}/Applications/QuickITKApps/BrainAtlas
)


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
)

ADD_EXECUTABLE("binary2mesh" BinaryToMesh.cxx)
ADD_EXECUTABLE("viewmesh" MeshViewer.cxx)
ADD_EXECUTABLE("binary2metis" ImageToMetis.cxx)
ADD_EXECUTABLE("metis2mask" MetisToImage.cxx)
ADD_EXECUTABLE("sulcusfilter" FilterSulci.cxx)
ADD_EXECUTABLE("maskimage" MaskImage.cxx)
ADD_EXECUTABLE(braintracer ${TRACER_SRCS})

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

# Add the option to use Metis
OPTION(USE_METIS "Use METIS graph partitioning library" OFF)
IF(USE_METIS)
  FIND_PATH(METIS_ROOT_DIR "Lib/metis.h")
  IF(METIS_ROOT_DIR)
    FIND_LIBRARY(METIS_LIBRARY "metis" ${METIS_ROOT_DIR} ${METIS_ROOT_DIR}/Release ${METIS_ROOT_DIR}/Debug)
    FIND_PATH(METIS_INCLUDE_DIR "metis.h" ${METIS_ROOT_DIR}/Lib)

    INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR})
    LINK_LIBRARIES(${METIS_LIBRARY})
  ENDIF(METIS_ROOT_DIR)
ENDIF(USE_METIS)

# Configure scripts
FIND_PATH(BATLAS_DATA_ROOT "atlas02" "~/data2004" "c:/data2004")
FIND_PATH(METIS_BINARY_PATH "libmetis.a")

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
