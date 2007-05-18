PROJECT(TRACER)

INCLUDE_DIRECTORIES(
  ${TRACER_SOURCE_DIR}/Code
)

# Create the GUI cource
SET(TRACER_FLUIDS Code/TracerUserInterface.fl)
FLTK_WRAP_UI(braintracer ${TRACER_FLUIDS})

# Tracer application GUI
SET(TRACER_SRCS
  Code/BinaryHeap.h
  Code/EventListenerModel.h
  Code/FLWidgetActivationManager.h
  Code/ShortestPath.cxx
  Code/ShortestPath.h
  Code/TracerCurves.cxx
  Code/TracerCurves.h
  Code/TracerData.cxx
  Code/TracerData.h
  Code/TracerMain.cxx
  Code/TracerMainWindow.h
  Code/TracerMainWindow.cxx
  Code/TracerUserInterfaceBase.h
  Code/TracerUserInterfaceLogic.cxx
  Code/TracerUserInterfaceLogic.h
  Code/Trackball.cxx
  Code/Trackball.h
  Code/VTKMeshShortestDistance.cxx
  Code/VTKMeshShortestDistance.h
  Code/VTKMeshVoronoiDiagram.cxx
  Code/VTKMeshVoronoiDiagram.h
  Code/VTKMeshHalfEdgeWrapper.h
  ${braintracer_FLTK_UI_SRCS}
)

ADD_EXECUTABLE(braintracer ${TRACER_SRCS})

TARGET_LINK_LIBRARIES(braintracer
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
  #  vtkPatented
  ${FLTK_LIBRARY}
)

# Set up the Windows install script
CONFIGURE_FILE(
  ${TRACER_SOURCE_DIR}/Utilities/InnoSetup/InstallBrainTracer.iss.in
  ${TRACER_BINARY_DIR}/Utilities/InnoSetup/InstallBrainTracer.iss)
