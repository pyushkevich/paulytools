SET(OVOIDS_SRC
  blobmodl.cpp
  blobmodl.h
  fracviewer.cpp
  fracviewer.h
  glcode.cpp
  glcode.h
  implicit.c
  implicit.h
  imagesurface.cpp
  imagesurface.h
  plib.cpp
  plib.h
)

INCLUDE_DIRECTORIES(
  ${PAULY_SOURCE_DIR}/Applications/QuickITKApps/Ovoids
)

ADD_EXECUTABLE(ovoids ${OVOIDS_SRC} main.cpp)
ADD_EXECUTABLE(bettermesh ${OVOIDS_SRC} bettermesh.cpp)

FLTK_WRAP_UI(ovoids ui.fl)

LINK_LIBRARIES(
  ITKCommon 
  ITKBasicFilters 
  ITKIO 
  ITKNumerics
  vtkCommon
  ${PAULY_CORE_LIBS}
  ${FLTK_LIBRARY}
)

IF(WIN32)
  LINK_LIBRARIES(wsock32 comctl32)
ENDIF(WIN32)
