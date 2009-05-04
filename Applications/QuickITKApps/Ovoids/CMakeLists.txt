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

LINK_LIBRARIES(
  ${ITK_LIBRARIES}
  ${VTK_COMMON_LIBS}
  ${PAULY_CORE_LIBS}
  ${FLTK_GL_LIBRARY}
  ${FLTK_BASE_LIBRARY}
)

FLTK_WRAP_UI(ovoids ui.fl)
ADD_EXECUTABLE(ovoids 
  ${OVOIDS_SRC}
  ${ovoids_FLTK_UI_SRCS}
  main.cpp)

ADD_EXECUTABLE(bettermesh 
  ${OVOIDS_SRC} 
  bettermesh.cpp)

IF(WIN32)
  LINK_LIBRARIES(wsock32 comctl32)
ENDIF(WIN32)
