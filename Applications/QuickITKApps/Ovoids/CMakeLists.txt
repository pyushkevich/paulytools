SET(OVOIDS_SRC
  blobmodl.cpp
  blobmodl.h
  fracviewer.cpp
  fracviewer.h
  glcode.cpp
  glcode.h
  implicit.c
  implicit.h
  main.cpp
  plib.cpp
  plib.h
)

INCLUDE_DIRECTORIES(
  ${PAULY_SOURCE_DIR}/Applications/QuickITKApps/Ovoids
)

ADD_EXECUTABLE(ovoids ${OVOIDS_SRC})

FLTK_WRAP_UI(ovoids ui.fl)

LINK_LIBRARIES(
  ITKNumerics
  ${PAULY_CORE_LIBS}
  ${FLTK_GL_LIBRARY}
  ${FLTK_BASE_LIBRARY}
)


