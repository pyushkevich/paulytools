SET(CREPS_SRCS
  bspline.cpp
  interp2D.cpp
  mreps2D.cpp
  phquintic.cpp
  phspline.cpp
  bspline.h
  cboundariness.h
  crep2d.h
  interp2D.h
  mreps2D.h
  phquintic.h
  phspline.h
  SMLVec4f.h
)

IF(USE_FLTK)
  SUBDIRS(GLTools)
ENDIF(USE_FLTK)

ADD_LIBRARY(MReps2D ${CREPS_SRCS})
