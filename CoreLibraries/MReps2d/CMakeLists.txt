SET(CREPS_SRCS
  bspline.cpp
  interp2D.cpp
  mreps2D.cpp
  phquintic.cpp
  phspline.cpp
)

IF(USE_OPENGL)
  SUBDIRS(GLTools)
ENDIF(USE_OPENGL)

ADD_LIBRARY(MReps2D ${CREPS_SRCS})
