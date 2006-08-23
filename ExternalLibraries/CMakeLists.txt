# Define sub-options for each package
OPTION(BUILD_SOPLEX "Build the external SoPlex library" OFF)
OPTION(BUILD_LASPACK "Build the external LasPack library" OFF)

SUBDIRS(Netlib/toms/toms611)

IF(BUILD_SOPLEX)
  SUBDIRS(SoPlex)
ENDIF(BUILD_SOPLEX)

IF(BUILD_LASPACK)
  SUBDIRS(LASPack)
ENDIF(BUILD_LASPACK)

IF(USE_FLTK)
  SUBDIRS(FLTKExtras)
ENDIF(USE_FLTK)
