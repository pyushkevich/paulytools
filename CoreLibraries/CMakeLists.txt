PROJECT(MYLIBS)

# Set a single directory for outputting all libraries
SET (LIBRARY_OUTPUT_PATH ${MYLIBS_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all libraries.")

# Set up the subdirectories as the include directories
SET(MYLIBS_INCLUDE_DIRS 
    ${MYLIBS_SOURCE_DIR}
    ${MYLIBS_SOURCE_DIR}/MReps2d
    ${MYLIBS_SOURCE_DIR}/Matrix
    ${MYLIBS_SOURCE_DIR}/Optima
    ${MYLIBS_SOURCE_DIR}/Registry
)

# Say which directories to include in this project
INCLUDE_DIRECTORIES(${MYLIBS_INCLUDE_DIRS})

# Include some sub-directories
SUBDIRS(
  Matrix
  Optima
  Registry
)

OPTION(USE_INTRINSICS "Use Intel intrinsics?  This option lets you use mreps and creps" ON)
IF(USE_INTRINSICS)
  SUBDIRS(
    MReps2d
  )
ENDIF(USE_INTRINSICS)

IF(USE_FLTK)
  SUBDIRS(
    FLTKViewers
  )
ENDIF(USE_FLTK)

# Pass on the library and include paths to config file in binary directory
#CONFIGURE_FILE(
#	"${MYLIBS_SOURCE_DIR}/Config.cmake.in"
#	"${MYLIBS_BINARY_DIR}/Config.cmake"
#	IMMEDIATE @ONLY
#)


