PROJECT(CMREPS)

# Don't build anything unless the version of CMake is high enough.
# The matching ELSE/ENDIF should be the last lines in the file.
IF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 1.5)

# Find the ITK Package
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build CMReps without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

# Find MyLibs (for matrix code)
FIND_PATH(MYLIBS_BINARY_PATH Config.cmake)
IF(MYLIBS_BINARY_PATH) 
    INCLUDE(${MYLIBS_BINARY_PATH}/Config.cmake)
    INCLUDE_DIRECTORIES(${MYLIBS_INCLUDE_DIRS})
    LINK_DIRECTORIES(${MYLIBS_LIBRARY_DIRS})
ENDIF(MYLIBS_BINARY_PATH)

# Find SML Include path
# FIND_PATH(SML_INLCUDE_PATH smlmath.h)
# FIND_LIBRARY(SML_LIBRARY ${SML_INCLUDE_PATH}/Release)

# Find GLUT and other GL components
INCLUDE (${CMAKE_ROOT}/Modules/FindGLUT.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/FindGLU.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/FindOpenGL.cmake)

# Include the PDE code
SUBDIRS(Code/PDE)

# Include the various directories
INCLUDE_DIRECTORIES(
  ${CMREPS_SOURCE_DIR}/Code 
#  ${SML_INLCUDE_PATH} 
  ${GLUT_INCLUDE_DIR})

# Create a source listing 
SET(CMREPS_SOURCES
  Code/curveaxis.cpp
  Code/Danielsson.cpp
  Code/DrawFillInside.cpp
  Code/DrawTriangles.cpp
  Code/fracviewer.cpp
  Code/geodesic.cpp
  Code/glengine.cpp
  Code/imaging.cpp
  Code/imatch.cpp
  Code/implicit.c
  Code/main.cpp
  Code/mspline.cpp
  Code/optimization.cpp
  Code/optmodes.cpp
  Code/splgeo.cpp
  Code/uimodes.cpp
)

# Create the include files
SET(CMREPS_HEADERS
  Code/array2d.h
  Code/colors.h
  Code/Danielsson.h
  Code/DrawFillInside.h
  Code/DrawTriangles.h
  Code/fracviewer.h
  Code/geodesic.h
  Code/glengine.h
  Code/glui.h
  Code/imaging.h
  Code/imatch.h
  Code/implicit.h
  Code/mspline.h
  Code/optimization.h
  Code/smlmath.h
  Code/splgeo.h
)

# Create the executable
ADD_EXECUTABLE(cmreps ${CMREPS_SOURCES} ${CMREPS_HEADERS})
ADD_EXECUTABLE(testglut Code/testglut.cpp)

# Add the libraries
TARGET_LINK_LIBRARIES(cmreps
  ${PAULY_CORE_LIBS}
  ${GLUT_LIBRARY}
  ITKIO
  ITKNumerics
  ITKAlgorithms
  ITKCommon
  ITKBasicFilters
)

# All of this contigent on CMake version
ENDIF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 1.5)

# Copy over the icon files
#SET(CMREP_RESOURCES 
#	robot.tga 
#	robot_act.tga 
#	woild.tga 
#	woild_act.tga 
#	mspline.ini
#)
#
#FOREACH(FILE ${CMREP_RESOURCES})  
#	CONFIGURE_FILE(
#		${PAULY_SOURCE_DIR}/Applications/CMRepTool3d/Code/Data/${FILE}
#		${PAULY_BINARY_DIR}/Applications/CMRepTool3d/${FILE} COPY_ONLY)
#ENDFOREACH(FILE)
