PROJECT(CMREPS)

# Set the compiler flags on GCC
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse -mmmx")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

# Find MyLibs (for matrix code)
#FIND_PATH(MYLIBS_BINARY_PATH Config.cmake)
#IF(MYLIBS_BINARY_PATH) 
#    INCLUDE(${MYLIBS_BINARY_PATH}/Config.cmake)
#    INCLUDE_DIRECTORIES(${MYLIBS_INCLUDE_DIRS})
#    LINK_DIRECTORIES(${MYLIBS_LIBRARY_DIRS})
#ENDIF(MYLIBS_BINARY_PATH)

# Find SML Include path
# FIND_PATH(SML_INLCUDE_PATH smlmath.h)
# FIND_LIBRARY(SML_LIBRARY ${SML_INCLUDE_PATH}/Release)

# Include the various directories
INCLUDE_DIRECTORIES(
  ${CMREPS_SOURCE_DIR}/Code 
  ${GLUT_INCLUDE_DIR})

# Create a source listing 
SET(CMREPS_SOURCES
  Code/curveaxis.cpp
  Code/Danielsson.cpp
  Code/DrawFillInside.cpp
  Code/DrawTriangles.cpp
  Code/main.cpp
  Code/geodesic.cpp
  Code/glengine.cpp
  Code/imaging.cpp
  Code/imatch.cpp
  Code/implicit.c
  Code/mspline.cpp
  Code/optimization.cpp
  Code/optmodes.cpp
  Code/splgeo.cpp
  Code/uimodes.cpp
  Code/fracviewer.cpp
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

TARGET_LINK_LIBRARIES(testglut ${GLUT_LIBRARY})

# Copy over the icon files
SET(CMREP_RESOURCES 
	robot.tga 
	robot_act.tga 
	woild.tga 
	woild_act.tga 
	mspline.ini
)

FOREACH(FILE ${CMREP_RESOURCES})  
	CONFIGURE_FILE(
		${CMREPS_SOURCE_DIR}/Code/Data/${FILE}
		${CMREPS_BINARY_DIR}/${FILE} COPYONLY)
ENDFOREACH(FILE)
