PROJECT(CMREPS)

# Option to build the new PDE code
OPTION(BUILD_PDE "Build the PDE CM-Reps code" OFF)
VARIABLE_REQUIRES(BUILD_PDE DO_BUILD_PDE USE_SWIG) 

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

# Set the compiler flags on GCC
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse -mmmx")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

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

# Include the various directories
INCLUDE_DIRECTORIES(
  ${CMREPS_SOURCE_DIR}/Code 
  ${CMREPS_SOURCE_DIR}/Code/PDE
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

# Building the PDE executable
IF(BUILD_PDE)

  # Sources for the PDE executable
  SET(PDE_COMMON_SRCS 
    Code/PDE/MedialPDESolver.cxx 
    Code/PDE/FourierSurface.cxx
    Code/PDE/ScriptInterface.cxx
    Code/Danielsson.cpp
    Code/DrawFillInside.cpp
    Code/DrawTriangles.cpp
    Code/glengine.cpp
    Code/imaging.cpp
    Code/ImageCubeITK.cxx
    Code/mspline.cpp
    Code/fracviewer.cpp)

  # Apply SWIG to the .i file
  SWIG_WRAP_PYTHON_SHADOW(
    ${CMREPS_SOURCE_DIR}/Code/PDE/medialpde.i
    ${CMREPS_BINARY_DIR}/Code/PDE/ScriptInterface_wrap.cxx
    ${CMREPS_BINARY_DIR}
    medialpde)

  # Add python to the includes
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

  # The sources for the script-wrapped library
  SET(PDE_LIB_SRCS ${PDE_COMMON_SRCS}
    ${CMREPS_BINARY_DIR}/Code/PDE/ScriptInterface_wrap.cxx)

  # Create the executable
  ADD_EXECUTABLE(testpde ${PDE_COMMON_SRCS} Code/PDE/MedialPDETest.cxx)

  # Create the shared library
  ADD_LIBRARY("_medialpdemodule" SHARED ${PDE_LIB_SRCS})
  SET_TARGET_PROPERTIES("_medialpdemodule" PROPERTIES PREFIX "")

  # Find and read the MKL package
  INCLUDE(${PAULY_SOURCE_DIR}/CMake/FindMKL.cmake)

  # Load Intel MKL libraries, etc
  FIND_LIBRARY(PARDISO_LIBRARY NAMES pardiso DOC "Location of the PARDISO library")
  FIND_LIBRARY(INTEL_FORTRAN_LIB ifcore /opt/intel_fc_80/lib DOC "Location of Intel FC library")

  # Set the library set
  SET(PDE_LIBS
    ${PARDISO_LIBRARY} 
    mkl_lapack mkl_def 
    ${INTEL_FORTRAN_LIB} 
    ITKNumerics
    ITKIO
    ITKCommon
    ITKBasicFilters
    ${GLUT_LIBRARY})

  # Create the executable
  TARGET_LINK_LIBRARIES(testpde ${PDE_LIBS})
  TARGET_LINK_LIBRARIES("_medialpdemodule" ${PDE_LIBS})
  
  # Copy the script files
  SET(SCRIPT_FILES example.py)
  FOREACH(FILE ${SCRIPT_FILES})  
 	CONFIGURE_FILE(
 		${CMREPS_SOURCE_DIR}/Code/PDE/Python/${FILE}
 		${CMREPS_BINARY_DIR}/${FILE} COPY_ONLY)
  ENDFOREACH(FILE ${SCRIPT_FILES})  
  
ENDIF(BUILD_PDE)

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
