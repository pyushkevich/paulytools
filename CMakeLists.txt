###########################################################################
# PaulY Tools root CMake lists file 
###########################################################################
PROJECT(PAULY)

# Define optional configuration (USE) variables
OPTION(USE_ITK "Use ITK Libraries" OFF)
OPTION(USE_VTK "Use VTK Libraries" OFF)
OPTION(USE_FLTK "Use FLTK Libraries" OFF)
OPTION(USE_OPENGL "Use OpenGL Libraries" OFF)

# Define the optional inclusion (BUILD) variables 
OPTION(BUILD_APPLICATIONS "Should I build standalone applications?" OFF)
OPTION(BUILD_EXTERNAL "Should I build external libraries?" ON)

# Define the inclusion path
SUBDIRS(
  CoreLibraries
)

IF(BUILD_EXTERNAL)
  SUBDIRS(ExternalLibraries)
ENDIF(BUILD_EXTERNAL)

IF(BUILD_APPLICATIONS)
  SUBDIRS(Applications)
ENDIF(BUILD_APPLICATIONS)

# Set up OpenGL
IF(USE_OPENGL)
  INCLUDE (${CMAKE_ROOT}/Modules/FindOpenGL.cmake)
  INCLUDE (${CMAKE_ROOT}/Modules/FindGLU.cmake)
  LINK_LIBRARIES(${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})
ENDIF(USE_OPENGL)

# Set up FLTK
IF(USE_FLTK)
  INCLUDE (${CMAKE_ROOT}/Modules/FindFLTK.cmake)
  INCLUDE_DIRECTORIES(${FLTK_INCLUDE_PATH})
ENDIF(USE_FLTK)

# Set up ITK
IF(USE_ITK)
  FIND_PACKAGE(ITK)
  IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
      "Cannot build without ITK.  Please set ITK_DIR.")
  ENDIF(ITK_FOUND)
ENDIF(USE_ITK)

# Set up VTK
IF(USE_VTK)
  FIND_PACKAGE(VTK)
  IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
  ELSE(VTK_FOUND)
    MESSAGE(FATAL_ERROR
      "Cannot build without VTK.  Please set VTK_DIR.")
  ENDIF(VTK_FOUND)
ENDIF(USE_VTK)  

# Set up the includes 
SET(PAULY_INCLUDE_DIRS
  ${PAULY_SOURCE_DIR}/CoreLibraries/
  ${PAULY_SOURCE_DIR}/CoreLibraries/BasicImaging
  ${PAULY_SOURCE_DIR}/CoreLibraries/FLTKViewers
  ${PAULY_SOURCE_DIR}/CoreLibraries/Matrix
  ${PAULY_SOURCE_DIR}/CoreLibraries/MReps2d
  ${PAULY_SOURCE_DIR}/CoreLibraries/MReps2d/GLTools
  ${PAULY_SOURCE_DIR}/CoreLibraries/Optima
  ${PAULY_SOURCE_DIR}/CoreLibraries/Registry/
  ${PAULY_SOURCE_DIR}/ExternalLibraries/FLTKExtras/flIDE-0.3/src
  ${PAULY_SOURCE_DIR}/ExternalLibraries/SoPlex/soplex-1.2.1/src
)

# Set up the libraries to link
SET(PAULY_LIBRARY_DIRS
  ${PAULY_BINARY_DIR}/CoreLibraries/bin
  ${PAULY_BINARY_DIR}/ExternalLibraries/SoPlex
  ${PAULY_BINARY_DIR}/ExternalLibraries/FLTKExtras
)

# Set up the libraries
SET(PAULY_CORE_LIBS
  BasicImaging
  Matrix
  Optima
  Registry
)
                                                                                
# Say which directories to include in this project
INCLUDE_DIRECTORIES(${PAULY_INCLUDE_DIRS})
                                                                                

