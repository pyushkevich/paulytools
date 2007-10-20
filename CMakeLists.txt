###########################################################################
# PaulY Tools root CMake lists file 
###########################################################################
PROJECT(PAULY)

# Define optional configuration (USE) variables
OPTION(USE_ITK "Use ITK Libraries" OFF)
OPTION(USE_VTK "Use VTK Libraries" OFF)
OPTION(USE_FLTK "Use FLTK Libraries" OFF)
OPTION(USE_OPENGL "Use OpenGL Libraries" OFF)
OPTION(USE_GLUT "Use GLUT Library for OpenGL" OFF)
OPTION(USE_MMA "Use Mathematica MathLink Libraries" OFF)
OPTION(USE_CLUTO "Use Cluto Clustering Libraries" OFF)
OPTION(USE_SWIG "Use the SWIG wrapper" OFF)
OPTION(USE_INTEL_MKL "Use Intel's MKL Numerics Library (Alternative to Lapack)" OFF)
OPTION(USE_LAPACK "Use the Lapack Numerics Library" OFF)
OPTION(USE_PARDISO "Use the PARDISO sparse solver library" OFF)
OPTION(USE_METIS "Use the METIS Graph Partitioning Library" OFF)
OPTION(USE_DART "Whether this project should be configured to use Dart" OFF)
OPTION(USE_GSL "Use GSL (GNU Scientific Library)" OFF)

# Define the optional inclusion (BUILD) variables 
OPTION(BUILD_APPLICATIONS "Should I build standalone applications?" OFF)
OPTION(BUILD_EXTERNAL "Should I build external libraries?" OFF)
OPTION(BUILD_TESTING "Should I build the test code?" OFF)

# Set the path for the CMake modules
FIND_PATH(PAULY_CMAKE_DIR "FindMathematica.cmake" PATHS ${PAULY_SOURCE_DIR}/CMake)
MARK_AS_ADVANCED(PAULY_CMAKE_DIR)

# By default, include the core libraries
IF(USE_ITK)
  SUBDIRS( CoreLibraries )
ENDIF(USE_ITK)
  
# Optionally, include the external libraries
IF(BUILD_EXTERNAL)
  SUBDIRS(ExternalLibraries)
ENDIF(BUILD_EXTERNAL)

# Optionally, include the applications
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
# ELSE(ITK_FOUND)
#   MESSAGE(FATAL_ERROR
#     "Cannot build without ITK.  Please set ITK_DIR.")
  ENDIF(ITK_FOUND)
ENDIF(USE_ITK)

# Set up VTK
IF(USE_VTK)
  FIND_PACKAGE(VTK)
  IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})

    # Set flags depending on VTK version
    IF("${VTK_MAJOR_VERSION}" EQUAL 5)
      SET(HAVE_VTK_5 1)
      SET(VTK_COMMON_LIBS vtkCommon vtkIO vtkGraphics vtkRendering
        vtkFiltering vtkImaging vtkHybrid)
    ELSE("${VTK_MAJOR_VERSION}" EQUAL 5)
      SET(HAVE_VTK_5 0)
      SET(VTK_COMMON_LIBS vtkCommon vtkIO vtkGraphics vtkRendering vtkPatented
        vtkFiltering vtkImaging vtkHybrid)
    ENDIF("${VTK_MAJOR_VERSION}" EQUAL 5)

  ELSE(VTK_FOUND)
    MESSAGE(FATAL_ERROR
      "Cannot build without VTK.  Please set VTK_DIR.")
  ENDIF(VTK_FOUND)
ENDIF(USE_VTK)  

# Set up Mathematica
IF(USE_MMA)
  INCLUDE(${PAULY_CMAKE_DIR}/FindMathematica.cmake)
  INCLUDE_DIRECTORIES(${MMA_INCLUDE_DIR})
ENDIF(USE_MMA)

# CLUTO Clustering package
IF(USE_CLUTO)
  INCLUDE(${PAULY_CMAKE_DIR}/FindCluto.cmake)
  INCLUDE_DIRECTORIES(${CLUTO_INCLUDE_DIR})
ENDIF(USE_CLUTO)

# Set up SWIG (not the ITK version)
IF(USE_SWIG)
  INCLUDE(${PAULY_CMAKE_DIR}/FindSWIG.cmake)
ENDIF(USE_SWIG)

# Set up GLUT library
IF(USE_GLUT)
  INCLUDE(${CMAKE_ROOT}/Modules/FindGLUT.cmake)
  INCLUDE_DIRECTORIES(${GLUT_INCLUDE_DIR})
ENDIF(USE_GLUT)

# Set up LAPACK
IF(USE_LAPACK)
  INCLUDE(${PAULY_CMAKE_DIR}/FindLapack.cmake)
  IF(LAPACK_FOUND)
    SET(LAPACK_OR_MKL_FOUND 1 CACHE INTERNAL "Lapack or MKL are available")
    SET(LAPACK_OR_MKL_LIBS ${LAPACK_LIBS})
  ENDIF(LAPACK_FOUND)
ENDIF(USE_LAPACK)

# Set up INTEL MKL
IF(USE_INTEL_MKL)
  INCLUDE(${PAULY_CMAKE_DIR}/FindMKL.cmake)
  IF(MKL_FOUND)
    SET(LAPACK_OR_MKL_FOUND 1 CACHE INTERNAL "Lapack or MKL are available")
    SET(LAPACK_OR_MKL_LIBS ${MKL_LIBS})
  ENDIF(MKL_FOUND)
ENDIF(USE_INTEL_MKL)

# Set up PARDISO
IF(USE_PARDISO)
  INCLUDE(${PAULY_CMAKE_DIR}/FindPARDISO.cmake)
ENDIF(USE_PARDISO)

# Set up PARDISO
IF(USE_METIS)
  INCLUDE(${PAULY_CMAKE_DIR}/FindMETIS.cmake)
ENDIF(USE_METIS)

# Set up GSL
IF(USE_GSL)
  INCLUDE(${PAULY_CMAKE_DIR}/FindGSL.cmake)
ENDIF(USE_GSL)

# Dart Support
IF(USE_DART)
  # Include the standard Dart testing module
  INCLUDE (${CMAKE_ROOT}/Modules/Dart.cmake)
ENDIF(USE_DART)

# Set up the includes 
SET(PAULY_INCLUDE_DIRS
  ${PAULY_SOURCE_DIR}/CoreLibraries/
  ${PAULY_SOURCE_DIR}/CoreLibraries/BasicImaging
  ${PAULY_SOURCE_DIR}/CoreLibraries/BSpline
  ${PAULY_SOURCE_DIR}/CoreLibraries/FLTKViewers
  ${PAULY_SOURCE_DIR}/CoreLibraries/ITKExtensions/WindowedSincFilter/
  ${PAULY_SOURCE_DIR}/CoreLibraries/ITKExtensions/VoxBoIO/
  ${PAULY_SOURCE_DIR}/CoreLibraries/ITKExtensions/PovRayIO/
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
                                                                                
# This is necessary for scripting
IF(WIN32)
  SET(BINARY_SUBDIR "/Release")
ELSE(WIN32)
  SET(BINARY_SUBDIR "")
ENDIF(WIN32)

