PROJECT(LASPACK)

# Option to build the LASPACK example programs
OPTION(BUILD_LASPACK_EXAMPLES "Build example programs supplied with Laspack" OFF)

INCLUDE_DIRECTORIES(${LASPACK_SOURCE_DIR})

SET(LASP_SOURCES
  laspack/eigenval.c
  laspack/errhandl.c
  laspack/factor.c
  laspack/itersolv.c
  laspack/matrix.c
  laspack/mlsolv.c
  laspack/operats.c
  laspack/precond.c
  laspack/qmatrix.c
  laspack/rtc.c
  laspack/vector.c)

SET(LASP_HEADERS
  laspack/copyrght.h
  laspack/eigenval.h
  laspack/elcmp.h
  laspack/errhandl.h
  laspack/factor.h
  laspack/itersolv.h
  laspack/lastypes.h
  laspack/matrix.h
  laspack/mlsolv.h
  laspack/operats.h
  laspack/precond.h
  laspack/qmatrix.h
  laspack/rtc.h
  laspack/vector.h
  laspack/version.h)

# Create a laspack based library 
ADD_LIBRARY(laspack ${LASP_SOURCES} ${LASP_HEADERS})

# Build the example programs if asked to
IF(BUILD_LASPACK_EXAMPLES)
  
  # Build the XC library
  ADD_LIBRARY(xc xc/getopts.c xc/xstring.c) 

  # Link the new libs
  LINK_LIBRARIES(laspack xc)

  # Build the lastest program
  ADD_EXECUTABLE(lastest ${LASPACK_SOURCE_DIR}/laspack/examples/lastest/lastest.c)

  # Build the matropt prog
  SET(MATROPT_DIR ${LASPACK_SOURCE_DIR}/laspack/examples/matropt)
  INCLUDE_DIRECTORIES(${MATROPT_DIR})
  ADD_EXECUTABLE(matropt ${MATROPT_DIR}/matropt.c ${MATROPT_DIR}/testproc.c)

  # Build the vectopt prog
  SET(VECTOPT_DIR ${LASPACK_SOURCE_DIR}/laspack/examples/vectopt)
  INCLUDE_DIRECTORIES(${VECTOPT_DIR})
  ADD_EXECUTABLE(vectopt ${VECTOPT_DIR}/vectopt.c ${VECTOPT_DIR}/testproc.c)

  # Build the MLSTest prog
  SET(MLSTEST_DIR ${LASPACK_SOURCE_DIR}/laspack/examples/mlstest)
  INCLUDE_DIRECTORIES(${MLSTEST_DIR})
  ADD_EXECUTABLE(mlstest ${MLSTEST_DIR}/extrsolv.c ${MLSTEST_DIR}/mlstest.c)

ENDIF(BUILD_LASPACK_EXAMPLES)

# Create a CMake configuration file for using LasPack in other programs
# Set the include directory
SET(LASPACK_INCLUDE_DIR ${LASPACK_SOURCE_DIR})

# Set the library directory
SET(LASPACK_LIBRARY_DIR ${LASPACK_BINARY_DIR})

# Pass on the library and include paths to config file in binary directory
CONFIGURE_FILE(
	"${LASPACK_SOURCE_DIR}/LASPACK_Config.cmake.in"
	"${LASPACK_BINARY_DIR}/LASPACK_Config.cmake"
	IMMEDIATE @ONLY
)


