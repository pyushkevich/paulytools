PROJECT(MedialPDE)

# Use the LASPACK library
INCLUDE(${PAULY_SOURCE_DIR}/CMake/FindLASPACK.cmake)
IF(LASPACK_FOUND)
  USE_LASPACK()
ENDIF(LASPACK_FOUND)

# Create the executable
ADD_EXECUTABLE(testpde MedialPDESolver.cxx MedialPDETest.cxx)
