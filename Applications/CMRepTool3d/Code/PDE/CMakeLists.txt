PROJECT(MedialPDE)

# Find and read the MKL package
INCLUDE(${PAULY_SOURCE_DIR}/CMake/FindMKL.cmake)

# Load Intel MKL libraries, etc
FIND_LIBRARY(PARDISO_LIBRARY NAMES pardiso DOC "Location of the PARDISO library")
FIND_LIBRARY(INTEL_FORTRAN_LIB ifcore /opt/intel_fc_80/lib DOC "Location of Intel FC library")

# Create the executable
LINK_LIBRARIES(${PARDISO_LIBRARY} mkl_lapack mkl_def ${INTEL_FORTRAN_LIB})
ADD_EXECUTABLE(testpde MedialPDESolver.cxx MedialPDETest.cxx)
