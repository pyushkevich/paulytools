# Find the directory where MKL lives
FIND_PATH(MKL_ROOT_DIR "include/mkl.h" "/opt/intel/mkl72" 
  DOC "Location of the Intel MKL package")

FIND_LIBRARY(INTEL_FORTRAN_LIB "ifcore" "/opt/intel_fc_80/lib" 
  DOC "Location of Intel FC (Fortran) library")

# Read the config file
IF(MKL_ROOT_DIR AND INTEL_FORTRAN_LIB)

  # Read the configuration file
  LINK_DIRECTORIES(${MKL_ROOT_DIR}/lib/32)
  INCLUDE_DIRECTORIES(${MKL_ROOT_DIR}/include)

  # Set the name of the libraries
  SET(MKL_LIBS "mkl_lapack" "mkl_def" ${INCLUDE_DIRECTORIES})

  # Set the found flag
  SET(MKL_FOUND 1 CACHE INTERNAL "MKL library is available")

ENDIF(MKL_ROOT_DIR AND INTEL_FORTRAN_LIB)
