# Find the directory where MKL lives
FIND_PATH(MKL_ROOT_DIR "include/mkl.h" "/opt/intel/mkl72" DOC "Location of the MKL package")

# Read the config file
IF(MKL_ROOT_DIR)

  # Read the configuration file
  LINK_DIRECTORIES(${MKL_ROOT_DIR}/lib/32)
  INCLUDE_DIRECTORIES(${MKL_ROOT_DIR}/include)

  # Set the found flag
  SET(MKL_FOUND 1 CACHE INTERNAL "LASPACK library is available")

ENDIF(MKL_ROOT_DIR)
