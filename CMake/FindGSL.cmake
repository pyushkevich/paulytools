# Find the root directory for GSL
FIND_PATH(GSL_ROOT_DIR "include/gsl/gsl_const.h" "/usr/local" "/usr" 
  DOC "Specify the root directory for the GSL installation")

# Find the GSL library and includes
IF(GSL_ROOT_DIR)

  INCLUDE_DIRECTORIES(${GSL_ROOT_DIR}/include)
  
  FIND_LIBRARY(GSL_LIBRARY NAMES "gsl" PATHS ${GSL_ROOT_DIR}/lib)
  FIND_LIBRARY(GSL_CBLAS_LIBRARY NAMES "gslcblas" PATHS ${GSL_ROOT_DIR}/lib)

ENDIF(GSL_ROOT_DIR)

# Set the GSL-found flag
IF(GSL_LIBRARY)
  SET (GSL_FOUND 1 CACHE INTERNAL "GSL library is available!")
ENDIF(GSL_LIBRARY)


