# Find the root directory for METIS
FIND_PATH(METIS_ROOT_DIR "Lib/metis.h"
  DOC "Specify the root directory for the METIS installation")

# Find the METIS library and includes
IF(METIS_ROOT_DIR)

  FIND_LIBRARY(METIS_LIBRARIES "metis" PATHS 
    ${METIS_ROOT_DIR} ${METIS_ROOT_DIR}/Release ${METIS_ROOT_DIR}/Debug)

  FIND_PATH(METIS_INCLUDE_DIR "metis.h" ${METIS_ROOT_DIR}/Lib)

ENDIF(METIS_ROOT_DIR)

# Set the METIS-found flag
IF(METIS_ROOT_DIR)
  IF(METIS_INCLUDE_DIR)
    IF(METIS_LIBRARIES)
      SET (METIS_FOUND 1 CACHE INTERNAL "METIS library is available!")
    ENDIF(METIS_LIBRARIES)
  ENDIF(METIS_INCLUDE_DIR)
ENDIF(METIS_ROOT_DIR)

