# Find the root directory for CLUTO
FIND_PATH(CLUTO_INCLUDE_DIR "metis.h"
  DOC "Specify the root directory for the CLUTO installation")

# Find the CLUTO library and includes
IF(CLUTO_INCLUDE_DIR)

  FIND_LIBRARY(CLUTO_LIBRARIES NAMES "cluto" "libcluto" PATHS 
    ${CLUTO_INCLUDE_DIR} ${CLUTO_INCLUDE_DIR}/Win32 ${CLUTO_INCLUDE_DIR}/Linux)

ENDIF(CLUTO_INCLUDE_DIR)

# Set the CLUTO-found flag
IF(CLUTO_INCLUDE_DIR)
  IF(CLUTO_LIBRARIES)
    SET (CLUTO_FOUND 1 CACHE INTERNAL "CLUTO library is available!")
  ENDIF(CLUTO_LIBRARIES)
ENDIF(CLUTO_INCLUDE_DIR)


