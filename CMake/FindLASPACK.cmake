# Find the directory where LASPACK config file lives
FIND_PATH(LASPACK_ROOT_DIR "LASPACK_Config.cmake" DOC "Specify the directiry where LasPack was built")

# Read the config file
IF(LASPACK_ROOT_DIR)
  
  # Read the configuration file
  INCLUDE(${LASPACK_ROOT_DIR}/LASPACK_Config.cmake)

  # Set the found flag
  SET (LASPACK_FOUND 1 CACHE INTERNAL "LASPACK library is available")

ENDIF(LASPACK_ROOT_DIR)
