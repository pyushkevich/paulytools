# Find the root directory for OPENMESH
FIND_PATH(OPENMESH_ROOT_DIR "Core/Mesh/BaseMesh.hh"
  DOC "Specify the root directory of the OpenMesh source tree")

# Find the OPENMESH library and includes
IF(OPENMESH_ROOT_DIR)
  
  # Compute the parent directory, since openmesh includes carry full path
  GET_FILENAME_COMPONENT(OPENMESH_PARENT_DIR ${OPENMESH_ROOT_DIR} PATH)

  # Find the include directory
  FIND_PATH(OPENMESH_INCLUDE_DIR "OpenMesh/Core/Mesh/BaseMesh.hh"
    PATHS ${OPENMESH_PARENT_DIR})

  # Find the openmesh library
  FIND_LIBRARY(OPENMESH_LIBRARY "Core" PATHS
    ${OPENMESH_ROOT_DIR}/Core 
    ${OPENMESH_ROOT_DIR}/Win/msvc7/Core/Release)
	
  # Add compiler definitions
	ADD_DEFINITIONS(-D_USE_MATH_DEFINES)

ENDIF(OPENMESH_ROOT_DIR)

# Set the OPENMESH-found flag
IF(OPENMESH_ROOT_DIR)
  IF(OPENMESH_INCLUDE_DIR)
    IF(OPENMESH_LIBRARY)
      SET (OPENMESH_FOUND 1 CACHE INTERNAL "OPENMESH library is available!")
    ENDIF(OPENMESH_LIBRARY)
  ENDIF(OPENMESH_INCLUDE_DIR)
ENDIF(OPENMESH_ROOT_DIR)
    
