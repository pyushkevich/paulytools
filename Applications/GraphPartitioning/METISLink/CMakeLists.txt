PROJECT(METISLINK)

# Set up the METIS link
IF(METIS_FOUND)
  # Use our custom MPREP macro to custom-compile the .tm file
  MMA_MPREP(
 	${METISLINK_SOURCE_DIR}/metis_template.tm
	${METISLINK_BINARY_DIR}/metis_template.cpp)

  ADD_EXECUTABLE(metisML MetisMMA.cpp ${METISLINK_BINARY_DIR}/metis_template.cpp)
  TARGET_LINK_LIBRARIES(metisML ${MMA_LIBRARIES} ${METIS_LIBRARIES})
ENDIF(METIS_FOUND)

# Set up the CLUTO link
IF(CLUTO_FOUND)
  # Use our custom MPREP macro to custom-compile the .tm file
  MMA_MPREP(
 	${METISLINK_SOURCE_DIR}/cluto_template.tm
	${METISLINK_BINARY_DIR}/cluto_template.cpp)

  ADD_EXECUTABLE(clutoML ClutoMMA.cpp ${METISLINK_BINARY_DIR}/cluto_template.cpp)
  TARGET_LINK_LIBRARIES(clutoML ${MMA_LIBRARIES} ${CLUTO_LIBRARIES})
ENDIF(CLUTO_FOUND)


# Copy over the MMA Notebook file that launches the link
CONFIGURE_FILE(
  ${METISLINK_SOURCE_DIR}/ClutoLink.m
  ${METISLINK_BINARY_DIR}/ClutoLink.m)
  
