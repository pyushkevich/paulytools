# Optionally include the Mathematica directory
OPTION(BUILD_METIS_MMA "Should I build the Mathematica interface to METIS?" OFF)

VARIABLE_REQUIRES(BUILD_METIS_MMA DO_BUILD_METIS_MMA
  MMA_FOUND
  METIS_FOUND
  CLUTO_FOUND
)

IF(BUILD_METIS_MMA)
  SUBDIRS(METISLink)  
ENDIF(BUILD_METIS_MMA)

# Optionally include the Graph Cut directory
OPTION(BUILD_METIS_GCUT "Should I build the image/graph cut application based on METIS?" ON)
VARIABLE_REQUIRES(BUILD_METIS_GCUT DO_BUILD_METIS_GCUT METIS_FOUND)

IF(BUILD_METIS_GCUT)
  SUBDIRS(ImageGraphCut)  
ENDIF(BUILD_METIS_GCUT)
