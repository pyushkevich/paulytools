PROJECT(TESTSINC)

ADD_EXECUTABLE("sincfilter2d"
  itkTestWindowedSincInterpolateImageFunction.cxx)

INCLUDE_DIRECTORIES(${TESTSINC_SOURCE_DIR})

# Set up ITK
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
    "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

# Link Libraries
TARGET_LINK_LIBRARIES(sincfilter2d
  ITKCommon 
  ITKBasicFilters 
  ITKIO 
  ITKNumerics
)
