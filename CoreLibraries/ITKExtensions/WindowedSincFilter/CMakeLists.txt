PROJECT(TEST_SINC)

ADD_EXECUTABLE("sincfilter2d"
  itkTestWindowedSincInterpolateImageFunction.cxx)

# Set up ITK
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
    "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

# Link Libraries
LINK_LIBRARIES(
  ITKCommon 
  ITKBasicFilters 
  ITKIO 
  ITKNumerics
)
