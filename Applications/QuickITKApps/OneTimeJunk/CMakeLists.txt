PROJECT(ONE_TIME_JUNK)

ADD_EXECUTABLE(unc_mart2sean UNCConvertMartinImageToSeanImage.cxx)
ADD_EXECUTABLE(convert3d ConvertImage3D.cxx)

LINK_LIBRARIES(
  ITKCommon 
  ITKBasicFilters 
  ITKIO 
  ITKNumerics
  ITKStatistics)
