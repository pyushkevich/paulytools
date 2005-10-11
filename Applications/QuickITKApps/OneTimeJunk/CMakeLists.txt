PROJECT(ONE_TIME_JUNK)

ADD_EXECUTABLE(unc_mart2sean UNCConvertMartinImageToSeanImage.cxx)
ADD_EXECUTABLE(convert3d ConvertImage3D.cxx)
ADD_EXECUTABLE(randfx RandomEffects.cxx)

SET(COMMON_LIBS
  ITKCommon 
  ITKBasicFilters 
  ITKIO 
  ITKNumerics
  ITKStatistics)

TARGET_LINK_LIBRARIES(unc_mart2sean ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(convert3d ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(randfx ${COMMON_LIBS})
