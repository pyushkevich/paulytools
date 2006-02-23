PROJECT(ONE_TIME_JUNK)

ADD_EXECUTABLE(unc_mart2sean UNCConvertMartinImageToSeanImage.cxx)
ADD_EXECUTABLE(convert3d ConvertImage3D.cxx)
ADD_EXECUTABLE(blur3d BlurImage3D.cxx)
ADD_EXECUTABLE(randfx RandomEffects.cxx)
ADD_EXECUTABLE(tempfile TempFile.c)
ADD_EXECUTABLE(maskmunge MaskMunge.cxx)
ADD_EXECUTABLE(findpeaks FindPeaks.cxx)

SET(COMMON_LIBS
  ITKCommon 
  ITKBasicFilters 
  ITKIO 
  ITKNumerics
  ITKStatistics
  ITKVoxBoIO)

TARGET_LINK_LIBRARIES(unc_mart2sean ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(convert3d ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(blur3d ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(randfx ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(maskmunge ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(findpeaks ${COMMON_LIBS})

INSTALL_TARGETS("/bin" convert3d blur3d randfx tempfile maskmunge findpeaks)
