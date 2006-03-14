PROJECT(ONE_TIME_JUNK)

ADD_EXECUTABLE(convert3d ConvertImage3D.cxx ConvertImage3DTemplates.cxx)
ADD_EXECUTABLE(identify3d IdentifyImage.cxx)

SET(COMMON_LIBS
  ITKCommon 
  ITKBasicFilters 
  ITKIO 
  ITKNumerics
  ITKStatistics
  ITKPovRayIO
  ITKVoxBoIO)

TARGET_LINK_LIBRARIES(convert3d ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(identify3d ${COMMON_LIBS})

INSTALL_TARGETS("/bin" convert3d identify3d)
