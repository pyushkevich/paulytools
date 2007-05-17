PROJECT(ONE_TIME_JUNK)

ADD_LIBRARY(cndtemplates
  ConvertImage2DTemplates.cxx
  ConvertImage3DTemplates.cxx)

ADD_EXECUTABLE(c3d 
  ConvertImage3D.cxx 
  Convert3DMain.cxx)

ADD_EXECUTABLE(c2d 
  ConvertImage3D.cxx 
  Convert2DMain.cxx)

SET(COMMON_LIBS
  ITKCommon 
  ITKNumerics
  ITKBasicFilters 
  ITKIO 
  ITKPovRayIO
  ITKVoxBoIO)

TARGET_LINK_LIBRARIES(c3d cndtemplates ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(c2d cndtemplates ${COMMON_LIBS})

INSTALL_TARGETS("/bin" c3d c2d)
