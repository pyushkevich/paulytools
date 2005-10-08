PROJECT(STACKEM)

ADD_EXECUTABLE(stackem Stackem.cxx)

TARGET_LINK_LIBRARIES(stackem
  ITKCommon
  ITKIO
  ITKNumerics
)
