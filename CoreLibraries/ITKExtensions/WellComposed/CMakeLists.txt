PROJECT(WellComposed)

INCLUDE_DIRECTORIES(${WellComposed_SOURCE_DIR})

ADD_LIBRARY(ITKWellComposed
  itkBinaryWellComposed3DImageFilter.h
  itkBinaryWellComposed3DImageFilter.txx)

ADD_EXECUTABLE(wellcomtest BinaryWellComposed3DImageFilter.cxx)
TARGET_LINK_LIBRARIES(wellcomtest ITKWellComposed ${ITK_LIBRARIES})

