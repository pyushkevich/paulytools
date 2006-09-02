PROJECT(WellComposed)

INCLUDE_DIRECTORIES(${WellComposed_SOURCE_DIR})

ADD_EXECUTABLE(wellcomtest BinaryWellComposed3DImageFilter.cxx)
TARGET_LINK_LIBRARIES(wellcomtest ${ITK_LIBRARIES})

