PROJECT(AntiAlias)

OPTION(USE_BOOST "Use BOOST graph library" ON)

ADD_EXECUTABLE(aalias AntiAlias.cpp)
ADD_EXECUTABLE(lsplit LabelledToBinary.cpp)
ADD_EXECUTABLE(lcombine BinaryToLabelled.cpp)
ADD_EXECUTABLE(labelribbons LabelRibbons.cpp)

SET(COMMON_LIBS
  ITKCommon 
  ITKBasicFilters 
  ITKIO 
  ${PAULY_CORE_LIBS}
  vtkCommon 
  vtkFiltering
  vtkImaging
  vtkIO 
  vtkGraphics 
  vtkHybrid
  vtkRendering
  vtkPatented)

TARGET_LINK_LIBRARIES(aalias ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(lsplit ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(lcombine ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(labelribbons ${COMMON_LIBS})

IF(USE_BOOST)
  ADD_EXECUTABLE(bpaint BrainPaint.cpp)
ENDIF(USE_BOOST)
