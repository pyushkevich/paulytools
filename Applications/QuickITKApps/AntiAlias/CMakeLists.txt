PROJECT(AntiAlias)

OPTION(USE_BOOST "Use BOOST graph library" ON)

ADD_EXECUTABLE(aalias AntiAlias.cpp)
ADD_EXECUTABLE(lsplit LabelledToBinary.cpp)
LINK_LIBRARIES(ITKCommon ITKBasicFilters ITKIO ${PAULY_CORE_LIBS}
  vtkCommon 
  vtkFiltering
  vtkImaging
  vtkIO 
  vtkGraphics 
  vtkHybrid
  vtkRendering
  vtkPatented
)

IF(USE_BOOST)
  ADD_EXECUTABLE(bpaint BrainPaint.cpp)
  INCLUDE_DIRECTORIES("/usr/local/include/boost-1_31")
ENDIF(USE_BOOST)
