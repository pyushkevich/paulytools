PROJECT(AntiAlias)

ADD_EXECUTABLE(aalias AntiAlias.cpp)
ADD_EXECUTABLE(bpaint BrainPaint.cpp)
ADD_EXECUTABLE(lsplit LabelledToBinary.cpp)
LINK_LIBRARIES(ITKCommon ITKBasicFilters ITKIO ${PAULY_CORE_LIBS}
  vtkCommon 
  vtkFiltering
  vtkImaging
  vtkIO 
  vtkGraphics 
  vtkHybrid
  vtkRendering
)

