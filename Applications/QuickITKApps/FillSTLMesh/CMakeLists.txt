PROJECT(FillSTLMesh)

SET(FILLSTLMESH_SRCS
  FillSTLMesh.cxx
)

ADD_EXECUTABLE(stltoimg ${FILLSTLMESH_SRCS})

TARGET_LINK_LIBRARIES(stltoimg 
  BasicImaging
  ITKCommon
  ITKIO
  vtkCommon 
  vtkFiltering
  vtkImaging
  vtkIO 
  vtkGraphics 
  vtkHybrid
  vtkRendering
)

