PROJECT(FillSTLMesh)

SET(FILLSTLMESH_SRCS
  FillSTLMesh.cxx
)

ADD_EXECUTABLE(mesh2img ${FILLSTLMESH_SRCS})

TARGET_LINK_LIBRARIES(mesh2img 
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

