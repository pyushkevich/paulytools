PROJECT(FillSTLMesh)

SET(FILLSTLMESH_SRCS
  FillSTLMesh.cxx
)

ADD_EXECUTABLE(meshdiff CompareMeshes.cxx)
ADD_EXECUTABLE(mesh2img ${FILLSTLMESH_SRCS})
ADD_EXECUTABLE(avgmesh AverageMeshes.cxx)
ADD_EXECUTABLE(identimg IdentifyImage.cxx)
ADD_EXECUTABLE(scalemesh ScaleMesh.cxx)

LINK_LIBRARIES(
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

