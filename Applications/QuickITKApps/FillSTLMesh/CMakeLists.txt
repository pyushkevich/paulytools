PROJECT(FillSTLMesh)

SET(FILLSTLMESH_SRCS
  FillSTLMesh.cxx
)

ADD_EXECUTABLE(meshdiff CompareMeshes.cxx)
ADD_EXECUTABLE(mesh2img ${FILLSTLMESH_SRCS})
ADD_EXECUTABLE(avgmesh AverageMeshes.cxx)
ADD_EXECUTABLE(identimg IdentifyImage.cxx)
ADD_EXECUTABLE(scalemesh ScaleMesh.cxx)

SET(ITK_LIBS
  BasicImaging
  ITKCommon
  ITKIO
)

SET(COMMON_LIBS
  ${ITK_LIBS}
  vtkCommon 
  vtkFiltering
  vtkImaging
  vtkIO 
  vtkGraphics 
  vtkHybrid
  vtkRendering
)

TARGET_LINK_LIBRARIES(meshdiff ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(avgmesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(identimg ${ITK_LIBS})
TARGET_LINK_LIBRARIES(scalemesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(mesh2img ${COMMON_LIBS})
