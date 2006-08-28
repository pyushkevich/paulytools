PROJECT(FillSTLMesh)

SET(FILLSTLMESH_SRCS
  FillSTLMesh.cxx
)

ADD_EXECUTABLE(meshdiff CompareMeshes.cxx)
ADD_EXECUTABLE(mesh2img ${FILLSTLMESH_SRCS})
ADD_EXECUTABLE(avgmesh AverageMeshes.cxx)
ADD_EXECUTABLE(avgmesharr AverageMeshArrays.cxx)
ADD_EXECUTABLE(identimg IdentifyImage.cxx)
ADD_EXECUTABLE(scalemesh ScaleMesh.cxx)
ADD_EXECUTABLE(dumpmeshpoints DumpMeshPoints.cxx)

ADD_EXECUTABLE(importqvoronoi 
  ImportQVoronoi.cxx
  ${PAULY_SOURCE_DIR}/Applications/BrainTracer/Code/VTKMeshShortestDistance.cxx
)

SET(ITK_LIBS
  BasicImaging
  ITKCommon
  ITKIO
  ITKVoxBoIO
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

INCLUDE_DIRECTORIES(${MYLIBS_INCLUDE_DIRS}
  ${PAULY_SOURCE_DIR}/Applications/BrainTracer/Code)

TARGET_LINK_LIBRARIES(meshdiff ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(avgmesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(avgmesharr ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(identimg ${ITK_LIBS})
TARGET_LINK_LIBRARIES(scalemesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(dumpmeshpoints ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(importqvoronoi ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(mesh2img ${COMMON_LIBS})

INSTALL_TARGETS("/bin" mesh2img meshdiff avgmesh identimg scalemesh
  dumpmeshpoints)
