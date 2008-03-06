PROJECT(FillSTLMesh)

SET(FILLSTLMESH_SRCS
  FillSTLMesh.cxx
)

INCLUDE_DIRECTORIES(${PAULY_SOURCE_DIR}/Applications/MedialPDE/Code)

ADD_EXECUTABLE(meshcluster SurfaceClusterAnalysis.cxx ReadWriteVTK.cxx
  ${PAULY_SOURCE_DIR}/Applications/MedialPDE/Code/Registry.cxx)
ADD_EXECUTABLE(meshdiff CompareMeshes.cxx)
ADD_EXECUTABLE(mesh2img ${FILLSTLMESH_SRCS})
ADD_EXECUTABLE(avgmesh AverageMeshes.cxx)
ADD_EXECUTABLE(submesh SubdivideMesh.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(quadcluster QuadClustering.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(smoothmesh SmoothMesh.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(copymesh CopyMesh.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(avgmesharr AverageMeshArrays.cxx)
ADD_EXECUTABLE(identimg IdentifyImage.cxx)
ADD_EXECUTABLE(scalemesh ScaleMesh.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(dumpmeshpoints DumpMeshPoints.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(dumpmeshattr DumpMeshAttr.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(real2mesh RealImageToMesh.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(WarpMeshBackward WarpMeshBackwards.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(flipmesh FlipMesh.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(cmrep2img CMRepMeshToImage.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(cmrep2warp CMRepToWarp.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(mesh_image_sample MeshImageSample.cxx ReadWriteVTK.cxx)

ADD_EXECUTABLE(importqvoronoi 
  ImportQVoronoi.cxx
  ReadWriteVTK.cxx
  ${PAULY_SOURCE_DIR}/Applications/BrainTracer/Code/VTKMeshShortestDistance.cxx
)

SET(ITK_LIBS
  BasicImaging
  ITKCommon
  ITKIO
  ITKVoxBoIO
  ITKBasicFilters
  ITKAlgorithms
)

SET(COMMON_LIBS
  ${ITK_LIBS}
  ${VTK_COMMON_LIBS}
)

INCLUDE_DIRECTORIES(${MYLIBS_INCLUDE_DIRS}
  ${PAULY_SOURCE_DIR}/Applications/BrainTracer/Code)

TARGET_LINK_LIBRARIES(meshdiff ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(avgmesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(submesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(quadcluster ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(smoothmesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(copymesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(avgmesharr ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(identimg ${ITK_LIBS})
TARGET_LINK_LIBRARIES(scalemesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(dumpmeshpoints ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(meshcluster ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(dumpmeshattr ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(importqvoronoi ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(mesh2img ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(flipmesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(real2mesh ${COMMON_LIBS} ${ITK_LIBS})
TARGET_LINK_LIBRARIES(WarpMeshBackward ${COMMON_LIBS} ${ITK_LIBS})
TARGET_LINK_LIBRARIES(cmrep2img ${COMMON_LIBS} ${ITK_LIBS})
TARGET_LINK_LIBRARIES(cmrep2warp ${COMMON_LIBS} ${ITK_LIBS})
TARGET_LINK_LIBRARIES(mesh_image_sample ${COMMON_LIBS} ${ITK_LIBS})

INSTALL_TARGETS("/bin" meshdiff avgmesh submesh quadcluster 
  smoothmesh copymesh avgmesharr identimg scalemesh 
  dumpmeshpoints meshcluster dumpmeshattr importqvoronoi 
  mesh2img flipmesh real2mesh WarpMeshBackward cmrep2img cmrep2warp
  mesh_image_sample)
