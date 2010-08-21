PROJECT(FillSTLMesh)

SET(FILLSTLMESH_SRCS
  FillSTLMesh.cxx
)

INCLUDE_DIRECTORIES(${PAULY_SOURCE_DIR}/Applications/MedialPDE/Code)

ADD_EXECUTABLE(tetdumpattr TetDumpAttr.cxx)
ADD_EXECUTABLE(vtklevelset RealImageToMesh.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(meshcluster SurfaceClusterAnalysis.cxx ReadWriteVTK.cxx
  ${PAULY_SOURCE_DIR}/Applications/MedialPDE/Code/Registry.cxx)
ADD_EXECUTABLE(meshdiff CompareMeshes.cxx)
ADD_EXECUTABLE(mesh2img ${FILLSTLMESH_SRCS})
ADD_EXECUTABLE(avgmesh AverageMeshes.cxx)
ADD_EXECUTABLE(submesh SubdivideMesh.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(quadcluster QuadClustering.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(smoothmesh SmoothMesh.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(copymesh CopyMesh.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(identimg IdentifyImage.cxx)
ADD_EXECUTABLE(scalemesh ScaleMesh.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(dumpmeshpoints DumpMeshPoints.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(importqvoronoi ImportQVoronoi.cxx ReadWriteVTK.cxx
  ${PAULY_SOURCE_DIR}/Applications/BrainTracer/Code/VTKMeshShortestDistance.cxx)
ADD_EXECUTABLE(dumpmeshattr DumpMeshAttr.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(warpmesh WarpMeshBackwards.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(warpmesh2 WarpMeshBackwards2.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(flipmesh FlipMesh.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(cmrep2img CMRepMeshToImage.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(cmrep2warp CMRepToWarp.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(mesh_image_sample MeshImageSample.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(avgmesharr AverageMeshArrays.cxx)
ADD_EXECUTABLE(mesharrstat  MeshArrayStats.cxx ReadWriteVTK.cxx)
# ADD_EXECUTABLE(mesh_array_average  MeshArrayAverage.cxx ReadWriteVTK.cxx)
# ADD_EXECUTABLE(vtk2off  ConvertVTKtoOFF.cxx ReadWriteVTK.cxx)
ADD_EXECUTABLE(surf2tet SurfMesh2VolMesh.cxx ReadWriteVTK.cxx TetraLib.cxx)
ADD_EXECUTABLE(contour2surf ContourMesh2SurfMesh.cxx ReadWriteVTK.cxx) 
ADD_EXECUTABLE(tetgen2vtk TetGenToVTK.cxx ReadWriteVTK.cxx TetraLib.cxx)
ADD_EXECUTABLE(tetfill TetraFill.cxx)
ADD_EXECUTABLE(tetsample TetSampleImage.cxx)
ADD_EXECUTABLE(tetjac TetJacobian.cxx)
ADD_EXECUTABLE(tetmergearr MergeVTKArrays.cxx)

### ADD_EXECUTABLE(importqvoronoi 
###  ImportQVoronoi.cxx
###  ReadWriteVTK.cxx
###  ${PAULY_SOURCE_DIR}/Applications/BrainTracer/Code/VTKMeshShortestDistance.cxx
###)

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

INCLUDE_DIRECTORIES(${PAULY_SOURCE_DIR}/ExternalLibraries/tetgen1.4.2)
INCLUDE_DIRECTORIES(${PAULY_SOURCE_DIR}/ExternalLibraries/triangle)
INCLUDE_DIRECTORIES(${PAULY_SOURCE_DIR}/CoreLibraries/ITKExtensions)

TARGET_LINK_LIBRARIES(tetfill ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(meshdiff ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(avgmesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(submesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(quadcluster ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(smoothmesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(copymesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(avgmesharr ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(mesharrstat ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(identimg ${ITK_LIBS})
TARGET_LINK_LIBRARIES(scalemesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(dumpmeshpoints ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(meshcluster ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(dumpmeshattr ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(importqvoronoi ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(mesh2img ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(flipmesh ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(vtklevelset ${COMMON_LIBS} ${ITK_LIBS})
TARGET_LINK_LIBRARIES(warpmesh ${COMMON_LIBS} ${ITK_LIBS})
TARGET_LINK_LIBRARIES(warpmesh2 ${COMMON_LIBS} ${ITK_LIBS})
TARGET_LINK_LIBRARIES(cmrep2img ${COMMON_LIBS} ${ITK_LIBS})
TARGET_LINK_LIBRARIES(cmrep2warp ${COMMON_LIBS} ${ITK_LIBS})
TARGET_LINK_LIBRARIES(mesh_image_sample ${COMMON_LIBS} ${ITK_LIBS})
# TARGET_LINK_LIBRARIES(mesh_array_average ${COMMON_LIBS} ${ITK_LIBS})
# TARGET_LINK_LIBRARIES(vtk2off ${COMMON_LIBS} ${ITK_LIBS})
TARGET_LINK_LIBRARIES(contour2surf ${COMMON_LIBS} triangle)
TARGET_LINK_LIBRARIES(surf2tet ${COMMON_LIBS} tetgen)
TARGET_LINK_LIBRARIES(tetgen2vtk ${COMMON_LIBS} tetgen)
TARGET_LINK_LIBRARIES(tetsample ${COMMON_LIBS} ${ITK_LIBS})
TARGET_LINK_LIBRARIES(tetjac ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(tetdumpattr ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(tetmergearr ${COMMON_LIBS})
# SET_TARGET_PROPERTIES(contour2surf PROPERTIES LINK_FLAGS -ltriangle)

INSTALL_TARGETS("/bin" meshdiff avgmesh submesh quadcluster 
  smoothmesh copymesh avgmesharr identimg scalemesh 
  dumpmeshpoints meshcluster dumpmeshattr 
  mesh2img flipmesh vtklevelset warpmesh warpmesh2 cmrep2img cmrep2warp
  mesh_image_sample contour2surf surf2tet tetgen2vtk tetfill tetjac tetdumpattr tetsample
tetmergearr)
