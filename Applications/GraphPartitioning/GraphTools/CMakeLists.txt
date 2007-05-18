PROJECT(GraphTools)

INCLUDE_DIRECTORIES(
  ${GraphTools_SOURCE_DIR}
  ${PAULY_SOURCE_DIR}/Applications/BrainTracer/Code
  ${PAULY_SOURCE_DIR}/Applications/QuickITKApps/FillSTLMesh
  ${PAULY_SOURCE_DIR}/Applications/MedialPDE/Code
  ${PAULY_SOURCE_DIR}/Applications/GraphPartitioning/ImageGraphCut)

LINK_LIBRARIES(${ITK_LIBRARIES} vtkCommon vtkIO)

ADD_EXECUTABLE(imglandgraph
  BinaryImageLandmarkGraph.cxx
  ../../BrainTracer/Code/ShortestPath.cxx)

ADD_EXECUTABLE(meshlandgraph
  MeshLandmarkGraph.cxx
  ../../BrainTracer/Code/ShortestPath.cxx
  ../../QuickITKApps/FillSTLMesh/ReadWriteVTK.cxx
  ../../MedialPDE/Code/SparseMatrix.cxx
  )
