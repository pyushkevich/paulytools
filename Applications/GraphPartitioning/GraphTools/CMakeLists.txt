PROJECT(GraphTools)

INCLUDE_DIRECTORIES(
  ${GraphTools_SOURCE_DIR}
  ${PAULY_SOURCE_DIR}/Applications/BrainTracer/Code
  ${PAULY_SOURCE_DIR}/Applications/GraphPartitioning/ImageGraphCut)

LINK_LIBRARIES(${ITK_LIBRARIES})
ADD_EXECUTABLE(imglandgraph
  BinaryImageLandmarkGraph.cxx
  ../../BrainTracer/Code/ShortestPath.cxx)
