PROJECT(FillSTLMesh)

SET(FILLSTLMESH_SRCS
  FillSTLMesh.cxx
  DrawFillInside.cpp
  DrawTriangles.cpp
)

SET(FILLSTLMESH_HEADERS
  DrawFillInside.h
  DrawTriangles.h
)

ADD_EXECUTABLE(stltoimg ${FILLSTLMESH_SRCS} ${FILLSTLMESH_HEADERS})

TARGET_LINK_LIBRARIES(stltoimg 
  ITKCommon
  ITKIO
  VTKCommon 
  VTKFiltering
  VTKImaging
  VTKIO 
  VTKGraphics 
  VTKHybrid
  VTKRendering
)

