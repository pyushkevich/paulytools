#include "VTKMeshShortestDistance.h"

VTKMeshShortestDistance
::VTKMeshShortestDistance(vtkPolyData *mesh)
{
  // Store the input
  m_SourcePolys = mesh;

  // Compute the edge map
  fltEdge = vtkExtractEdges::New();
  fltEdge->SetInput(m_SourcePolys);

  cout << "   extracting edges from the mesh" << endl;
  fltEdge->Update();

  // Got the new poly data
  m_EdgePolys = fltEdge->GetOutput();
  m_EdgePolys->BuildCells();
  unsigned int nEdges = m_EdgePolys->GetNumberOfLines();
  cout << "      number of edges (lines) : " << nEdges << endl;

  // Construct a locator
  fltLocator = vtkPointLocator::New();
  fltLocator->SetDataSet(m_EdgePolys);
  fltLocator->BuildLocator();

  // Construct the cell locator
  fltCellLocator = vtkCellLocator::New();
  fltCellLocator->SetDataSet(m_EdgePolys);
  fltCellLocator->BuildLocator();

  // Create an edge list
  m_Edges.resize(nEdges);
  m_EdgeWeights.resize(nEdges);

  vtkIdType nPoints = 0; vtkIdType *xPoints = NULL;
  for(unsigned int i=0;i<nEdges;i++)
    {
    // Get the next edge
    m_EdgePolys->GetCellPoints(i, nPoints, xPoints);

    // Place the edge into the Edge structure
    assert(nPoints == 2);
    m_Edges[i].first = xPoints[0];
    m_Edges[i].second = xPoints[1];

    // Compute the associated weight
    vnl_vector_fixed<double,3> p1(m_EdgePolys->GetPoint(m_Edges[i].first));
    vnl_vector_fixed<double,3> p2(m_EdgePolys->GetPoint(m_Edges[i].second));
    m_EdgeWeights[i] = (int) (1000 * (p1-p2).two_norm());
    }

  // Declare the graph object
  cout << "   constructing the graph" << endl;
  m_Graph = new GraphType(
    m_Edges.begin(), m_Edges.end(), 
    m_EdgeWeights.begin(), m_EdgePolys->GetNumberOfPoints());
}


VTKMeshShortestDistance
::~VTKMeshShortestDistance()
{
  delete m_Graph;
  fltLocator->Delete();
  fltCellLocator->Delete();
  fltEdge->Delete();
}


bool 
VTKMeshShortestDistance
::PickPoint(Vec xStart, Vec xEnd, vtkIdType &point)
{
  Vec ptLine, pCoords;
  double t; 
  vtkIdType subId;

  // Compute the intersection with the line
  int rc = fltCellLocator->IntersectWithLine(
    xStart.data_block(),xEnd.data_block(),
    0.001, t, ptLine.data_block(), pCoords.data_block(), subId);

  // cout << "Tracing from " << xStart << " to " << xEnd << endl;
  // cout << "Intersection at t = " << t << ", ptline " << ptLine << endl;
  // cout << "Vertex ID = " << subId << endl;

  // Find the vertex closest to the intersection
  point = fltLocator->FindClosestPoint(ptLine.data_block());
  return subId == 0;
}
  

void 
VTKMeshShortestDistance
::ComputeDistances(vtkIdType iStartNode)
{
  m_Distance.resize(num_vertices(*m_Graph));
  m_Predecessor.resize(num_vertices(*m_Graph));
  VertexDescriptor start = vertex(iStartNode, *m_Graph);

  dijkstra_shortest_paths(
    *m_Graph,start,predecessor_map(&m_Predecessor[0]).distance_map(&m_Distance[0]));
}

