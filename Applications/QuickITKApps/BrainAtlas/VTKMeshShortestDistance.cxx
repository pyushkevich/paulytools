#include "VTKMeshShortestDistance.h"

VTKMeshShortestDistance
::VTKMeshShortestDistance(vtkPolyData *mesh)
{
  // Store the input
  m_SourcePolys = mesh;

  cout << "  input mesh has " << mesh->GetNumberOfPoints() << " points" << endl;
  cout << "  input mesh has " << mesh->GetNumberOfCells() << " cells" << endl;

  // Construct simple triangles
  fltTriangle = vtkTriangleFilter::New();
  fltTriangle->SetInput(m_SourcePolys);

  cout << "   converting mesh to triangles " << endl;
  fltTriangle->Update();

  cout << "  this mesh has " << fltTriangle->GetOutput()->GetNumberOfPoints() << " points" << endl;
  cout << "  this mesh has " << fltTriangle->GetOutput()->GetNumberOfCells() << " cells" << endl;

  // Clean the data
  fltCleaner = vtkCleanPolyData::New();
  fltCleaner->SetInput(fltTriangle->GetOutput());
  fltCleaner->SetTolerance(0);
  fltCleaner->ConvertPolysToLinesOn();
  
  cout << "   cleaning up triangle mesh " << endl;
  fltCleaner->Update();

  // Go through and delete the cells that are of the wrong type
  vtkPolyData *p = fltCleaner->GetOutput();
  for(vtkIdType i = p->GetNumberOfCells();i > 0;i--)
    {
    if(p->GetCellType(i-1) != VTK_TRIANGLE)
      p->DeleteCell(i-1);
    }
  p->BuildCells();
      
  cout << "  this mesh has " << p->GetNumberOfPoints() << " points" << endl;
  cout << "  this mesh has " << p->GetNumberOfCells() << " cells" << endl;
  
  // Compute the loop scheme
  //fltLoop = vtkLoopSubdivisionFilter::New();
  //fltLoop->SetInput(p);
  //fltLoop->SetNumberOfSubdivisions(1);

  //cout << "   subdividing the triangle mesh " << endl;
  //fltLoop->Update();
  //p = fltLoop->GetOutput();
  
  //cout << "  this mesh has " << p->GetNumberOfPoints() << " points" << endl;
  //cout << "  this mesh has " << p->GetNumberOfCells() << " cells" << endl;

  // Compute the edge map
  fltEdge = vtkFeatureEdges::New();
  fltEdge->BoundaryEdgesOff();
  fltEdge->FeatureEdgesOff();
  fltEdge->NonManifoldEdgesOff();
  fltEdge->ManifoldEdgesOn();
  fltEdge->ColoringOff();
  fltEdge->SetInput(p);

  cout << "   extracting edges from the mesh" << endl;
  fltEdge->Update();

  // Got the new poly data
  m_EdgePolys = fltEdge->GetOutput();
  m_EdgePolys->BuildCells();
  unsigned int nEdges = m_EdgePolys->GetNumberOfLines();
  cout << "      number of edges (lines) : " << nEdges << endl;
  cout << "      number of cells : " << m_EdgePolys->GetNumberOfCells() << endl;
  cout << "      number if points : " << m_EdgePolys->GetNumberOfPoints() << endl;

  // Construct a locator
  fltLocator = vtkPointLocator::New();
  fltLocator->SetDataSet(m_EdgePolys);
  fltLocator->BuildLocator();

  // Construct the cell locator
  fltCellLocator = vtkCellLocator::New();
  fltCellLocator->SetDataSet(p);
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
    m_EdgeWeights[i] = (int) (10000000 * (p1-p2).two_norm());
    if(m_EdgeWeights[i] == 0)
      cout << "Bad length at " << i << " equal to " << (p1-p2).two_norm() << endl;
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
  fltTriangle->Delete();
  fltCleaner->Delete();
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

