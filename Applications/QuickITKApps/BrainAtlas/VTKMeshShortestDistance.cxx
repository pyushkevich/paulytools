#include "VTKMeshShortestDistance.h"
#include "vtkIdList.h"

const double
VTKMeshShortestDistance
::EDGE_WEIGHT_FACTOR = 100000.0;

VTKMeshShortestDistance
::VTKMeshShortestDistance()
{
  // Set the distance function
  m_DefaultWeightFunction2.SetPitchFactor(1.0);
  m_WeightFunctionPtr = &m_DefaultWeightFunction;

  // Initialize the filters
  fltLocator = vtkPointLocator::New();
  fltCellLocator = vtkCellLocator::New();

  // Initialize the graph to NULL
  m_Graph = NULL;
};

void
VTKMeshShortestDistance
::SetInputMesh(vtkPolyData *mesh)
{
  // Store the input mesh
  m_SourceMesh = mesh;

  // Make sure that the links and the cells are built
  m_SourceMesh->BuildCells();
  m_SourceMesh->BuildLinks();
  
  // Construct a locator
  fltLocator->SetDataSet(m_SourceMesh);
  fltLocator->BuildLocator();

  // Construct the cell locator
  fltCellLocator->SetDataSet(m_SourceMesh);
  fltCellLocator->BuildLocator();
  
/*
  // Store the input
  m_SourceMesh = mesh;
  cout << "  input mesh has " << m_SourceMesh->GetNumberOfPoints() << " points" << endl;
  cout << "  input mesh has " << m_SourceMesh->GetNumberOfCells() << " cells" << endl;

  // Construct simple triangles
  fltTriangle->SetInput(m_SourceMesh);

  cout << "   converting mesh to triangles " << endl;
  fltTriangle->Update();

  cout << "  this mesh has " << fltTriangle->GetOutput()->GetNumberOfPoints() << " points" << endl;
  cout << "  this mesh has " << fltTriangle->GetOutput()->GetNumberOfCells() << " cells" << endl;

  // Clean the data
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
  fltLocator->SetDataSet(m_EdgePolys);
  fltLocator->BuildLocator();

  // Construct the cell locator
  fltCellLocator->SetDataSet(p);
  fltCellLocator->BuildLocator();
*/  
}

void
VTKMeshShortestDistance
::ComputeGraph()
{
  list<EdgeType> lEdges;
  list<unsigned int> lWeights;
  
  // Clear the edge arrays
  // m_Edges.clear();
  // m_EdgeWeights.clear();

  // List used to keep cells adjacent to each point
  vtkIdList *lCells = vtkIdList::New();
  vtkIdList *lAdjacent = vtkIdList::New();

  // Go trough all the points and find their neighbors
  for(vtkIdType iPoint = 0; iPoint < m_SourceMesh->GetNumberOfPoints(); iPoint++)
    {
    // Get all cells containing the point
    m_SourceMesh->GetPointCells(iPoint, lCells);

    // Reset the adjacency list
    lAdjacent->Reset();

    // Go though the cells
    for(unsigned int iCellIdx = 0; iCellIdx < lCells->GetNumberOfIds(); iCellIdx++)
      {
      // Get the cell id
      vtkIdType iCell = lCells->GetId(iCellIdx);

      // Get the list of points in the cell
      vtkIdType nCellPoints;
      vtkIdType *pCellPoints;
      m_SourceMesh->GetCellPoints(iCell, nCellPoints, pCellPoints);

      // Parse the points
      for(unsigned int iNbrIdx = 0; iNbrIdx < nCellPoints; iNbrIdx++)
        {
        // Get the neighbor ID
        vtkIdType iNbr = pCellPoints[iNbrIdx];
        
        // Determine if this is an 'actionable' edge
        if(iNbr != iPoint /* && m_SourceMesh->IsEdge(iNbr, iPoint) */ )
          {
          // Add the tartget to the adjacency list
          lAdjacent->InsertUniqueId(iNbr);
          }
        } // go through neighbor points
      } // go through cells

    // Now we have an adjacency list for this vertex. Add the edges
    for(unsigned int iAdjIdx = 0; iAdjIdx < lAdjacent->GetNumberOfIds(); iAdjIdx++)
      {
      vtkIdType iNbr = lAdjacent->GetId(iAdjIdx);

      // Add the edge to the list of edges
      if(!m_SourceMesh->IsEdge(iPoint, iNbr))
        cout << "Bad edge " << iPoint << " " << iNbr <<  endl;
      EdgeType edge(iPoint, iNbr);
      lEdges.push_back(edge);

      // Compute the associated weight
      double xWeight = 
        m_WeightFunctionPtr->GetEdgeWeight(
          m_SourceMesh, edge.first, edge.second);

      // cout << "Edge [" << edge.first << "," << edge.second << "] = " << xWeight << endl;

      // Convert the weight to integer
      lWeights.push_back((int) (EDGE_WEIGHT_FACTOR * xWeight));

      // Check that there are no zero weights (?)
      // if(m_EdgeWeights.back() <= 0)
      //   m_EdgeWeights.back() = 1;

      //if(m_EdgeWeights.back() <= 0)
      //  cout << "Edge [" << edge.first << ","  << edge.second 
      //    << "] has zero integer length (" << xWeight << ")" << endl;
      }

    } // go through points

  cout << "Number of edges: " << m_Edges.size() << endl;

  // Clean up
  lCells->Delete();
  lAdjacent->Delete();

  // Create arrays
  m_Edges.resize(lEdges.size());
  m_EdgeWeights.resize(lEdges.size());
  unsigned int iEdge = 0;
  list<EdgeType>::const_iterator itEdge = lEdges.begin();
  list<unsigned int>::const_iterator itWeight = lWeights.begin();
  while(itEdge != lEdges.end())
    {
    m_Edges[iEdge] = *itEdge;
    m_EdgeWeights[iEdge] = *itWeight;
    ++iEdge; ++itEdge; ++itWeight;
    }

  // Construct the graph from the edge information
  cout << "   constructing the graph" << endl;

  if(m_Graph)
    delete m_Graph;

  m_Graph = new GraphType(
    m_Edges.begin(), m_Edges.end(), 
    m_EdgeWeights.begin(), m_SourceMesh->GetNumberOfPoints());
}

VTKMeshShortestDistance
::~VTKMeshShortestDistance()
{
  if(m_Graph) delete m_Graph;
  fltLocator->Delete();
  fltCellLocator->Delete();
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

