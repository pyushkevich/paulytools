#include "VTKMeshShortestDistance.h"
#include "vtkIdList.h"

VTKMeshShortestDistance
::VTKMeshShortestDistance()
{
  // Set the distance function
  m_WeightFunctionPtr = &m_DefaultWeightFunction;

  // Initialize the filters
  fltLocator = vtkPointLocator::New();
  fltCellLocator = vtkCellLocator::New();

  // Initialize the graph to NULL
  m_ShortestPath = NULL;
  m_Adjacency = m_AdjacencyIndex = NULL;
  m_EdgeWeights = NULL;
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

  // Set the number of vertices
  m_NumberOfVertices = m_SourceMesh->GetNumberOfPoints();
}

void
VTKMeshShortestDistance
::ComputeGraph()
{
  unsigned int iPoint, iCellIdx, iNbrIdx;  

  // Keep track of number of vertices, edges
  m_NumberOfEdges = 0;

  // Create an array of adjacency lists (is there a better way to do this?)
  vector<vtkIdList *> lAdjacency(m_NumberOfVertices, NULL);

  // List used to keep cells adjacent to each point
  vtkIdList *lCells = vtkIdList::New();
  vtkIdList *lNeighbors = vtkIdList::New();

  // Go trough all the points and find their neighbors
  for(iPoint = 0; iPoint < m_NumberOfVertices; iPoint++)
    {
    // Initialize the adjacency list
    lAdjacency[iPoint] = vtkIdList::New();

    // Get all cells containing the point
    m_SourceMesh->GetPointCells(iPoint, lCells);

    // Go though the cells
    for(iCellIdx = 0; iCellIdx < lCells->GetNumberOfIds(); iCellIdx++)
      {
      // Get the cell id
      vtkIdType iCell = lCells->GetId(iCellIdx);

      // Get the list of points in the cell
      m_SourceMesh->GetCellPoints(iCell, lNeighbors);

      // Parse the points
      for(iNbrIdx = 0; iNbrIdx < lNeighbors->GetNumberOfIds(); iNbrIdx++)
        {
        // Get the neighbor ID
        vtkIdType iNbr = lNeighbors->GetId(iNbrIdx);
        
        // Determine if this is an 'actionable' edge
        if(iNbr != iPoint && m_SourceMesh->IsEdge(iNbr, iPoint))
          {
          // Add the tartget to the adjacency list
          lAdjacency[iPoint]->InsertUniqueId(iNbr);

          // Increment the edge counter
          m_NumberOfEdges++;
          }
        } 
      }
    }

  // We have traversed all the points and constructed an adjacency matrix. 
  // Now we must compute the edge weights associated with the edges, 
  // and place the graph info into a compact structure.

  // Clean up the old graph info
  DeleteGraphData();

  // Allocate the graph data
  m_AdjacencyIndex = new unsigned int[m_NumberOfVertices + 1];
  m_Adjacency = new unsigned int[m_NumberOfEdges];
  m_EdgeWeights = new float[m_NumberOfEdges];

  // Populate the edge data
  unsigned int iLastEdge = 0;
  for(iPoint = 0; iPoint < m_NumberOfVertices; iPoint++)
    {
    // Point the adjacency index to the right place
    m_AdjacencyIndex[iPoint] = iLastEdge;

    // Fill the adjacencies
    for(iNbrIdx = 0; iNbrIdx < lAdjacency[iPoint]->GetNumberOfIds(); iNbrIdx++)
      {
      // Get the adjacent vertex
      vtkIdType iAdj = lAdjacency[iPoint]->GetId(iNbrIdx);

      // Compute the edge weight
      m_EdgeWeights[iLastEdge] = 
        (float) m_WeightFunctionPtr->GetEdgeWeight(m_SourceMesh, iPoint, iAdj);

      // Store the edge itself
      m_Adjacency[iLastEdge] = iAdj;

      // Increment the edge counter
      iLastEdge++;
      }

    // Delete the adjacency list
    lAdjacency[iPoint]->Delete();
    }

  // Set the last adjacency
  m_AdjacencyIndex[iPoint] = iLastEdge;

  cout << "Number of edges: " << m_NumberOfEdges << endl;

  // Clean up temporary data 
  lCells->Delete();
  lNeighbors->Delete();

  // Create the shortest path object
  m_ShortestPath = new DijkstraAlgorithm(
    m_NumberOfVertices, m_AdjacencyIndex, m_Adjacency, m_EdgeWeights);
}

VTKMeshShortestDistance
::~VTKMeshShortestDistance()
{
  DeleteGraphData();
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
  m_ShortestPath->ComputePathsFromSource(iStartNode);
}

void 
VTKMeshShortestDistance
::DeleteGraphData()
{
  if(m_ShortestPath)
    {
    delete m_ShortestPath; m_ShortestPath = NULL;
    delete m_Adjacency;
    delete m_AdjacencyIndex;
    delete m_EdgeWeights;
    }
}
