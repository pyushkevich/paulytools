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
  m_EdgeWeights = NULL;
};

void
VTKMeshShortestDistance
::SetInputMesh(VTKMeshHalfEdgeWrapper *wrapper)
{
  // Store the input mesh
  m_HalfEdge = wrapper;
  m_SourceMesh = wrapper->GetPolyData();

  // Construct a locator
  fltLocator->SetDataSet(m_SourceMesh);
  fltLocator->BuildLocator();

  // Construct the cell locator
  fltCellLocator->SetDataSet(m_SourceMesh);
  fltCellLocator->BuildLocator();

  // Set the number of vertices
  m_NumberOfVertices = m_SourceMesh->GetNumberOfPoints();
  m_NumberOfEdges = m_HalfEdge->GetNumberOfHalfEdges();
}

void
VTKMeshShortestDistance
::ComputeGraph()
{
  // Clean up the old graph info
  DeleteGraphData();

  // Allocate the graph data
  m_EdgeWeights = new float[m_NumberOfEdges];

  // Compute the length of each half-edge in the graph
  for(unsigned int i=0;i<m_NumberOfEdges;i++)
    {
    vtkIdType iHead = m_HalfEdge->GetHalfEdgeVertex(i);
    vtkIdType iTail = m_HalfEdge->GetHalfEdgeTailVertex(i);

    // Compute the edge weight
    m_EdgeWeights[i] = 
      (float) m_WeightFunctionPtr->GetEdgeWeight(m_SourceMesh, iHead, iTail);
    }

  // Create the shortest path object
  m_ShortestPath = new DijkstraAlgorithm( 
    m_NumberOfVertices, 
    m_HalfEdge->GetAdjacencyIndex(), 
    m_HalfEdge->GetAdjacency(), 
    m_EdgeWeights);
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
::PickCell(Vec xStart, Vec xEnd, vtkIdType &point) const
{
  Vec ptLine, pCoords;
  vtkIdType subId;
  double t; 

  // Compute the intersection with the line
  int rc = fltCellLocator->IntersectWithLine(
    xStart.data_block(),xEnd.data_block(),
    0.001, t, ptLine.data_block(), pCoords.data_block(), subId, point);

  return subId == 0;
}

bool 
VTKMeshShortestDistance
::PickPoint(Vec xStart, Vec xEnd, vtkIdType &point, ICellChecher *cbCell) const
{
  Vec ptLine, pCoords;
  double t; 
  int subId; vtkIdType iCell;

  do
    {
    // cout << "Searching ray " << xStart << "  to  " << xEnd << endl;
    
    // Compute the intersection with the line
    int rc = fltCellLocator->IntersectWithLine(
      xStart.data_block(),xEnd.data_block(),
      0.001, t, ptLine.data_block(), pCoords.data_block(), subId, iCell);

    // cout << "   RC: " << rc << " \t T = " << t << " \t iCell " << iCell << endl;

    // If no intersection found, return false
    if(!rc) return false;

    // Increase the starting vector by t + epsilon
    xStart += (t + 0.001) * (xEnd - xStart);
    }
  while(cbCell && !cbCell->CheckCell(iCell));

  // cout << "Tracing from " << xStart << " to " << xEnd << endl;
  // cout << "Intersection at t = " << t << ", ptline " << ptLine << endl;
  // cout << "Vertex ID = " << subId << endl;

  // Find the vertex closest to the intersection
  point = fltLocator->FindClosestPoint(ptLine.data_block());
  return subId == 0;
}

void 
VTKMeshShortestDistance
::ComputeDistances(const list<vtkIdType> &lSources)
{
  // Create an array from the list
  unsigned int *lSourceArray = new unsigned int[lSources.size()];
  unsigned int iSource = 0;

  list<vtkIdType>::const_iterator it = lSources.begin();
  while(it != lSources.end())
    lSourceArray[iSource++] = *it++;

  m_ShortestPath->ComputePathsFromManySources(iSource, lSourceArray);
  
  delete[] lSourceArray;
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
    delete m_EdgeWeights;
    }
}
