#include "VTKMeshVoronoiDiagram.h"
#include "vtkIdList.h"
#include <vector>
#include <map>
#include <set>

using namespace std;

VTKMeshVoronoiDiagram
::VTKMeshVoronoiDiagram()
{
  // Initialize the graph to NULL
  m_Diagram = NULL;
  m_Adjacency = m_AdjacencyIndex = NULL;
  m_DualEdgeMap = NULL;
  m_RawEdgeWeights = m_FullEdgeWeights = NULL;
  m_IsDiagramValid = false;
}

VTKMeshVoronoiDiagram
::~VTKMeshVoronoiDiagram()
{
  DeleteGraphData();
}

void
VTKMeshVoronoiDiagram
::SetInputMesh(VTKMeshHalfEdgeWrapper *halfedge)
{
  // Store the input mesh
  m_HalfEdge = halfedge;
  m_Mesh = m_HalfEdge->GetPolyData();

  // Set the number of vertices to the number of cells in the mesh
  m_NumberOfVertices = m_Mesh->GetNumberOfCells();
  
  // There is no valid diagram now
  m_IsDiagramValid = false;
}

void
VTKMeshVoronoiDiagram
::ComputeGraph()
{
  // Clean up the old graph info
  DeleteGraphData();

  // Allocate the graph adjacency index info
  m_AdjacencyIndex = new unsigned int[m_NumberOfVertices + 1];
  m_AdjacencyIndex[0] = 0;
  
  // Allocate the edge arrays
  m_Adjacency = new unsigned int[m_HalfEdge->GetNumberOfHalfEdges()];
  m_DualEdgeMap = new unsigned int[m_HalfEdge->GetNumberOfHalfEdges()];

  // Create an array to hold cell center points
  Vec *xCenter = new Vec[m_NumberOfVertices];

  // Iterate over all cells in the mesh
  unsigned int iDualEdge = 0;
  for(unsigned int iCell = 0; iCell < m_Mesh->GetNumberOfCells(); iCell++)
    {
    // Initialize the center
    xCenter[iCell].fill(0.0);

    // Find all faces adjacent to this one
    unsigned int iEdge = m_HalfEdge->GetFaceHalfEdge(iCell);
    unsigned int iTest = iEdge;
    do 
      {
      // Find the adjacent face to this one
      m_Adjacency[iDualEdge] = m_HalfEdge->GetHalfEdgeFace(
        m_HalfEdge->GetHalfEdgeOpposite(iTest));

      // Record the dual edge
      m_DualEdgeMap[iTest] = iDualEdge;

      // Add point coordinate to center
      double *xPoint = m_Mesh->GetPoint(m_HalfEdge->GetHalfEdgeVertex(iTest));
      xCenter[iCell][0] += xPoint[0];
      xCenter[iCell][1] += xPoint[1];
      xCenter[iCell][2] += xPoint[2];

      // Go to the next edge
      iTest = m_HalfEdge->GetHalfEdgeNext(iTest);

      // Increment the adjacency index
      iDualEdge++;
      } 
    while(iTest != iEdge);

    // Adjust the center value
    xCenter[iCell] *= 1.0 / (iDualEdge - m_AdjacencyIndex[iCell]);

    // Set the adjacency index for the next cell
    m_AdjacencyIndex[iCell+1] = iDualEdge;
    }
  
  // Set the number of edges
  m_NumberOfEdges = m_AdjacencyIndex[m_NumberOfVertices];

  // Compute the edge weights 
  m_RawEdgeWeights = new float[m_NumberOfEdges];
  m_FullEdgeWeights = new float[m_NumberOfEdges];

  for(unsigned int iEdge = 0; iEdge < m_NumberOfEdges; iEdge++)
    {
    // Find the faces at the two sides of the edge
    unsigned int iDual = m_DualEdgeMap[iEdge];
    vtkIdType iFaceLeft = m_HalfEdge->GetHalfEdgeFace(iDual);
    vtkIdType iFaceRight = m_HalfEdge->GetHalfEdgeFace(
      m_HalfEdge->GetHalfEdgeOpposite(iDual));

    // Compute the distance between centers
    m_FullEdgeWeights[iEdge] = m_RawEdgeWeights[iEdge] = 
      (float) (xCenter[iFaceLeft] - xCenter[iFaceRight]).two_norm();
    }

  // Clean up
  delete[] xCenter;

  // Create the diagram
  m_Diagram = new VoronoiDiagram(
    m_NumberOfVertices, m_AdjacencyIndex, m_Adjacency, m_FullEdgeWeights);
}

void 
VTKMeshVoronoiDiagram
::SetBarrierEdges(const list<Edge> &edges)
{
  // Copy data from raw edges
  memcpy(m_FullEdgeWeights, m_RawEdgeWeights, sizeof(float) * m_NumberOfEdges);

  // Assign the final edge weights to infinity
  list<Edge>::const_iterator it = edges.begin();
  while(it != edges.end())
    {
    unsigned int iDual = m_DualEdgeMap[*it];
    unsigned int iDualOpp = m_DualEdgeMap[m_HalfEdge->GetHalfEdgeOpposite(*it)];

    m_FullEdgeWeights[iDual] = m_FullEdgeWeights[iDualOpp] = 
      VoronoiDiagram::INFINITE_WEIGHT;
    }
}

void 
VTKMeshVoronoiDiagram
::ComputeVoronoiDiagram(const list<vtkIdType> &lSources) 
{
  // Create an array from the list
  unsigned int *lSourceArray = new unsigned int[lSources.size()];
  unsigned int iSource = 0;

  list<vtkIdType>::const_iterator it = lSources.begin();
  while(it != lSources.end())
    lSourceArray[iSource++] = *it++;

  m_Diagram->ComputePathsFromManySources(iSource, lSourceArray);
  
  delete[] lSourceArray;
  
  // There is now a valid diagram
  m_IsDiagramValid = true;
}

void 
VTKMeshVoronoiDiagram
::DeleteGraphData()
{
  if(m_Diagram)
    {
    delete m_Diagram; m_Diagram = NULL;
    delete m_Adjacency;
    delete m_AdjacencyIndex;
    delete m_RawEdgeWeights;
    delete m_FullEdgeWeights;
    }
}

