#include "VTKMeshVoronoiDiagram.h"
#include "vtkIdList.h"
#include <vector>
#include <map>

using namespace std;

VTKMeshVoronoiDiagram
::VTKMeshVoronoiDiagram()
{
  // Initialize the graph to NULL
  m_Diagram = NULL;
  m_Adjacency = m_AdjacencyIndex = NULL;
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
::SetInputMesh(vtkPolyData *mesh)
{
  // Store the input mesh
  m_Mesh = mesh;

  // Make sure that the links and the cells are built
  m_Mesh->BuildCells();
  m_Mesh->BuildLinks();
  
  // Set the number of vertices to the number of cells in the mesh
  m_NumberOfVertices = m_Mesh->GetNumberOfCells();
  
  // There is no valid diagram now
  m_IsDiagramValid = false;
}

void
VTKMeshVoronoiDiagram
::ComputeGraph()
{
  vtkIdType iCell, iEdge;
  
  // Clean up the old graph info
  DeleteGraphData();

  // Allocate the graph adjacency index info
  m_AdjacencyIndex = new unsigned int[m_NumberOfVertices + 1];

  // Create an array to hold cell center points
  Vec *xCenter = new Vec[m_NumberOfVertices];

  // Create a vector to collect adjacencies
  vector<vtkIdType> lAdjacency;

  // Create a structure to hold neighbor info
  vtkIdList *lNeighbors = vtkIdList::New();

  // Iterate over all cells in the mesh
  for(iCell = 0; iCell < m_Mesh->GetNumberOfCells(); iCell++)
    {
    // Record the next adjacency index
    m_AdjacencyIndex[iCell] = lAdjacency.size();

    // Get the points belonging to the cell
    vtkIdType nPoints, *lPoints;
    m_Mesh->GetCellPoints(iCell,nPoints,lPoints);

    // For each point, look at the cells that the point belongs to and 
    // add them to a counter.
    typedef list<vtkIdType> IdList;
    typedef map<vtkIdType, IdList> CellMap;
    CellMap xCellCount;

    // Compute the center of the cell and add adjacencies
    xCenter[iCell].fill(0.0);
    for(unsigned int iPoint = 0; iPoint < nPoints; iPoint++)
      {
      // Add point coordinate to center
      double *xPoint = m_Mesh->GetPoint(lPoints[iPoint]);
      xCenter[iCell][0] += xPoint[0];
      xCenter[iCell][1] += xPoint[1];
      xCenter[iCell][2] += xPoint[2];

      // Find cells around the current point
      unsigned short nCells;
      vtkIdType *lCells;
      m_Mesh->GetPointCells(lPoints[iPoint], nCells, lCells);

      // Add these cells to the multi-map
      for(unsigned int iNbr = 0; iNbr < nCells; iNbr++)
        xCellCount[lCells[iNbr]].push_back(lPoints[iPoint]);
      }

    // At this point, the map xCellCount has a number associated with each cell
    // in the neighborhood, telling us how many points a nearby cell shares with
    // this cell. 

    // Find the cells that have two points shared with the current cell
    CellMap::iterator it = xCellCount.begin();
    while(it != xCellCount.end())
      {
      if(it->second.size() == 2)
        {
        Edge edge;
        if(it->second.front() < it->second.back())
          {
          edge.first = it->second.front();
          edge.second = it->second.back();
          }
        else
          {
          edge.second = it->second.front();
          edge.first = it->second.back();
          }
        m_EdgeMap.insert(make_pair(edge,lAdjacency.size()));
        lAdjacency.push_back(it->first);
        }
      ++it;
      }

    // Scale the center by the number of points
    xCenter[iCell] *= ( 1.0 / nPoints );
    }
  
  // Record the last adjacency index
  m_AdjacencyIndex[m_NumberOfVertices] = lAdjacency.size();
  cout << " Number of edges (2) : " << m_AdjacencyIndex[m_NumberOfVertices] << endl;

  // Compute the edge weights and store adjacencies as an array
  m_NumberOfEdges = lAdjacency.size();
  m_Adjacency = new unsigned int[m_NumberOfEdges];
  m_RawEdgeWeights = new float[m_NumberOfEdges];
  m_FullEdgeWeights = new float[m_NumberOfEdges];

  // Populate the edge data
  for(iCell = 0; iCell < m_NumberOfVertices; iCell++)
    {
    for(iEdge = m_AdjacencyIndex[iCell]; 
      iEdge < m_AdjacencyIndex[iCell+1]; iEdge++)
      {
      // Copy the adjacency value
      m_Adjacency[iEdge] = lAdjacency[iEdge];

      // Compute the distance between centers
      m_RawEdgeWeights[iEdge] = (float)
        (xCenter[iCell] - xCenter[m_Adjacency[iEdge]]).two_norm();
      }
    }

  // Clean up
  delete[] xCenter;
  lNeighbors->Delete();

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
    EdgeMap::iterator itMap = m_EdgeMap.find(*it);
    if(itMap != m_EdgeMap.end())
      {
      while(itMap != m_EdgeMap.upper_bound(*it))
        {
        cout << " Barrier edge [" << it->first << "," << it->second << "]" << " maps to " << itMap->second << endl;
        m_FullEdgeWeights[itMap->second] = VoronoiDiagram::INFINITE_WEIGHT;
        ++itMap;
        }
      }
    ++it;
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

