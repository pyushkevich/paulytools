#ifndef __VTKMeshShortestDistance_h_
#define __VTKMeshShortestDistance_h_

#include <vtkCellLocator.h>
#include <vtkPolyData.h>
#include <vtkPointLocator.h>
#include <vtkExtractEdges.h>

#include <vnl/vnl_vector_fixed.h>

#include <vector>
#include <utility>
#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace boost;
using namespace std;

/***************************************************************************
 * Create a BOOST graph structure from the white matter mesh
 **************************************************************************/
class VTKMeshShortestDistance 
{
public:
  // type definitions
  typedef vnl_vector_fixed<double,3> Vec;
  
  /** Set the input mesh */
  VTKMeshShortestDistance(vtkPolyData *mesh)
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
      m_Edges.begin(), m_Edges.end(), m_EdgeWeights.begin(), m_EdgePolys->GetNumberOfPoints());
    }

  ~VTKMeshShortestDistance()
    {
    delete m_Graph;
    fltLocator->Delete();
    fltCellLocator->Delete();
    fltEdge->Delete();
    }

  /** Compute shortest distances from a vertex on a mesh to other vertices */
  void ComputeDistances(vtkIdType iStartNode)
    {
    m_Distance.resize(num_vertices(*m_Graph));
    m_Predecessor.resize(num_vertices(*m_Graph));
    VertexDescriptor start = vertex(iStartNode, *m_Graph);

    dijkstra_shortest_paths(
      *m_Graph,start,predecessor_map(&m_Predecessor[0]).distance_map(&m_Distance[0]));
    }

  /** Get the distance between start node and given node */
  unsigned int GetVertexDistance(vtkIdType iNode)
    {
    return m_Distance[iNode];
    }

  /** Use this to get the path between start node and given node */
  vtkIdType GetVertexPredecessor(vtkIdType iNode)
    {
    return m_Predecessor[iNode];
    }

  /** This is a helper method: find vertex whose Euclidean distance to a
    given point is minimal */
  vtkIdType FindClosestVertexInSpace(Vec vec)
    {
    return fltLocator->FindClosestPoint(vec.data_block());
    }

  /** Get the edge mesh to which the indices map */
  vtkPolyData *GetEdgeMesh() 
    {
    return m_EdgePolys;
    }
  
  /** Given a ray, find a point closest to that ray */
  vtkIdType PickPoint(Vec xStart, Vec xEnd)
    {
    Vec ptLine, pCoords;
    double t; 
    vtkIdType subId;
    
    // Compute the intersection with the line
    fltCellLocator->IntersectWithLine(
      xStart.data_block(),xEnd.data_block(),
      0.001, t, ptLine.data_block(), pCoords.data_block(), subId);

    // Find the vertex closest to the intersection
    return fltLocator->FindClosestPoint(ptLine.data_block());
    }

private:

  // Graph-related typedefs (boost)
  typedef property<edge_weight_t, unsigned int> WeightType;
  typedef adjacency_list<vecS, vecS, bidirectionalS, no_property, WeightType > GraphType;
  typedef std::pair<vtkIdType,vtkIdType> EdgeType;
  typedef graph_traits<GraphType>::vertex_descriptor VertexDescriptor;

  // Graph-related attributes: list of edges
  vector<EdgeType> m_Edges;
  vector<unsigned int> m_EdgeWeights;

  // The graph based on the mesh
  GraphType *m_Graph;

  // The distance and predecessor maps
  vector<int> m_Distance;
  vector<VertexDescriptor> m_Predecessor;

  // VTK filters
  vtkExtractEdges *fltEdge;
  vtkPointLocator *fltLocator;
  vtkCellLocator *fltCellLocator;

  // VTK edge poly-data
  vtkPolyData *m_EdgePolys;

  // VTK source poly-data
  vtkPolyData *m_SourcePolys;
};

#endif
