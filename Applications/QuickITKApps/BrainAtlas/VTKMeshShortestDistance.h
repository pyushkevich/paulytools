#ifndef __VTKMeshShortestDistance_h_
#define __VTKMeshShortestDistance_h_

#include <vtkCellLocator.h>
#include <vtkPolyData.h>
#include <vtkPointLocator.h>
#include <vtkFeatureEdges.h>
#include <vtkTriangleFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkLoopSubdivisionFilter.h>

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
  VTKMeshShortestDistance(vtkPolyData *mesh);

  ~VTKMeshShortestDistance();

  /** Compute shortest distances from a vertex on a mesh to other vertices */
  void ComputeDistances(vtkIdType iStartNode);

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
  bool PickPoint(Vec xStart, Vec xEnd, vtkIdType &point);

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
  vtkFeatureEdges *fltEdge;
  vtkTriangleFilter *fltTriangle;
  vtkPointLocator *fltLocator;
  vtkCellLocator *fltCellLocator;
  vtkCleanPolyData *fltCleaner;
  vtkLoopSubdivisionFilter *fltLoop;

  // VTK edge poly-data
  vtkPolyData *m_EdgePolys;

  // VTK source poly-data
  vtkPolyData *m_SourcePolys;
};

#endif
