#ifndef __VTKMeshShortestDistance_h_
#define __VTKMeshShortestDistance_h_

#include <vtkCellLocator.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
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

/** 
 * Hierarchy of function objects used to compute edge weights
 * in the shortest path application
 */
class MeshEdgeWeightFunction 
{
public:
  virtual double GetEdgeWeight(
    vtkPolyData *mesh, vtkIdType x1, vtkIdType x2) const = 0;
};

class EuclideanDistanceMeshEdgeWeightFunction : public MeshEdgeWeightFunction
{
public:
  virtual double 
    GetEdgeWeight(vtkPolyData *mesh, vtkIdType x1, vtkIdType x2) const
    {
    // Get the two points
    vnl_vector_fixed<double,3> p1(mesh->GetPoint(x1));
    vnl_vector_fixed<double,3> p2(mesh->GetPoint(x2));

    // Compute the associated weight
    return (p1 - p2).two_norm();
    }
};
  
class PitchBasedMeshEdgeWeightFunction : public MeshEdgeWeightFunction
{
public:
  PitchBasedMeshEdgeWeightFunction() : m_PitchFactor(1.0) {}

  void SetPitchFactor(double c)
    { this->m_PitchFactor = c; };

  double GetPitchFactor() const
    { return m_PitchFactor; }
  
  virtual double 
    GetEdgeWeight(vtkPolyData *mesh, vtkIdType x1, vtkIdType x2) const
    {
    // Get the two points 
    vnl_vector_fixed<double,3> p1(mesh->GetPoint(x1));
    vnl_vector_fixed<double,3> p2(mesh->GetPoint(x2));

    // Get the normals associated with the points
    vnl_vector_fixed<double,3> n1( 
      mesh->GetPointData()->GetNormals()->GetTuple3(x1));
    vnl_vector_fixed<double,3> n2( 
      mesh->GetPointData()->GetNormals()->GetTuple3(x2));

    // Compute the 'pitch', i.e., the change in the normal w.r.t. step
    // along the edge, as projected on the edge
    double xPitch = fabs(dot_product(n2 - n1,p2 - p1));
    double xDistance = (p2 - p1).two_norm();

    // Compute the associated weight of the edge
    return (m_PitchFactor * xPitch + xDistance);
    }
private:
  double m_PitchFactor;
};
  
/***************************************************************************
 * Create a BOOST graph structure from the white matter mesh
 *
 * This class takes as an input a mesh and computes a graph that can
 * be used to calculate Dijkstra's shortest path on the surface. 
 *
 * The newer version of this class does not preprocess the mesh. Do all the
 * preprocessing on the mesh externally!
 *
 * In addition to providing graph shortest path, this class can also be 
 * used to find the closest vertex to a given point in space and to intersect
 * the graph by a ray. Later this functionality oughtta be moved to another 
 * class
 **************************************************************************/
class VTKMeshShortestDistance 
{
public:
  // type definitions
  typedef vnl_vector_fixed<double,3> Vec;
  
  /** Constructor */
  VTKMeshShortestDistance();
  ~VTKMeshShortestDistance();

  /** Specify the mesh to use for computing distances */
  void SetInputMesh(vtkPolyData *mesh);

  /** Specify the edge weight function object to use in order to 
   * compute the costs of edge traversal */
  void SetEdgeWeightFunction(MeshEdgeWeightFunction *function)
    { m_WeightFunctionPtr = function; }

  /** Compute the edge graph */
  void ComputeGraph();

  /** Compute shortest distances from a vertex on a mesh to other vertices */
  void ComputeDistances(vtkIdType iStartNode);

  /** Get the distance between start node and given node */
  double GetVertexDistance(vtkIdType iNode)
    { return m_Distance[iNode] / EDGE_WEIGHT_FACTOR; }

  /** Use this to get the path between start node and given node */
  vtkIdType GetVertexPredecessor(vtkIdType iNode)
    { return m_Predecessor[iNode]; }

  /** This is a helper method: find vertex whose Euclidean distance to a
    given point is minimal */
  vtkIdType FindClosestVertexInSpace(Vec vec)
    { return fltLocator->FindClosestPoint(vec.data_block()); }

  /** Get the edge mesh to which the indices map */
  vtkPolyData *GetInputMesh() 
    { return m_SourceMesh; }

  /** Get the weight associated with an edge */
  double GetEdgeWeight(vtkIdType x1, vtkIdType x2)
    { return m_WeightFunctionPtr->GetEdgeWeight(m_SourceMesh,x1,x2); }
  
  /** Given a ray, find a point closest to that ray */
  bool PickPoint(Vec xStart, Vec xEnd, vtkIdType &point);

  /** Get the number of edges in the graph */
  unsigned int GetNumberOfEdges() const
    { return m_Edges.size(); }

  /** Get the given edge */
  vtkIdType GetEdgeStart(unsigned int iEdge) const
    { return m_Edges[iEdge].first; }

  /** Get the given edge */
  vtkIdType GetEdgeEnd(unsigned int iEdge) const
    { return m_Edges[iEdge].second; }

private:

  // Constant mapping Euclidean (and other) mesh distances to integer
  // values that are used for Dijkstra's algorithm
  static const double EDGE_WEIGHT_FACTOR;

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
  vtkPointLocator *fltLocator;
  vtkCellLocator *fltCellLocator;

  // VTK source poly-data
  vtkPolyData *m_SourceMesh;

  // Function object used to compute edge distances
  MeshEdgeWeightFunction *m_WeightFunctionPtr;

  // Default function object used to compute edge distances
  EuclideanDistanceMeshEdgeWeightFunction m_DefaultWeightFunction;
  PitchBasedMeshEdgeWeightFunction m_DefaultWeightFunction2;
};

#endif
