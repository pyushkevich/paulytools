#ifndef __VTKMeshVoronoiDiagram_
#define __VTKMeshVoronoiDiagram_

#include <vnl/vnl_vector_fixed.h>
#include "vtkPolyData.h"
#include "VTKMeshHalfEdgeWrapper.h"
#include "ShortestPath.h"

#include <list>
#include <map>

using namespace std;

/**
 * This class is used to construct a Voronoi diagram on a mesh. The voronoi
 * diagram is computed on a graph whose vertices are cells. The weight 
 * associated with each edge is the distance between cell centers
 */
class VTKMeshVoronoiDiagram
{
public:
  typedef vnl_vector_fixed<double, 3> Vec;
  typedef GraphVoronoiDiagram<float> VoronoiDiagram;
  // typedef pair<vtkIdType, vtkIdType> Edge;
  typedef unsigned int Edge;
  
  /** Constructor */
  VTKMeshVoronoiDiagram();
  ~VTKMeshVoronoiDiagram();

  /** Set the input mesh (must be triangular, clean, etc) */
  void SetInputMesh(VTKMeshHalfEdgeWrapper *halfEdge);

  /** Compute the graph structure */
  void ComputeGraph();

  /** Add barrier edges to the graph structure */
  void SetBarrierEdges(const list<Edge> &edges);

  /** Set the set of cells that are the 'sources' */
  void ComputeVoronoiDiagram(const list<vtkIdType> &lSources);

  /** Get the id of the closest source cell for a given cell */
  vtkIdType GetVertexSource(vtkIdType iCell) const
    { return m_Diagram->GetVertexSource(iCell); }

  /** Get the voronoi diagram computer */
  VoronoiDiagram *GetDiagram() const
    { return m_Diagram; }

  /** Check whether the diagram is valid (has been computed) */
  bool IsDiagramValid() const
    { return m_IsDiagramValid; }
    
private:
  // The VTK mesh
  vtkPolyData *m_Mesh;

  // The half-edge wrapper
  VTKMeshHalfEdgeWrapper *m_HalfEdge;

  // Clean up graph structures
  void DeleteGraphData();

  // Adjacency table used to represent the edge graph
  unsigned int *m_AdjacencyIndex, *m_Adjacency; 

  // Mesh dimensions
  unsigned int m_NumberOfEdges, m_NumberOfVertices;

  // The array of weights of the graph edges based on distance between
  // the centers of triangles
  float *m_RawEdgeWeights;

  // The array of adjusted edge weights, based on partition curves
  float *m_FullEdgeWeights;

  // The structure used to compute the shortest paths on the mesh
  VoronoiDiagram *m_Diagram;

  // Whether the diagram has been computed and is valid
  bool m_IsDiagramValid;

  // A mapping from edges in dual graph to the edges in the primal graph
  unsigned int *m_DualEdgeMap;
};



#endif // __VTKMeshVoronoiDiagram_
