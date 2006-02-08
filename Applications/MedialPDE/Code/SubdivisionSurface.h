#ifndef __SubdivisionSurface_h_
#define __SubdivisionSurface_h_

#include <vector>
#include "Registry.h"
#include "SparseMatrix.h"

class vtkPolyData;

using namespace std;

class SubdivisionSurface 
{
public:

  struct Triangle
  {
    size_t vertices[3];
    size_t neighbors[3];
    short nedges[3];

    // Initializes to dummy values
    Triangle();
  };

  /**
   * Each vertex is associated with a walk. A walk starts at some triangle and
   * ends at a different triangle. We store the starting and ending points of 
   * the walk, as well as the number of steps in the walk
   */
  struct Vertex
  {
    // First triangle in the walk and this vertex's index in it
    size_t t0; short v0;
    
    // Second triangle in the walk and this vertex's index in it
    size_t t1; short v1;

    // The size of the walk (number of triangles) - equal to valence for
    // internal vertices but not for boundary ones!!!
    size_t n;

    // Whether the vertex is on the boundary.
    bool bnd;

    // Is this a boundary vertex
    bool IsBoundary() { return bnd; }
    
    // Constructors
    Vertex(size_t it0, short iv0, size_t it1, short iv1, size_t in, bool ibnd) 
      : t0(it0), v0(iv0), t1(it1), v1(iv1), n(in), bnd(ibnd) {}
    Vertex() : t0(NOID), t1(NOID), v0(-1), v1(-1), n(0), bnd(false) {}
  };

  struct MeshLevel
  {
    typedef vector<Triangle>::iterator TriangleIt;

    // List of triangles in this mesh level
    vector<Triangle> triangles;

    // List of vertices in this mesh level
    vector<Vertex> vertices;

    // Number of vertices at this mesh level
    size_t nVertices;

    // Parent mesh level
    MeshLevel *parent;

    // Representation of vertex walks in the mesh level. This representation
    // associates each vertex with a set of vertices that form a loop around
    // it or, if the vertex is on the boundary, it gives an open walk
    size_t *xWalkIndex, *xWalk;
    
    // Sparse matrix mapping the mesh level to the parent level
    ImmutableSparseMatrix<double> weights;

    // Constructor
    MeshLevel()
      { parent = NULL; nVertices = 0; }
  };

  // A vertex walk is a walk around a vertex 'over the faces' of triangles.
  // The walk has a consistant direction on the mesh and it terminates before
  // reaching a boundary or making a loop. The 'center of the walk' is the 
  // vertex that is being walked about. The other two vertices in the current
  // triangle are called to and from, indicating which is visited first in the
  // walk direction.
  class VertexWalkIterator {
  public:
    // Constructor: takes, mesh starting triangle and vertex
    VertexWalkIterator(MeshLevel *inMesh, size_t iVertex)
      : vtx(inMesh->vertices[iVertex]), mesh(inMesh)
      { 
      T = &mesh->triangles[vtx.t0]; n = vtx.n - 1; t = vtx.t0; v = vtx.v0;
      }

    // Check whether the vertex is internal
    bool IsBoundary()
      { return vtx.bnd; }

    // Move the iterator to the end
    void GoToEnd()
      { n = 0; t = vtx.t1; v = vtx.v1; T = &mesh->triangles[t]; }

    // Move the iterator to the beginning
    void GoToBegin()
      { n = vtx.n - 1; } 

    // Get the size of the walk
    size_t Size()
      { return vtx.n; }
    
    // Walk forward to the next triangle in the walk. Check end condition
    // first; walks around internal vertices do not loop around!
    VertexWalkIterator &operator++()
      { 
      n--;
      short ifrom = FromVertexIndex();
      t = T->neighbors[ifrom];
      v = (T->nedges[ifrom] + 1) % 3;
      T = &mesh->triangles[t];
      return *this;
      }

    // Check if we've reached the end of the walk
    bool IsAtEnd()
      { return n < 0; }

    // Get the triangle that we are standing on right now
    size_t Triangle()
      { return t; }
    
    // Get the index of the center vertex in the current triangle (not its
    // vertex id, this is a number between 0 and 2
    short CenterVertexIndex()
      { return v; }

    // Get the vertex that we are walking from
    short FromVertexIndex()
      { return (v + 1) % 3; }

    // Get the vertex that we are walking towards
    short ToVertexIndex()
      { return (v + 2) % 3; }

    // Get the vertex ids (1 .. n)
    size_t CenterVertexId() { return T->vertices[v]; }
    size_t FromVertexId() { return T->vertices[FromVertexIndex()]; }
    size_t ToVertexId() { return T->vertices[ToVertexIndex()]; }
    
    // Get the next triangle in the walk (can be NOID)
    size_t NextTriangle()
      { return T->neighbors[FromVertexIndex()]; }
      
    // Get the previous triangle in the walk (can be NOID)
    size_t PreviousTriangle()
      { return T->neighbors[ToVertexIndex()]; }
    
  private:
    // A reference to the vertex about which we are walking
    Vertex &vtx;
    
    // The current triangle in the walk
    size_t t;

    // The index of the center vertex in the current triangle
    size_t v;

    // The number of triangles left in the walk
    int n;

    // A pointer to the triangle we are on
    SubdivisionSurface::Triangle *T;

    // A mesh level
    MeshLevel *mesh;
  };

  /** Subdivide a mesh level once */
  void Subdivide(MeshLevel *src, MeshLevel *dst);

  /** Import a mesh from a VTK mesh that's been wrapped in a half-edge structure */
  void ImportLevelFromVTK(vtkPolyData *, MeshLevel &dest);

  /** Apply a subdivision to a vtk mesh */
  void ApplySubdivision(
    vtkPolyData *src, vtkPolyData *target, MeshLevel &m);

  /** Load a mesh level from a registry */
  void LoadMeshLevel(Registry &registry);

  /** Test the correctness of a mesh level */
  bool CheckMeshLevel(MeshLevel &mesh);

private:
  // Mutable sparse matrix
  typedef vnl_sparse_matrix<double> MutableSparseMatrix;

  // Vertex assignment function (visit vertices to assign labels)
  void RecursiveAssignVertexLabel(MeshLevel *mesh, size_t t, size_t v, size_t id);

  // Set the weights for an even vertex
  void SetEvenVertexWeights(MutableSparseMatrix &W,
    MeshLevel *parent, MeshLevel *child, size_t t, size_t v);

  // Set the weights for an odd vertex
  void SetOddVertexWeights(MutableSparseMatrix &W,
    MeshLevel *parent, MeshLevel *child, size_t t, size_t v);

  // Compute the walks in a mesh level. This associates each vertex with a
  // triangle, which makes subsequent processing a lot easier
  void ComputeWalks(MeshLevel *mesh);

  const static size_t NOID;
};

ostream &operator << (ostream &out, const SubdivisionSurface::Triangle &t)
{
  out << "[V=(" << t.vertices[0] << "," << t.vertices[1] << "," << t.vertices[2] << ")";
  out << " N=(" << t.neighbors[0] << "," << t.neighbors[1] << "," << t.neighbors[2] << ")";
  out << " NE=(" << t.nedges[0] << "," << t.nedges[1] << "," << t.nedges[2] << ")]";
}

#endif
