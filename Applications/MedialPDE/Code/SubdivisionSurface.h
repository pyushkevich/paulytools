#ifndef __SubdivisionSurface_h_
#define __SubdivisionSurface_h_

#include <vector>
#include "smlmath.h"
#include "Registry.h"
#include "SparseMatrix.h"

class vtkPolyData;

using namespace std;

// Useful global-level inline routines.
// TODO: are there faster ops for this?
inline short ror(short x) { return (x + 1) % 3; }
inline short rol(short x) { return (x + 2) % 3; }

// Constant used for unreasonable size_t objects
#define NOID 0xffffffff

/**
 * Subdivision surface representation 
 */
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
    bool IsBoundary() const { return bnd; }

    // Get the valence of the vertex. For internal vertices it is equal to the
    // number of triangles that share the vertex, but for boundary vertices it
    // is equal to 1 plus that (when the mesh has a disk topology)
    size_t Valence() const { return bnd ? n + 1 : n; } 
    
    // Constructors
    Vertex(size_t it0, short iv0, size_t it1, short iv1, size_t in, bool ibnd) 
      : t0(it0), v0(iv0), t1(it1), v1(iv1), n(in), bnd(ibnd) {}
    Vertex() : t0(NOID), t1(NOID), v0(-1), v1(-1), n(0), bnd(false) {}
  };

  // Information describing vertex neighborhood relationship
  struct NeighborInfo
  {
    // Triangle in front and behind
    size_t tFront, tBack;

    // Index of the vertex relative to these triangles
    short vFront, vBack;

    // Constructor and destructor
    NeighborInfo() : tFront(NOID), tBack(NOID), vFront(-1), vBack(-1) {}
    NeighborInfo(size_t tf, short vf, size_t tb, short vb)
      : tFront(tf), vFront(vf), tBack(tb), vBack(vb) {}
  };

  struct MeshLevel
  {
    typedef vector<Triangle>::iterator TriangleIt;

    // List of triangles in this mesh level
    vector<Triangle> triangles;

    // List of vertices in this mesh level
    // vector<Vertex> vertices;

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

    // Sparse matrix representing the neighborhood relationship between
    // vertices in the mesh. This information allows easy walk-line iteration
    // around vertices.
    ImmutableSparseArray<NeighborInfo> nbr;

    // Return the valence of a vertex
    size_t GetVertexValence(size_t ivtx)
      { return nbr.GetRowIndex()[ivtx+1] - nbr.GetRowIndex()[ivtx]; }

    // Check if the vertex is on the boundary
    bool IsVertexInternal(size_t ivtx)
      { return nbr.GetSparseData()[nbr.GetRowIndex()[ivtx]].tBack != NOID; }
    
    // Constructor
    MeshLevel()
      { parent = NULL; nVertices = 0; }
  };

  /** Subdivide a mesh level once */
  static void Subdivide(MeshLevel *src, MeshLevel *dst);

  /** Import a mesh from a VTK mesh that's been wrapped in a half-edge structure */
  static void ImportLevelFromVTK(vtkPolyData *, MeshLevel &dest);

  /** Apply a subdivision to a vtk mesh */
  static void ApplySubdivision(
    vtkPolyData *src, vtkPolyData *target, MeshLevel &m);

  /** Apply subdivision to raw data */
  static void ApplySubdivision(
    SMLVec3d *xSrc, double *rhoSrc, 
    SMLVec3d *xTrg, double *rhoTrg, MeshLevel &m);

  /** Load a mesh level from a registry */
  static void LoadMeshLevel(Registry &registry);

  /** Test the correctness of a mesh level */
  static bool CheckMeshLevel(MeshLevel &mesh);

private:
  // Mutable sparse matrix
  typedef vnl_sparse_matrix<double> MutableSparseMatrix;

  // Vertex assignment function (visit vertices to assign labels)
  static void RecursiveAssignVertexLabel(MeshLevel *mesh, size_t t, size_t v, size_t id);

  // Set the weights for an even vertex
  static void SetEvenVertexWeights(MutableSparseMatrix &W,
    MeshLevel *parent, MeshLevel *child, size_t t, size_t v);

  // Set the weights for an odd vertex
  static void SetOddVertexWeights(MutableSparseMatrix &W,
    MeshLevel *parent, MeshLevel *child, size_t t, size_t v);

  // Compute the walks in a mesh level. This associates each vertex with a
  // triangle, which makes subsequent processing a lot easier
  static void ComputeWalks(MeshLevel *mesh);

};

/**
 * This is an iterator that represents a walk around the vertex where at each
 * point in time you are standing on one edge. To your right is a vertex that
 * you are walking around, to your left is a neighbor vertex. In front of you
 * or behind of you may be a triangle (unless you are at the beginning). The
 * number of steps is equal to the valence of the vertex
 */
class EdgeWalkAroundVertex
{
public:
  /** Constructor: takes a vertex in a fully initialized mesh level */
  EdgeWalkAroundVertex(SubdivisionSurface::MeshLevel *mesh, size_t iVertex)
    { 
    this->mesh = mesh;
    this->iVertex = iVertex;
    this->pos = mesh->nbr.GetRowIndex()[iVertex];
    this->end = mesh->nbr.GetRowIndex()[iVertex + 1];
    }

  /** Check if this walk is closed or open (boundary vertex) */
  bool IsOpen()
    { return mesh->nbr.GetSparseData()[end - 1].tFront == NOID; }

  /** Return the valence of the vertex in the center */
  size_t Valence()
    { return mesh->nbr.GetRowIndex()[iVertex + 1] - mesh->nbr.GetRowIndex()[iVertex]; }

  /** 
   * Are we at the end of the walk? The end is reached when you have made a
   * full circle or moved past the end of the mesh. No operations should be
   * performed once you are at the end of the walk
   */
  bool IsAtEnd() const
    { return pos >= end; }

  /** Increment operator */
  EdgeWalkAroundVertex &operator++ ()
    { pos++; }

  /** Go to the last edge in the walk (next step puts you at end) */
  void GoToLastEdge()
    { pos = end - 1; }

  /** Go to the first edge in the walk (next step puts you at end) */
  void GoToFirstEdge()
    { pos = mesh->nbr.GetRowIndex()[iVertex]; }

  /** Get the id of the 'fixed' vertex. This is always the same id */
  size_t FixedVertexId() 
    { return iVertex; }

  /** Get the moving vertex id */
  size_t MovingVertexId() 
    { return mesh->nbr.GetColIndex()[pos]; }

  /** Get the triangle ahead */
  size_t TriangleAhead() 
    { return mesh->nbr.GetSparseData()[pos].tFront; }

  /** Get the triangle behind */
  size_t TriangleBehind() 
    { return mesh->nbr.GetSparseData()[pos].tBack; }

  /** Get the index of the fixed vertex in the triangle ahead */
  short FixedVertexIndexInTriangleAhead() 
    { return mesh->nbr.GetSparseData()[pos].vFront; }
  
  /** Get the index of the fixed vertex in the triangle behind */
  short FixedVertexIndexInTriangleBehind() 
    { return mesh->nbr.GetSparseData()[pos].vBack; }

  /** Get the index of the moving vertex in the triangle ahead */
  short MovingVertexIndexInTriangleAhead() 
    { return ror(FixedVertexIndexInTriangleAhead()); }

  /** Get the index of the moving vertex in the triangle behind */
  short MovingVertexIndexInTriangleBehind() 
    { return rol(FixedVertexIndexInTriangleBehind()); }

  /** Get the index of the face-connected vertex in triangle ahead */
  short OppositeVertexIndexInTriangleAhead() 
    { return rol(FixedVertexIndexInTriangleAhead()); }

  /** Get the index of the face-connected vertex in triangle behind */
  short OppositeVertexIndexInTriangleBehind() 
    { return ror(FixedVertexIndexInTriangleBehind()); }

  /** Get the id of the face-connected vertex ahead */
  size_t VertexIdAhead()
    { 
    if(pos == end - 1) 
      if(mesh->nbr.GetSparseData()[pos].tFront == NOID)
        return NOID;
      else 
        return mesh->nbr.GetColIndex()[mesh->nbr.GetRowIndex()[iVertex]];
    else 
      return mesh->nbr.GetColIndex()[pos + 1];
    }
    
  /** Get the id of the face-connected vertex behind */
  size_t VertexIdBehind()
    { 
    if(pos == mesh->nbr.GetRowIndex()[iVertex]) 
      if(mesh->nbr.GetSparseData()[pos].tBack == NOID)
        return NOID;
      else 
        return mesh->nbr.GetColIndex()[end - 1];
    else 
      return mesh->nbr.GetColIndex()[pos - 1];
    }

  /** Get the index of the iterator's position in the sparse array */
  size_t GetPositionInMeshSparseArray() const
    { return pos; }

private:
  // The mesh being iterated around
  SubdivisionSurface::MeshLevel *mesh;
  
  // The index of the vertex walked around by this iterator
  size_t iVertex;
  
  // The position in the walk, points to the column entry in the mesh's nbr
  // sparse matrix
  size_t pos, end;
};

inline ostream &operator << (ostream &out, const SubdivisionSurface::Triangle &t)
{
  out << "[V=(" << t.vertices[0] << "," << t.vertices[1] << "," << t.vertices[2] << ")";
  out << " N=(" << t.neighbors[0] << "," << t.neighbors[1] << "," << t.neighbors[2] << ")";
  out << " NE=(" << t.nedges[0] << "," << t.nedges[1] << "," << t.nedges[2] << ")]";
}

#endif
