#ifndef __MeshMedialPDESolver_h_
#define __MeshMedialPDESolver_h_

#include "SparseMatrix.h"
#include "SubdivisionSurface.h"
#include <smlmath.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

/**
 * Mesh Medial PDE solver is similar to the regular medial PDE solver with
 * some substantial differences. The domain of the representation is a
 * triangular mesh of arbitrary shape. Finite difference expressions are
 * computed using mesh-based operators rather than differential geometry
 *
 * Solving an equation involves solving a series of linear Newton step 
 * partial differential equations:
 * \Delta \epsilon_i = \rho - \Delta \phi_i, subject to
 * 2 \Nabla \epsilon_i \cdot \Nabla \phi_i - 4 \epsilon_i = 
 *            = 4 \phi_i - \| \Nabla \phi_i \|^2 on \partial \Omega
 * To solve this system, we have to construct a matrix that contains the
 * weight of each \epsilon in each of the equations. This is a sparse matrix
 * where each vertex represents a row, and has non-zero values in the
 * columns that correspond to the adjacent vertices
 */
class MeshMedialPDESolver
{
public:
  typedef ImmutableSparseMatrix<double> SparseMat;
  typedef SubdivisionSurface::MeshLevel MeshLevel;
  typedef vnl_matrix<double> Mat;
  typedef vnl_vector<double> Vec;

  // Constructor
  MeshMedialPDESolver();

  // Destructor
  ~MeshMedialPDESolver();

  // Set the topology of the mesh. This determines the connectivity of the
  // vertices which, in turn, determines the structure of the sparse matrix
  // passed to the sparse solver 
  void SetMeshTopology(MeshLevel *topology);

  // Our first attempt at a solver method
  void SolveEquation(const SMLVec3d *X, const double *rho, const double *phi, double *soln);

private:

  // This method resets all the pointers associated with a mesh
  void Reset();

  // This method is used to compute the sparse matrix A for Newton's method
  void FillNewtonMatrix(const SMLVec3d *X, const double *phi);

  // This method is used to compute the right hand side B for Newton's method
  void FillNewtonRHS(const SMLVec3d *X, const double *rho, const double *phi);

  // This computes the geometry associated with a mesh before running the solver
  void ComputeMeshGeometry(const SMLVec3d *X);

  // An immutable sparse matrix used to represent the PDE. Each row
  // corresponds to a vertex, non-zero values are found between adjacent
  // vertices (those contributing to gradient and LBO computations)
  SparseMat A;

  // Arrays B and epsilon used in Newton's solver
  double *xRHS, *xEpsilon;
  
  // A pointer to the mesh topology
  MeshLevel *topology;

  // A structure representing triangle geometry in the mesh
  struct TriangleGeom
    {
    // Area of the triangle
    double xArea;

    // Cotangent of the three angles
    SMLVec3d xCotangent;
    };
  
  // A structure representing vertex geometry in the mesh (temp)
  struct VertexGeom
    { 
    // Area of the triangle fan around the vertex
    double xFanArea;

    // The two Loop tangent vectors at the vertex
    SMLVec3d t1, t2;

    // The first fundamental form and its inverse
    double gCovariant[2][2], gContravariant[2][2];
    };

  // Geometry arrays that store triangle-related and vertex-related info
  TriangleGeom *xTriangleGeom;
  VertexGeom *xVertexGeom;

  // An array of weights used to compute the tangent vectors at each vertex. These correspond 
  // to the sparse matrix A, i.e., include the neighbors of the vertex and the vertex itself
  double *xTangentWeights[2];
  
  // This utility array maps the entries in the nbr sparse array stored in the
  // mesh to elements of the sparse matrix A. This mapping is needed because
  // the entries of array A are sorted and include diagonal entries
  size_t *xMapVertexNbrToA, *xMapVertexToA;
};



#endif
