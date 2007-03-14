#ifndef __MeshMedialPDESolver_h_
#define __MeshMedialPDESolver_h_

#include "MedialAtom.h"
#include "MedialAtomGrid.h"
#include "SparseMatrix.h"
#include "SubdivisionSurface.h"
#include "GenericMedialModel.h"
#include "PardisoInterface.h"
#include <smlmath.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <gsl/gsl_multimin.h>

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
  // Matrix typedefs
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
  void SetMeshTopology(MeshLevel *topology, MedialAtom *managedAtoms);

  // Set the input data as an array of positions and rho values. Another way
  // to set the input data is to manipulate the values of the atoms directly
  // The third (optional) parameter is the 'guess' solution, which serves as
  // the initial solution for Newton iteration
  void SetInputData(const SMLVec3d *X, const double *rho, const double *phi = NULL);

  // Our first attempt at a solver method. The flag specifies whether the
  // terms involved in gradient computation should also be computed
  void SolveEquation(double *xInitSoln = NULL, bool flagGradient = false);

  // Compute the gradient of the solution with respect to some basis.
  void ComputeGradient(vector<MedialAtom *> dAtoms);

  // Test the accuracy of partial derivative computations in the gradient code
  // this method should be called after calling solve with some data
  int TestPartialDerivatives();

  // Test whether the Jacobian is computed correctly
  int TestJacobianAndGradient(double *xInitSoln);

  // Get the array of atoms
  MedialAtom *GetAtomArray() const
    { return xAtoms; }

  // GSL interaction
  static double gsl_f(const gsl_vector *v, void *params);
  static void gsl_df(const gsl_vector *v, void *params, gsl_vector *df);
  static void gsl_fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df);

private:

  // This method resets all the pointers associated with a mesh
  void Reset();

  // This method is used to compute the sparse matrix A for Newton's method
  void FillNewtonMatrix(const double *phi, bool flagInputChange);

  // This method is used to compute the right hand side B for Newton's method
  double FillNewtonRHS(const double *phi);

  // Compute the value of the function (Lap phi - rho at internal vertex, or
  // \|grad phi\|^2 - 4 phi at boundary vertex)
  double ComputeNodeF(size_t i, const double *phi);

  // This computes the geometry associated with a mesh before running the solver
  void ComputeMeshGeometry(bool flagGradient);

  // This computes the medial atoms once the phi has been solved for
  void ComputeMedialAtoms(const double *soln);

  // Compute the weight matrix used for gradient computations
  void ComputeRHSGradientMatrix();

  // Compute the condition number of the Jacobian
  void ComputeJacobianConditionNumber();

  // An immutable sparse matrix used to represent the PDE. Each row
  // corresponds to a vertex, non-zero values are found between adjacent
  // vertices (those contributing to gradient and LBO computations)
  SparseMat A;

  // Arrays B and epsilon used in Newton's solver
  vnl_vector<double> xRHS, xEpsilon;

  // A pointer to the mesh topology
  MeshLevel *topology;

  // A structure representing triangle geometry in the mesh
  struct TriangleGeom
    {
    // Area of the triangle
    double xArea;

    // Cotangent of the three angles
    double xCotangent[3];

    // The gradient of the area with respect to each vertex
    SMLVec3d xAreaGrad[3];

    // The gradient of the cotangent with respect to each vertex
    SMLVec3d xCotGrad[3][3];
    };

  // A structure representing vertex geometry in the mesh (temp)
  struct VertexGeom
    {
    // Area of the triangle fan around the vertex
    double xFanArea;

    // The two Loop tangent vectors at the vertex
    // SMLVec3d t1, t2;

    // The first fundamental form and its inverse
    // double gCovariant[2][2], gContravariant[2][2];
    };

  // Geometry arrays that store triangle-related and vertex-related info
  TriangleGeom *xTriangleGeom;
  VertexGeom *xVertexGeom;

  // Scheme for computing tangent vectors
  MedialAtomLoopScheme xLoopScheme;

  // This utility array maps the entries in the nbr sparse array stored in the
  // mesh to elements of the sparse matrix A. This mapping is needed because
  // the entries of array A are sorted and include diagonal entries
  size_t *xMapVertexNbrToA, *xMapVertexToA;

  // For gradient computation, there is an array W, which contains the
  // Jacobian of each finite difference equation with respect to each of the
  // neighbor atoms's positions. This is an array of vectors
  SMLVec3d *W;

  // Pardiso-compatible 1-based index into the matrix A (redundant, really)
  int *xPardisoRowIndex, *xPardisoColIndex;

  // Array of medial atoms managed by this solver (should it be managed
  // elsewhere?)
  MedialAtom *xAtoms;

  // Pardiso solver
  UnsymmetricRealPARDISO xPardiso;
};



#endif
