#ifndef _MedialPDESolver_h_
#define _MedialPDESolver_h_

#include <iostream>
#include <smlmath.h>
#include "MedialAtomGrid.h"

#include "BasisFunctions2D.h"

using namespace std;

struct GeometryDescriptor
{
  double xCovariantTensor[2][2];
  double xContravariantTensor[2][2];
  double xChristoffelFirst[2][2][2];
  double xChristoffelSecond[2][2][2];

  // The determinant of the covariant tensor and its inverse
  double g, gInv;

  // Initialize the descriptor using a Jet
  void SetJet(double *X, double *Xu, double *Xv, double *Xuu, double *Xuv, double *Xvv);

  // Dump information out
  void PrintSelf(ostream &str);
};

/** 
 * A Boundary atom
 */
struct BoundaryAtom 
{
  SMLVec3d X, N;
};

/**
 * The new medial atom representation
 */
struct MedialAtom
{
  // The coordinates of the atom in the domain
  double u, v;
  
  // The position on the medial surface and corresponding partial derivatives
  SMLVec3d X, Xu, Xv, Xuu, Xuv, Xvv;

  // The differential geometry descriptor of the medial surface
  GeometryDescriptor G;

  // The radius function and its partial derivatives
  double R, Ru, Rv, Ruu, Ruv, Rvv;

  // The normal vector and the Riemannian gradient of R on the surface
  SMLVec3d N, xGradR;

  // The Riemannian laplacian of R
  double xLapR;

  // Whether this is a 'crest' atom, and whether it's valid at all
  bool flagCrest, flagValid;

  // The two associated boundary 'atoms'
  BoundaryAtom xBnd[2];

  /** Compute the differential geometric quantities from X and its 1st and 2nd
   * derivatives. This does not compute the normal vector */
  void ComputeDifferentialGeometry();

  /** Compute the normal vector */
  void ComputeNormalVector();

  /** Given the differential geometry and the normal vector are computed,
   * compute GradR and the boundary sites. */
  bool ComputeBoundaryAtoms();
};

/** Helper function to access a boundary site in an atom array using a
 * boundary point iterator */
inline BoundaryAtom &
GetBoundaryPoint(MedialBoundaryPointIterator *itBoundary, MedialAtom *xAtoms)
{
  return 
    xAtoms[itBoundary->GetAtomIndex()].xBnd[itBoundary->GetBoundarySide()];
}

/** Helper function to access a boundary site in an atom array using a
 * boundary point iterator */
inline BoundaryAtom &
GetBoundaryPoint(MedialBoundaryQuadIterator *itBQuad, MedialAtom *xAtoms, 
  unsigned int i, unsigned int j)
{
  return 
    xAtoms[itBQuad->GetAtomIndex(i,j)].xBnd[itBQuad->GetBoundarySide()];
}

class FDAbstractSite
{
public:
  
  FDAbstractSite(unsigned int m, unsigned int n, unsigned int i, unsigned int j);
  virtual double ComputeEquation(double *Y) = 0;
  virtual void ComputeDerivative(double *Y, double *A, unsigned int iRow) = 0;
  virtual void SetGeometry(GeometryDescriptor *gd, double rho) = 0;

  // Returns the number of boundaries that a site touches (0 - int, 1 - border, 2 - corner)
  virtual bool IsBorderSite() = 0;
  virtual bool IsBorderSite(unsigned int dim) = 0;
  virtual double GetInitialValue();
  
  // Compute the R, Ru and Rv given the solution - used to compute atoms
  virtual void ComputeRJet(double *Y, double &R, double &Ru, double &Rv) = 0;

  // Return the number of columns in each row of the matrix A used in Newton's method
  virtual unsigned int GetNumberOfNetwonColumns() = 0;

  // Return the i-th column number that is non-zero in this equation
  virtual unsigned int GetNewtonColumnAtIndex(unsigned int index) = 0;

  // Compute a Jacobi iteration (this site's R as a function of the neighbors)
  virtual double ComputeJacobiIteration(double *Y) = 0;
  
protected:
  // Position in the grid and grid dimensions
  unsigned int i, j, m, n;

  // Step size (finite difference)
  double du, dv, _du, _dv;
};

/**
 * This site represents internal grid points where the Laplacial equation holds
 */
class FDInternalSite : public FDAbstractSite
{
public:
  FDInternalSite(unsigned int m, unsigned int n, unsigned int i, unsigned int j, 
    unsigned int iNode, unsigned int stride_u, unsigned int stride_v);

  void SetGeometry(GeometryDescriptor *gd, double rho);
  double ComputeEquation(double *Y);
  void ComputeDerivative(double *Y, double *A, unsigned int iRow);

  // Compute the R, Ru and Rv given the solution - used to compute atoms
  void ComputeRJet(double *Y, double &R, double &Ru, double &Rv);

  // This site touches no boundaries
  bool IsBorderSite() { return false; }
  bool IsBorderSite(unsigned int) { return false; }

  // Number of non-zero columns in newton's matrix
  unsigned int GetNumberOfNetwonColumns() { return 9; }

  // Get the column of the newton's matrix entry
  unsigned int GetNewtonColumnAtIndex(unsigned int index)
    { return xColumns[index]; }

  // Compute a Jacobi iteration (this site's R as a function of the neighbors)
  virtual double ComputeJacobiIteration(double *Y);
  
protected:
  /** Indices into the flat data array of all neighbors, indexed counterclockwise
   with the center at index 0 */
  unsigned int i0, i1, i2, i3, i4, i5, i6, i7, i8;

  /** Sorted list of 9 column numbers associated with the neighbors */
  unsigned int xColumns[9];

  /** Mapping from internal neighbor index to column number */
  unsigned int xEntry[9];

  /** Coefficients of the finite differences in the equation */
  double Cuu, Cuv, Cvv, Cu, Cv;

  /** Derivatives with respect to each site */
  double b0, b1, b2, b3, b4, b5, b6, b7, b8;

  /** The value of rho */
  double rho;
};

/**
 * This site represents border and corner grid points where the 
 * boundary condition holds
 */
class FDBorderSite : public FDAbstractSite
{
public:
  FDBorderSite(unsigned int m, unsigned int n, unsigned int i, unsigned int j, 
    unsigned int iNode, unsigned int stride_u, unsigned int stride_v);
  
  double ComputeEquation(double *Y);
  void ComputeDerivative(double *Y, double *A, unsigned int iRow);
  void SetGeometry(GeometryDescriptor *gd, double rho);
  
  // Compute the R, Ru and Rv given the solution - used to compute atoms
  void ComputeRJet(double *Y, double &R, double &Ru, double &Rv);
  
  // This site touches one or two boundaries
  bool IsBorderSite() { return true; } 
  bool IsBorderSite(unsigned int dim) 
    { return (dim == 0) ? (wu > 0.75) : (wv > 0.75); }

  // Number of non-zero columns in newton's matrix
  unsigned int GetNumberOfNetwonColumns() 
    { return nDistinctSites; }

  // Get the column of the newton's matrix entry
  unsigned int GetNewtonColumnAtIndex(unsigned int index)
    { return xDistinctSites[index]; }

  // Compute a Jacobi iteration (this site's R as a function of the neighbors)
  double ComputeJacobiIteration(double *Y);
  
protected:
  /** Is this a corner site? */
  bool corner;

  /** Indices into the flat data array of all neighbors */
  unsigned int xNeighbor[5];

  /** List of distinct sites involved in this site's equation */
  unsigned int nDistinctSites, xDistinctSites[4];

  /** Index of each neighbor into xDistinctSites */
  unsigned int xEntry[5];

  /** Coefficients of the finite differences in the equation */
  double CuCu, CuCv, CvCv, C0;

  /** Weight scalings in X and Y */
  double wu, wv;
};

// This is a measure that can be computed over a volume (just a R3 function)
class EuclideanFunction {
public:
  virtual double Evaluate(const SMLVec3d &x) = 0;
  virtual SMLVec3d ComputeGradient(const SMLVec3d &x)
    { return SMLVec3d(0.0); }
};

class MedialPDESolver
{
public:
  /**
   * Initialize the solver with a given number of points in n and m. This breaks up the 
   * unit square into (m-1) by (n-1) rectangles
   */
  MedialPDESolver(unsigned int m, unsigned int n);

  /** 
   * Store the current solution vector as the initialization vector
   */
  void SetSolutionAsInitialGuess();

  /** 
   * Generate the default initial guess
   */
  void SetDefaultInitialGuess(double xMagnitude = 0.01);

  /** Specify the surface for which the problem should be solved */
  void SetMedialSurface(IHyperSurface2D *xSurface);

  /**
   * Solve the PDE for a given surface and a given rho function up to the required 
   * level of accuracy.
   */
  void Solve(double delta = 1e-8);

  /** Alternative, very slow method to solve the equation */
  void SolveByJacobiMethod(double delta = 1e-8);
  
  /** Get the array of atoms */
  MedialAtom *GetAtomArray()
    { return xAtoms; }

  /** Get the atom grid on which the medial atoms are computed */
  MedialAtomGrid *GetAtomGrid()
    { return xGrid; }

  /** Integrate a boundary measure */
  virtual double IntegrateBoundaryMeasure(EuclideanFunction *m, double &area);

  /** Integrate a volume measure */
  virtual double IntegrateVolumeMeasure(
    EuclideanFunction *m, unsigned int nSamples, double &volume);

  /** Get the number of atoms */
  unsigned int GetNumberOfAtoms()
    { return nSites; }

  /** Get the number of atoms in u dimension */
  unsigned int GetNumberOfUPoints()
    { return m; }

  /** Get the number of atoms in u dimension */
  unsigned int GetNumberOfVPoints()
    { return n; }

private:
  /** Number of points in u and v */
  unsigned int m,n;

  /** Distance between sites in u and v */
  double du, dv;

  /** Number of 'sites' */
  unsigned int nSites;

  /** The hypersurface on which the PDE is solved */
  IHyperSurface2D *xSurface;

  /** Array of sites */
  FDAbstractSite **xSites;

  /** Index into the sites */
  unsigned int **xSiteIndex;

  /** The initial solution */
  double *xInitSoln;

  /** Representation for the sparse matrix */
  double *xSparseValues;
  unsigned int *xRowIndex, *xColIndex, nSparseEntries;

  /** A grid representing the structure of the atoms */
  CartesianMedialAtomGrid *xGrid;
  
  /** Array of medial atoms, final solution */
  MedialAtom *xAtoms;

  /** Three vectors used in the iteration */
  double *eps, *b, *y, *zTest;

  /** The values of parameters u and v at grid points */
  double *xGridU, *xGridV;

  /** Internal data for PARDISO */
  int PT[64], MTYPE, IPARM[64];

  /** Routine to compute medial atoms */
  void TestJacobi();
  void ReconstructAtoms(double *ySolution);
  void InitializeSiteGeometry();
};

/* *********************** Template Code ************************** */

/**
 * Compute an area of a triangle
 */
inline double TriangleArea(const SMLVec3d &A, const SMLVec3d &B, const SMLVec3d &C)
{
  return cross_3d(B - A, C - A).magnitude();
}

/**
 * Compute the volume of a prism formed by two triangles
 */
inline double PrismVolume( SMLVec3d A1, SMLVec3d B1, SMLVec3d C1, 
  SMLVec3d A2, SMLVec3d B2, SMLVec3d C2) 
{
  // TODO: Get/derive the correct formula!!!
  double xArea1 = TriangleArea(A1, B1, C1);
  double xArea2 = TriangleArea(A2, B2, C2);
  double d = ((A1 + B1 + C1) - (A2 + B2 + C2)).magnitude() / 3.0;
  return d * 0.5 * (xArea1 + xArea2);
}


#endif // _MedialPDESolver_h_
