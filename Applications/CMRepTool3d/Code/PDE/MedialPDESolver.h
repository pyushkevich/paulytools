#ifndef _MedialPDESolver_h_
#define _MedialPDESolver_h_

#include <iostream>
#include <smlmath.h>
#include "MedialAtom.h"
#include "CartesianMedialAtomGrid.h"
#include "BasisFunctions2D.h"

using namespace std;

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

  // Each site represents an equation F[N_i] = 0. This method computes F[N_i]
  double ComputeEquation(double *Y);

  // This function computes dF[N_i]/d[N_i]
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

  friend class MedialPDESolver;
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
  double SolveOnce(double delta);

  bool flagReuseLastSolution;
};

/* *********************** Template Code ************************** */


#endif // _MedialPDESolver_h_
