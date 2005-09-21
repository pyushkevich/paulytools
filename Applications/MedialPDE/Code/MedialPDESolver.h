#ifndef _MedialPDESolver_h_
#define _MedialPDESolver_h_

#include <iostream>
#include <smlmath.h>
#include "MedialPDEMasks.h"
#include "MedialPDESites.h"
#include "MedialAtom.h"
#include "CartesianMedialAtomGrid.h"
#include "BasisFunctions2D.h"
#include "PardisoInterface.h"

using namespace std;

class MedialPDESolver
{
public:
  // Vector and matrix typedefs
  typedef vnl_vector<double> Vec;
  typedef vnl_matrix<double> Mat;
  
  /** Use the parameter grid to initialize the medial PDE solver. */
  MedialPDESolver(const Vec &uGrid, const Vec &vGrid);

  /**
   * This helper constructor generates a medial PDE grid with given numbers
   * (nu, nv) of regularly spaced grid points. In addition, there is an option
   * to generate irregularly spaced grid points at the boundary of the domain.
   * These grid points have exponentially decreasing spacing as you approach
   * the boundary and are used to make sure that the m-rep boundary is sampled
   * densely near the endcaps. The exponent is currenly 2, but that may change
   */
  MedialPDESolver(size_t nu, size_t nv, 
    double xScale = 0.0, size_t pu = 0, size_t pv = 0);

  /** 
   * Store the current solution vector as the initialization vector
   */
  void SetSolutionAsInitialGuess();

  /** 
   * Generate the default initial guess
   */
  void SetDefaultInitialGuess(double xMagnitude = 0.01);

  /** Specify the surface for which the problem should be solved */
  void SetMedialSurface(IMutableHyperSurface2D *xSurface);

  /** Get the pointer to the internal surface */
  IMutableHyperSurface2D *GetMedialSurface() 
    { return xSurface; }

  /**
   * Solve the PDE for a given surface and a given rho function up to the required 
   * level of accuracy.
   */
  void Solve(double delta = 1e-12);

  /**
   * Compute the 'jet' of the equation with respect to the basis functions
   * that constitute the medial surface. This means solving the PDEs that define
   * the gradient of the phi function with respect to the basis functions, as 
   * well as computing the other partial derivatives */
  void ComputeVariationalDerivative(
    IHyperSurface2D *xVariation, bool flagRhoVariation, MedialAtom *dAtoms);

  /** Get the array of atoms */
  MedialAtom *GetAtomArray()
    { return xAtoms; }

  /** Get the atom grid on which the medial atoms are computed */
  MedialAtomGrid *GetAtomGrid()
    { return xGrid; }

  /** Get the number of atoms */
  unsigned int GetNumberOfAtoms()
    { return nSites; }

  /** Get the number of atoms in u dimension */
  unsigned int GetNumberOfUPoints()
    { return m; }

  /** Get the number of atoms in u dimension */
  unsigned int GetNumberOfVPoints()
    { return n; }

  /** Get the radius field */
  const Mat &GetPhiField() { return y; }

  /** Get the u and v finite difference grids */
  const Vec &GetGridU() const { return uGrid; }
  const Vec &GetGridV() const { return vGrid; }

  void TestFiniteDifferenceConvergence();

private:

  // Numbers of grid points
  size_t m, n;
  
  // U and V grids
  Vec uGrid, vGrid;

  /** Number of 'sites' */
  unsigned int nSites;

  /** The hypersurface on which the PDE is solved */
  IMutableHyperSurface2D *xSurface;

  // Array of mask pointers
  std::vector<FiniteDifferenceMask *> xMasks;

  /** Array of sites */
  std::vector<FDAbstractSite *> xSites;

  /** Index into the sites */
  vnl_matrix<size_t> xSiteIndex;

  /** Representation for the sparse matrix */
  double *xSparseValues;
  int *xRowIndex, *xColIndex, nSparseEntries;

  /** A grid representing the structure of the atoms */
  CartesianMedialAtomGrid *xGrid;
  
  /** Array of medial atoms, final solution */
  MedialAtom *xAtoms;

  /** Three vectors used in the iteration */
  Mat eps, b, y, zTest, xInitSoln, xDefaultInitSoln;
  // vnl_vector<double> eps, b, y, zTest;

  // Sparse linear matrix solver
  UnsymmetricRealPARDISO xPardiso;

  // Common initialization
  void Initialize(const Vec &uGrid, const Vec &vGrid);

  /** Routine to compute medial atoms */
  void TestJacobi();
  void ReconstructAtoms(const Mat &ySolution);
  void InitializeSiteGeometry();
  double SolveOnce(double delta);

  bool flagReuseLastSolution;
};

/* *********************** Template Code ************************** */


#endif // _MedialPDESolver_h_
