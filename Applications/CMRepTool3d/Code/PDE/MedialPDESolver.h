#ifndef _MedialPDESolver_h_
#define _MedialPDESolver_h_

#include <iostream>
extern "C" {  
  #include <laspack/qmatrix.h> 
  #include <laspack/vector.h> 
}

using namespace std;

class GeometryDescriptor
{
public:
  double xCovariantTensor[2][2];
  double xContravariantTensor[2][2];
  double xChristoffelFirst[2][2][2];
  double xChristoffelSecond[2][2][2];
  double xNormal;

  GeometryDescriptor(double *X, double *Xu, double *Xv, double *Xuu, double *Xuv, double *Xvv);

  void PrintSelf(ostream &str);

protected:
};

class FDAbstractSite
{
public:
  FDAbstractSite(unsigned int m, unsigned int n, unsigned int i, unsigned int j);
  virtual double ComputeEquation(double *Y) = 0;
  virtual void ComputeDerivative(double *Y, double *A, unsigned int iRow) = 0;
  virtual void SetGeometry(GeometryDescriptor *gd, double rho) = 0;

  // Returns the number of boundaries that a site touches (0 - int, 1 - border, 2 - corner)
  virtual bool IsBorderSite(unsigned int dim) = 0;
  virtual double GetInitialValue();

  // Return the number of columns in each row of the matrix A used in Newton's method
  virtual unsigned int GetNumberOfNetwonColumns() = 0;

  // Return the i-th column number that is non-zero in this equation
  virtual unsigned int GetNewtonColumnAtIndex(unsigned int index) = 0;
  
protected:
  // Position in the grid and grid dimensions
  unsigned int i, j, m, n;

  // Step size (finite difference)
  double du, dv;
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

  // This site touches no boundaries
  bool IsBorderSite(unsigned int) { return false; }

  // Number of non-zero columns in newton's matrix
  unsigned int GetNumberOfNetwonColumns() { return 9; }

  // Get the column of the newton's matrix entry
  unsigned int GetNewtonColumnAtIndex(unsigned int index)
    { return xColumns[index]; }

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
  
  // This site touches one or two boundaries
  bool IsBorderSite(unsigned int dim) 
    { return (dim == 0) ? (wu == 1) : (wv == 1); }

  // Number of non-zero columns in newton's matrix
  unsigned int GetNumberOfNetwonColumns() 
    { return nDistinctSites; }

  // Get the column of the newton's matrix entry
  unsigned int GetNewtonColumnAtIndex(unsigned int index)
    { return xDistinctSites[index]; }

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

/** Abstract problem to be fed to the solver */
class IMedialPDEProblem 
{
public:
  virtual void ComputeJet2(
    double u, double v, double *X, double *Xu, double *Xv, 
    double *Xuu, double *Xuv, double *Xvv) = 0;

  virtual double ComputeLaplacian(double u, double v) = 0;
};

class MedialPDESolver {
public:
  /**
   * Initialize the solver with a given number of points in n and m. This breaks up the 
   * unit square into (m-1) by (n-1) rectangles
   */
  MedialPDESolver(unsigned int m, unsigned int n);

  /**
   * Solve the PDE for a given surface and a given rho function up to the required 
   * level of accuracy.
   */
  void Solve(IMedialPDEProblem *problem, double delta = 1e-8);

private:
  /** Number of points in u and v */
  unsigned int m,n;

  /** Distance between sites in u and v */
  double du, dv;

  /** Number of 'sites' */
  unsigned int nSites;

  /** Array of sites */
  FDAbstractSite **xSites;

  /** Index into the sites */
  unsigned int **xSiteIndex;

  /** The initial solution */
  double *xInitSoln;

  /** Representation for the sparse matrix */
  double *xSparseValues;
  unsigned int *xRowIndex, *xColIndex, nSparseEntries;

  /** Three vectors used in the iteration */
  double *eps, *b, *y;
};



#endif // _MedialPDESolver_h_