#ifndef _MedialPDESolver_h_
#define _MedialPDESolver_h_

#include <iostream>
#include "mspline.h"

using namespace std;

template <typename TReal> class GeometryDescriptor
{
public:
  TReal xCovariantTensor[2][2];
  TReal xContravariantTensor[2][2];
  TReal xChristoffelFirst[2][2][2];
  TReal xChristoffelSecond[2][2][2];
  TReal xNormal;

  // Initialize the descriptor using a Jet
  void SetJet(TReal *X, TReal *Xu, TReal *Xv, TReal *Xuu, TReal *Xuv, TReal *Xvv);

  // Initialize the jet using medial point
  void SetJet(MedialPoint *mp);

  // Dump information out
  void PrintSelf(ostream &str);

protected:
};

class FDAbstractSite
{
public:
  typedef GeometryDescriptor<float> GDescriptor;
  
  FDAbstractSite(unsigned int m, unsigned int n, unsigned int i, unsigned int j);
  virtual double ComputeEquation(double *Y) = 0;
  virtual void ComputeDerivative(double *Y, double *A, unsigned int iRow) = 0;
  virtual void SetGeometry(GDescriptor *gd, double rho) = 0;

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
  typedef FDAbstractSite::GDescriptor GDescriptor;
  
  FDInternalSite(unsigned int m, unsigned int n, unsigned int i, unsigned int j, 
    unsigned int iNode, unsigned int stride_u, unsigned int stride_v);

  void SetGeometry(GDescriptor *gd, double rho);
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
  typedef FDAbstractSite::GDescriptor GDescriptor;

  FDBorderSite(unsigned int m, unsigned int n, unsigned int i, unsigned int j, 
    unsigned int iNode, unsigned int stride_u, unsigned int stride_v);
  
  double ComputeEquation(double *Y);
  void ComputeDerivative(double *Y, double *A, unsigned int iRow);
  void SetGeometry(GDescriptor *gd, double rho);
  
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
  // Compute the jet at a medial point
  virtual void ComputeJet2(MedialPoint *mp) = 0;

  // Compute the laplacian value
  virtual double ComputeLaplacian(double u, double v) = 0;
};

/** Interface that represents an atom array */
class IMedialSurfacePatch 
{
public:
  virtual unsigned int GetNumberOfAtoms() = 0; 
  virtual MedialPoint *GetAtom(unsigned int iAtom) = 0;
  virtual float IntegrateBoundaryMeasure(BoundaryMeasure *m, float &area) = 0;
  // virtual float IntegrateMedialMeasure(BoundaryMeasure *m, float &area) = 0;
};

class MedialPDESolver : public virtual IMedialSurfacePatch 
{
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
  
  /** Get one of the atoms from the solution */
  MedialPoint *GetAtom(unsigned int i, unsigned int j)
    { return xAtoms +  xSiteIndex[i][j]; }

  /** Get the atom, enumerated by a single index */
  MedialPoint *GetAtom(unsigned int iAtom)
    { return xAtoms + iAtom; }

  /** Integrate a boundary measure */
  virtual float IntegrateBoundaryMeasure(BoundaryMeasure *m, float &area);

  /** Integrate medial measure */
  // virtual float IntegrateMedialMeasure(BoundaryMeasure *m, float &area);

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

  /** Array of sites */
  FDAbstractSite **xSites;

  /** Index into the sites */
  unsigned int **xSiteIndex;

  /** The initial solution */
  double *xInitSoln;

  /** Representation for the sparse matrix */
  double *xSparseValues;
  unsigned int *xRowIndex, *xColIndex, nSparseEntries;
  
  /** Array of medial atoms, final solution */
  MedialPoint *xAtoms;

  /** Array of associated geometry descriptors (stupid) */
  typedef GeometryDescriptor<float> GDescriptor;
  GDescriptor *xGeometryDescriptors;

  /** Three vectors used in the iteration */
  double *eps, *b, *y, *zTest;

  /** Internal data for PARDISO */
  int PT[64], MTYPE, IPARM[64];

  /** Routine to compute medial atoms */
  bool ComputeMedialAtom(MedialPoint *p, GDescriptor *gd);
};

/* *********************** Template Code ************************** */
template<typename TReal>
void
GeometryDescriptor<TReal>
::SetJet(TReal *X, TReal *Xu, TReal *Xv, TReal *Xuu, TReal *Xuv, TReal *Xvv)
{
  // Compute the covariant tensor
  xCovariantTensor[0][0] = Xu[0] * Xu[0] + Xu[1] * Xu[1] + Xu[2] * Xu[2];
  xCovariantTensor[1][1] = Xv[0] * Xv[0] + Xv[1] * Xv[1] + Xv[2] * Xv[2];
  xCovariantTensor[1][0] = Xu[0] * Xv[0] + Xu[1] * Xv[1] + Xu[2] * Xv[2];
  xCovariantTensor[0][1] = xCovariantTensor[1][0];

  // Compute the determinant of the covariant tensor
  TReal g = xCovariantTensor[0][0] * xCovariantTensor[1][1] 
    - xCovariantTensor[0][1] * xCovariantTensor[0][1];
  TReal gInv = 1.0 / g;

  // Compute the contravariant tensor
  xContravariantTensor[0][0] = gInv * xCovariantTensor[1][1];
  xContravariantTensor[1][1] = gInv * xCovariantTensor[0][0];
  xContravariantTensor[0][1] = - gInv * xCovariantTensor[1][0];
  xContravariantTensor[1][0] = xContravariantTensor[0][1];

  // Compute the Christoffel symbols of the first kind
  xChristoffelFirst[0][0][0] = Xuu[0] * Xu[0] + Xuu[1] * Xu[1] + Xuu[2] * Xu[2];
  xChristoffelFirst[0][0][1] = Xuu[0] * Xv[0] + Xuu[1] * Xv[1] + Xuu[2] * Xv[2];
  xChristoffelFirst[0][1][0] = Xuv[0] * Xu[0] + Xuv[1] * Xu[1] + Xuv[2] * Xu[2];
  xChristoffelFirst[0][1][1] = Xuv[0] * Xv[0] + Xuv[1] * Xv[1] + Xuv[2] * Xv[2];
  xChristoffelFirst[1][1][0] = Xvv[0] * Xu[0] + Xvv[1] * Xu[1] + Xvv[2] * Xu[2];
  xChristoffelFirst[1][1][1] = Xvv[0] * Xv[0] + Xvv[1] * Xv[1] + Xvv[2] * Xv[2];
  xChristoffelFirst[1][0][0] = xChristoffelFirst[0][1][0];
  xChristoffelFirst[1][0][1] = xChristoffelFirst[0][1][1];

  // Compute the Christoffel symbols of the second kind
  xChristoffelSecond[0][0][0] = xContravariantTensor[0][0] * xChristoffelFirst[0][0][0] + 
    xContravariantTensor[1][0] * xChristoffelFirst[0][0][1];
  xChristoffelSecond[0][0][1] = xContravariantTensor[0][1] * xChristoffelFirst[0][0][0] + 
    xContravariantTensor[1][1] * xChristoffelFirst[0][0][1];
  
  xChristoffelSecond[0][1][0] = xContravariantTensor[0][0] * xChristoffelFirst[0][1][0] + 
    xContravariantTensor[1][0] * xChristoffelFirst[0][1][1];
  xChristoffelSecond[0][1][1] = xContravariantTensor[0][1] * xChristoffelFirst[0][1][0] + 
    xContravariantTensor[1][1] * xChristoffelFirst[0][1][1];

  xChristoffelSecond[1][1][0] = xContravariantTensor[0][0] * xChristoffelFirst[1][1][0] + 
    xContravariantTensor[1][0] * xChristoffelFirst[1][1][1];
  xChristoffelSecond[1][1][1] = xContravariantTensor[0][1] * xChristoffelFirst[1][1][0] + 
    xContravariantTensor[1][1] * xChristoffelFirst[1][1][1];

  xChristoffelSecond[1][0][0] = xChristoffelSecond[0][1][0];
  xChristoffelSecond[1][0][1] = xChristoffelSecond[0][1][1];
}

template<typename TReal>
void
GeometryDescriptor<TReal>
::SetJet(MedialPoint *mp)
{
  this->SetJet(
    mp->F.data_block(), mp->Fu.data_block(), mp->Fv.data_block(), 
    mp->Fuu.data_block(), mp->Fuv.data_block(), mp->Fvv.data_block()); 
}

template<typename TReal>
void 
GeometryDescriptor<TReal>
::PrintSelf(ostream &str)
{
  str << "CovariantTensor : {{" 
    << xCovariantTensor[0][0] << "," 
    << xCovariantTensor[0][1] << "}, {"
    << xCovariantTensor[1][0] << "," 
    << xCovariantTensor[1][1] << "}}" << endl;
  
  str << "ContravariantTensor : {{" 
    << xContravariantTensor[0][0] << "," 
    << xContravariantTensor[0][1] << "}, {"
    << xContravariantTensor[1][0] << "," 
    << xContravariantTensor[1][1] << "}}" << endl;

  str << "ChristoffelFirst : {{{" 
    << xChristoffelFirst[0][0][0] << ","
    << xChristoffelFirst[0][0][1] << "}, {"
    << xChristoffelFirst[0][1][0] << ","
    << xChristoffelFirst[0][1][1] << "}}, {{"
    << xChristoffelFirst[1][0][0] << ","
    << xChristoffelFirst[1][0][1] << "}, {"
    << xChristoffelFirst[1][1][0] << ","
    << xChristoffelFirst[1][1][1] << "}}}" << endl;

  str << "ChristoffelSecond : {{{" 
    << xChristoffelSecond[0][0][0] << ","
    << xChristoffelSecond[0][0][1] << "}, {"
    << xChristoffelSecond[0][1][0] << ","
    << xChristoffelSecond[0][1][1] << "}}, {{"
    << xChristoffelSecond[1][0][0] << ","
    << xChristoffelSecond[1][0][1] << "}, {"
    << xChristoffelSecond[1][1][0] << ","
    << xChristoffelSecond[1][1][1] << "}}}" << endl;
}



#endif // _MedialPDESolver_h_
