#ifndef _MedialPDESolver_h_
#define _MedialPDESolver_h_

class GeometryDescriptor
{
public:
  double xCovariantTensor[2][2];
  double xContravariantTensor[2][2];
  double xChristoffelFirst[2][2][2];
  double xChristoffelSecond[2][2][2];
  double xNormal;

  GeometryDescriptor(double *X, double *Xu, double *Xv, double *Xuu, double *Xuv, double *Xvv);
};

class FDAbstractSite
{
public:
  virtual double ComputeEquation(double *X) = 0;
  virtual void ComputeDerivative(double *X, double *A) = 0;
  virtual double GetInitialValue() = 0;
  
protected:
};

class FDInternalSite : public FDAbstractSite
{
public:
  FDInternalSite(unsigned int m, unsigned int n, unsigned int i, unsigned int j);
  double ComputeEquation(double *X);
  void ComputeDerivative(double *X, double *A);
  double GetInitialValue() { return 2.0; }

protected:
  /** Indices into the flat data array of all neighbors */
  unsigned int i0, i1, i2, i3, i4, i5, i6, i7, i8;

  /** Coefficients of the finite differences in the equation */
  double Cuu, Cuv, Cvv, Cu, Cv;

  /** Derivatives with respect to each site */
  double b0, b1, b2, b3, b4, b5, b6, b7, b8;

  /** The value of rho */
  double rho;
};

class FDBorderSite : public FDAbstractSite
{
public:
  FDBorderSite(unsigned int m, unsigned int n, unsigned int i, unsigned int j);
  double ComputeEquation(double *X) { return 0.0; }
  void ComputeDerivative(double *X, double *A) { return; }
  double GetInitialValue() { return 2.0; }
  
protected:
};

class FDCornerSite : public FDAbstractSite
{
public:
  FDCornerSite(unsigned int m, unsigned int n, unsigned int i, unsigned int j);
  double ComputeEquation(double *X) { return 0.0; }
  void ComputeDerivative(double *X, double *A) { return; }
  double GetInitialValue() { return 2.0; }
  
protected:
};

template <class TSurface, TLaplacian>
class MedialPDESolver {
public:
  /**
   * Initialize the solver with a given number of points in n and m. This breaks up the 
   * unit square into (m-1) by (n-1) rectangles
   */
  MedialPDESolver(unsigned int m, unsigned int n);

  /**
   * Solve the PDE for a given surface and a given rho function 
   */
  void Solve(TSurface &surface, TLaplacian &lap);

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
};


#endif // _MedialPDESolver_h_
