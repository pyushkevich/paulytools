#include "MedialPDESolver.h"
#include <cmath>
#include <algorithm>

using namespace std;

extern "C" {
  #include <laspack/itersolv.h> 
  #include <laspack/operats.h> 
  #include <laspack/precond.h> 
  #include <laspack/rtc.h> 
}

GeometryDescriptor
::GeometryDescriptor(double *X, double *Xu, double *Xv, double *Xuu, double *Xuv, double *Xvv)
{
  // Compute the covariant tensor
  xCovariantTensor[0][0] = Xu[0] * Xu[0] + Xu[1] * Xu[1] + Xu[2] * Xu[2];
  xCovariantTensor[1][1] = Xv[0] * Xv[0] + Xv[1] * Xv[1] + Xv[2] * Xv[2];
  xCovariantTensor[1][0] = Xu[0] * Xv[0] + Xu[1] * Xv[1] + Xu[2] * Xv[2];
  xCovariantTensor[0][1] = xCovariantTensor[1][0];

  // Compute the determinant of the covariant tensor
  double g = xCovariantTensor[0][0] * xCovariantTensor[1][1] 
    - xCovariantTensor[0][1] * xCovariantTensor[0][1];
  double gInv = 1.0 / g;

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

void 
GeometryDescriptor
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

FDAbstractSite::FDAbstractSite(unsigned int m, unsigned int n, unsigned int i, unsigned int j)
{
  // Store the position of the site
  this->m = m;
  this->n = n;
  this->i = i;
  this->j = j;

  // Compute du and dv - although it is wasteful to repeat this for each site
  du = 1.0 / (m - 1);
  dv = 1.0 / (n - 1);
}

double FDAbstractSite::GetInitialValue()
{
  double x = IsBorderSite(0) ? 0 : du;
  double y = IsBorderSite(1) ? 0 : dv;
  return sqrt(x * x + y * y);
}

/** Initialize the site with grid of size m and n */
FDInternalSite::FDInternalSite(
  unsigned int m, unsigned int n, 
  unsigned int i, unsigned int j,
  unsigned int iNode, 
  unsigned int stride_u, unsigned int stride_v)
: FDAbstractSite(m, n, i, j)
{
  // Set the indices of the eight neighbors in the array, go around a circle
  i0 = iNode;
  i1 = i0 + stride_u; 
  i3 = i0 + stride_v;   
  i5 = i0 - stride_u; 
  i7 = i0 - stride_v;
  i2 = i3 + stride_u; 
  i4 = i3 - stride_u; 
  i6 = i7 - stride_u; 
  i8 = i7 + stride_u;

  // Initialize the indices in sorted order
  if(stride_u < stride_v)
    {
    xColumns[0] = i6; xEntry[6] = 0;
    xColumns[1] = i7; xEntry[7] = 1;
    xColumns[2] = i8; xEntry[8] = 2;
    xColumns[3] = i5; xEntry[5] = 3;
    xColumns[4] = i0; xEntry[0] = 4;
    xColumns[5] = i1; xEntry[1] = 5;
    xColumns[6] = i4; xEntry[4] = 6;
    xColumns[7] = i3; xEntry[3] = 7;
    xColumns[8] = i2; xEntry[2] = 8;
    }
  else
    {
    xColumns[0] = i6; xEntry[5] = 0;
    xColumns[1] = i5; xEntry[4] = 1;
    xColumns[2] = i4; xEntry[3] = 2;
    xColumns[3] = i7; xEntry[7] = 3;
    xColumns[4] = i0; xEntry[0] = 4;
    xColumns[5] = i3; xEntry[3] = 5;
    xColumns[6] = i8; xEntry[8] = 6;
    xColumns[7] = i1; xEntry[1] = 7;
    xColumns[8] = i2; xEntry[2] = 8;
    }
}

/** Initialize the border site */
FDBorderSite::FDBorderSite(
  unsigned int m, unsigned int n, 
  unsigned int i, unsigned int j,
  unsigned int iNode, 
  unsigned int stride_u, unsigned int stride_v)
: FDAbstractSite(m, n, i, j)
{
  // Initialize the cites around a circle, forgetting for a moment that
  // we are at a border site. Some of these indices are off the grid
  xNeighbor[0] = iNode;
  xNeighbor[1] = iNode + stride_u; 
  xNeighbor[2] = iNode + stride_v; 
  xNeighbor[3] = iNode - stride_u; 
  xNeighbor[4] = iNode - stride_v;

  // Initialize the weights in u and v
  wu = 0.5; wv = 0.5;

  // Now, depending on which border we are on, correct the index and weights of finite
  // differences in u and v directions
  if(i == 0)
    { xNeighbor[3] = iNode; wu = 1; }
  else if(i == m - 1)
    { xNeighbor[1] = iNode; wu = 1; }
  if(j == 0)
    { xNeighbor[4] = iNode; wv = 1; }
  else if(j == n - 1)
    { xNeighbor[2] = iNode; wv = 1; }

  // Check whether we are at a corner
  corner = (wu == 1.0 && wv == 1.0);
  nDistinctSites = corner ? 3 : 4;

  // Now, sort the entries in xNeighbor
  unsigned int k, p;
  for(k = 0; k < 4; k++) xDistinctSites[k] = xNeighbor[k+1];
  sort(xDistinctSites,xDistinctSites+4);

  // Two of the distinct sites can be duplicates - they should be removed
  for(k = 0; k < 3; k++)
    if(xDistinctSites[k] == xDistinctSites[k+1])
      for(p = k+1; p < 3; p++)
        xDistinctSites[p] = xDistinctSites[p+1];

  // At this point, the first nDistinctSites in xDistinctSites hold sorted
  // unique column values. Now we need to compute a back-map taking each
  // neighbor to one of these values. For now we do this in a slow manner
  for(k = 0; k < 5; k++)
    for(p = 0; p < nDistinctSites; p++)
      if(xNeighbor[k] == xDistinctSites[p]) 
        { xEntry[k] = p; break; }
}

/** Compute the equation corresponding to an internal site. */
double FDInternalSite::ComputeEquation(double *Y)
{
  // Get the components
  double 
    F0 = Y[i0], F1 = Y[i1], F2 = Y[i2], 
    F3 = Y[i3], F4 = Y[i4], F5 = Y[i5], 
    F6 = Y[i6], F7 = Y[i7], F8 = Y[i8];
  
  // Compute the finite differences, without scaling factors
  double Fu = F1 - F5;
  double Fv = F3 - F7;
  double Fuu = F1 + F5 - F0 - F0;
  double Fvv = F3 + F7 - F0 - F0;
  double Fuv = F2 + F6 - F4 - F8;

  // Compute the generalized laplacian using precomputed coefficients
  return Cuu * Fuu + Cuv * Fuv + Cvv * Fvv + Cu * Fu + Cv * Fv - rho;
}

double FDBorderSite::ComputeEquation(double *Y)
{
  // Compute the finite differences, without scaling factors
  double Fu = Y[ xNeighbor[1] ] - Y[ xNeighbor[3] ];
  double Fv = Y[ xNeighbor[2] ] - Y[ xNeighbor[4] ];

  // Compute the generalized gradient magnitude using precomputed values
  return CuCu * Fu * Fu + CuCv * Fu * Fv + CvCv * Fv * Fv - C0 * Y[ xNeighbor[0] ];
}

/** Compute the derivatives for Newton's method */
void FDInternalSite::ComputeDerivative(double *, double *A, unsigned int iRow)
{
  // Simply copy the values of the b's into A
  A[ xEntry[0] ] = b0;
  A[ xEntry[1] ] = b1;
  A[ xEntry[2] ] = b2;
  A[ xEntry[3] ] = b3;
  A[ xEntry[4] ] = b4;
  A[ xEntry[5] ] = b5;
  A[ xEntry[6] ] = b6;
  A[ xEntry[7] ] = b7;
  A[ xEntry[8] ] = b8;
}

/** Compute the derivatives for Newton's method */
void FDBorderSite::ComputeDerivative(double *Y, double *A, unsigned int iRow)
{
  // This computes the derivative of the site's equation with respect to the
  // four variables involved in the equation.
  double Fu = Y[ xNeighbor[1] ] - Y[ xNeighbor[3] ];
  double Fv = Y[ xNeighbor[2] ] - Y[ xNeighbor[4] ];

  double d1 = (CuCu + CuCu) * Fu + CuCv * Fv;
  double d2 = CuCv * Fu + (CvCv + CvCv) * Fv;

  // Clear all the entries to zero
  for(unsigned int j = 0; j < nDistinctSites; j++)
    A[j] = 0.0;

  // Compute the derivatives
  A[ xEntry[0] ] -= C0;
  A[ xEntry[1] ] += d1;
  A[ xEntry[2] ] += d2;
  A[ xEntry[3] ] -= d1;
  A[ xEntry[4] ] -= d2;
}

/** Initialize internal site */
void FDInternalSite::SetGeometry(GeometryDescriptor *g, double rho)
{
  // Store the rho
  this->rho = rho;

  // Compute the finite difference premultipliers
  double _du = m - 1;
  double _dv = n - 1;
  double _du2 = _du * _du;
  double _dv2 = _dv * _dv;
  double _duv = _du * _dv;
  
  // Here we need to compute the coefficients for the site
  Cuu = g->xContravariantTensor[0][0] * _du2;
  Cvv = g->xContravariantTensor[1][1] * _dv2;
  Cuv = 0.5 * g->xContravariantTensor[0][1] * _duv;
  Cu  = - 0.5 * _du * (
    g->xContravariantTensor[0][0] * g->xChristoffelSecond[0][0][0] + 
    2.0 * g->xContravariantTensor[0][1] * g->xChristoffelSecond[0][1][0] +
    g->xContravariantTensor[1][1] * g->xChristoffelSecond[1][1][0]);
  Cv = - 0.5 * _dv * (
    g->xContravariantTensor[0][0] * g->xChristoffelSecond[0][0][1] +
    2.0 * g->xContravariantTensor[0][1] * g->xChristoffelSecond[0][1][1] +
    g->xContravariantTensor[1][1] * g->xChristoffelSecond[1][1][1]);

  // Compute the derivatives with respect to each neighbor
  b0 = -2.0 * (Cuu + Cvv);
  b1 = Cuu + Cu;
  b2 = Cuv;
  b3 = Cvv + Cv;
  b4 = -Cuv;
  b5 = Cuu - Cu;
  b6 = Cuv;
  b7 = Cvv - Cv;
  b8 = -Cuv;
}

/** Initialize the border site */
void FDBorderSite::SetGeometry(GeometryDescriptor *g, double)
{
  // Compute the finite difference premultipliers (they depend on the orientation)
  double _du = wu * (m - 1);
  double _dv = wv * (n - 1);
  double _du2 = _du * _du;
  double _dv2 = _dv * _dv;
  double _duv = _du * _dv;
  
  // Here we need to compute the coefficients for the site
  CuCu = g->xContravariantTensor[0][0] * _du2;
  CuCv = 2.0 * g->xContravariantTensor[0][1] * _duv;
  CvCv = g->xContravariantTensor[1][1] * _dv2;
  C0 = 4.0;
}

void DumpQMatrix(QMatrix *q)
{
  char ch = '{';
  unsigned int n = Q_GetDim(q);
  cout << "SparseArray[";
  for(unsigned int r = 1; r <= n; r++)
    {
    unsigned int k = Q_GetLen(q, r);
    for(unsigned int p = 0; p < k; p++)
      {
      unsigned int c = Q_GetPos(q, r, p);
      double v = Q_GetVal(q, r, p);
      if(v < 0.000001 && v > -0.000001) v = 0;
      cout << ch << " {" << r << "," << c << "}->" << v;
      ch = ',';
      }
    }
  cout << "}]" << endl;
}

/*
void DumpVector(Vector *q)
{
  char c = '{';
  unsigned int n = V_GetDim(q);
  
  cout << "SparseArray[";
  for(unsigned int i = 1; i <= n; i++)
    {
    double v = V_GetCmp(q,i);
    if(v < 0.000001 && v > -0.000001) v = 0;
    cout << c << i << "->" << v;
    c = ',';
    }
  cout << "}]" << endl;
}

double MaxAbs_VV(Vector *q)
{
  double cand = 0;
  unsigned int n = V_GetDim(q);
  for(unsigned int i = 1; i <= n; i++)
    if(cand == 0 || abs(V_GetCmp(q,i)) > cand)
      cand = abs(V_GetCmp(q,i));
  return cand;
}
*/

MedialPDESolver
::MedialPDESolver(unsigned int m, unsigned int n)
{
  // Copy the site dimensions
  this->m = m; this->n = n;
  du = 1.0 / (m - 1); dv = 1.0 / (n - 1);

  nSites = m * n;
  xSites = new FDAbstractSite*[nSites];
  xSiteIndex = new unsigned int*[m];

  // Initialize the sites and the site index
  unsigned int iSite = 0, i, j;
  for(i = 0; i < m; i++) 
    {
    xSiteIndex[i] = new unsigned int[n];
    for(j = 0; j < n; j++) 
      {
      // Decide whether to create a border site or a center site
      if(i == 0 || j == 0 || i == m-1 || j == n-1)
        xSites[iSite] = new FDBorderSite(m, n, i, j, iSite, n, 1);
      else
        xSites[iSite] = new FDInternalSite(m, n, i, j, iSite, n, 1);

      // Associate the raw offset with the u-v index
      xSiteIndex[i][j] = iSite++;
      }
    }

  // Initialize the sparse matrix row index (use Fortran-based indexing)
  xRowIndex = new unsigned int[nSites + 1];
  xRowIndex[0] = 1;
  for(iSite = 0; iSite < nSites; iSite++)
    xRowIndex[iSite+1] = xRowIndex[iSite] + xSites[iSite]->GetNumberOfNetwonColumns();

  // Initialize the column index
  nSparseEntries = xRowIndex[nSites] - 1;
  xColIndex = new unsigned int[nSparseEntries];
  for(iSite = 0; iSite < nSites; iSite++)
    {
    unsigned int nEntries = xRowIndex[iSite+1] - xRowIndex[iSite];
    for(unsigned int iEntry = 0; iEntry < nEntries; iEntry++)
      xColIndex[ xRowIndex[iSite] + iEntry - 1 ] = xSites[iSite]->GetNewtonColumnAtIndex(iEntry) + 1;
    }

  // Initialize the sparse data
  xSparseValues = new double[nSparseEntries];

  // Compute the initial solution
  xInitSoln = new double[nSites];
  for(iSite = 0; iSite < nSites; iSite++)
    xInitSoln[iSite] = xSites[iSite]->GetInitialValue();

  // Initialize our three vectors
  b = new double[nSites];
  y = new double[nSites];
  eps = new double[nSites];
}

extern "C" {
  typedef unsigned int CBLAS_INDEX;
  void cblas_dcopy(const int N, const double *X, const int incX, double *Y, const int incY);
  void cblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY);
  CBLAS_INDEX cblas_idamax(const int N, const double *X, const int incX);
}

void
MedialPDESolver
::Solve(IMedialPDEProblem *problem, double delta)
{
  // The Jet of a surface
  double X[3], Xu[3], Xv[3], Xuu[3], Xuv[3], Xvv[3];
  double u = 0.0, v = 0.0;
  unsigned int i, j, k;

  // Initialize each site with the current surface properties
  for(i = 0; i < m; i++, u+=du, v=0.0)
    {
    for(j = 0; j < n; j++, v+=dv)
      {
      // Get the index of the site
      unsigned int iSite = xSiteIndex[i][j];

      // Compute the surface jet and the laplacian
      problem->ComputeJet2(u,v,X,Xu,Xv,Xuu,Xuv,Xvv);

      // Compute the geometrical properties
      GeometryDescriptor gd(X,Xu,Xv,Xuu,Xuv,Xvv);

      // Compute the solution at this point
      xSites[iSite]->SetGeometry(&gd, problem->ComputeLaplacian(u,v));
      }
    }

  // Copy the initial solution to the current solution
  cblas_dcopy(nSites, xInitSoln, 1, y, 1);

  // Initialize epsilon to zero
  memset(eps, 0, sizeof(double) * nSites);

  // We are now ready to perform the Newton loop
  bool flagComplete = false;
  for(unsigned int iIter = 0; iIter < 16 && !flagComplete; iIter++)
    {
    // Compute the A matrix and the b vector at each node
    for(k = 0; k < nSites; k++)
      {
      // Compute the non-zero values of A for this row
      xSites[k]->ComputeDerivative(y, xSparseValues + xRowIndex[k] - 1, k+1);

      // Compute the value of b
      b[k] = - xSites[k]->ComputeEquation(y);
      }

    // Solve the equation for X
    // BiCGIter(&A, &eps, &b, nSites, NULL, 1.0);

    /*
    cout << "A = "; DumpQMatrix(&A);
    cout << "b = "; DumpVector(&b);
    cout << "eps = "; DumpVector(&eps);
    */

    // Get the largest error (eps)
    double epsMax = eps[cblas_idamax(nSites, eps, 1)];
    double bMax = b[cblas_idamax(nSites, b, 1)];

    // Append the epsilon vector to the result
    cblas_daxpy(nSites, 1.0, eps, 1, y, 1);
    
    // Print the statistics
    cout << "-----------" << endl;
    cout << "Step " << iIter << ": " << endl;
    cout << "  Largest Epsilon: " << epsMax << endl;
    cout << "  Largest Eqn Error: " << bMax << endl;

    // Convergence is defined when epsilon is smaller than some threshold
    flagComplete = (epsMax < delta);
    }
}
