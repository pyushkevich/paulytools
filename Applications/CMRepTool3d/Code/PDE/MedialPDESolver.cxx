#include "MedialPDESolver.h"
#include <cmath>
#include <algorithm>
#include "smlmath.h"

// BLAS/PARDISO references
extern "C" {
  typedef unsigned int CBLAS_INDEX;
  void cblas_dcopy(const int N, const double *X, const int incX, double *Y, const int incY);
  void cblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY);
  CBLAS_INDEX cblas_idamax(const int N, const double *X, const int incX);

  void pardisoinit_(int *, int *, int *);
  void pardiso_(int *, int *, int *, int *, int *, int *, double *, int *, int *, 
    int *, int *, int *, int *, double *, double *, int*);
}

using namespace std;

void
GeometryDescriptor
::SetJet(double *X, double *Xu, double *Xv, double *Xuu, double *Xuv, double *Xvv)
{
  // Compute the covariant tensor
  xCovariantTensor[0][0] = Xu[0] * Xu[0] + Xu[1] * Xu[1] + Xu[2] * Xu[2];
  xCovariantTensor[1][1] = Xv[0] * Xv[0] + Xv[1] * Xv[1] + Xv[2] * Xv[2];
  xCovariantTensor[1][0] = Xu[0] * Xv[0] + Xu[1] * Xv[1] + Xu[2] * Xv[2];
  xCovariantTensor[0][1] = xCovariantTensor[1][0];

  // Compute the determinant of the covariant tensor
  g = xCovariantTensor[0][0] * xCovariantTensor[1][1] 
    - xCovariantTensor[0][1] * xCovariantTensor[0][1];
  gInv = 1.0 / g;

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

void
MedialAtom
::ComputeDifferentialGeometry()
{
  G.SetJet( X.data_block(), Xu.data_block(), Xv.data_block(),
    Xuu.data_block(), Xuv.data_block(), Xvv.data_block());
}

void MedialAtom ::ComputeNormalVector()
{
  N = cross_3d(Xu, Xv) * sqrt(G.gInv);
}

bool MedialAtom::ComputeBoundaryAtoms()
{
  // Terms going into gradR
  double CXu = 
    Ru * G.xContravariantTensor[0][0] + Rv * G.xContravariantTensor[0][1]; 
  double CXv = 
    Ru * G.xContravariantTensor[0][1] + Rv * G.xContravariantTensor[1][1];
    
  // Compute Grad R
  xGradR[0] = Xv[0] * CXv + Xu[0] * CXu;
  xGradR[1] = Xv[1] * CXv + Xu[1] * CXu;
  xGradR[2] = Xv[2] * CXv + Xu[2] * CXu;

  // Compute squared length of GradR
  double xMagGradR2 = Ru * CXu + Rv * CXv;
  double sinTheta2 = 1.0f - xMagGradR2;
  
  // Correct the floating point / solver error
  if(fabs(sinTheta2) < 1e-8) sinTheta2 = 0.0;
  
  // Compute the boundary sites
  if (sinTheta2 > 0.0)
    {
    // We are not at a crest, valid atom
    flagCrest = false; flagValid = true;
    
    // Compute the normal and tangent components of the boundaries
    SMLVec3d CN = N * sqrt(sinTheta2);

    // Compute the position and normals of boundary points
    xBnd[0].N = - xGradR - CN;
    xBnd[1].N = - xGradR + CN;
    xBnd[0].X = X + xBnd[0].N * R;
    xBnd[1].X = X + xBnd[1].N * R;
    }
  else if (sinTheta2 < 0.0)
    { 
    flagValid = false; 
    // cerr << "Bad atom ("<< R << ", " << sinTheta2 << ") at " << u << ", " << v << endl;
    cout << "x" << flush;
    }
  else
    { 
    // We are at a crest, valid atom
    flagValid = true; flagCrest = true;

    // Simpler geometry, save a square root!
    xBnd[0].N = xBnd[1].N = -xGradR;
    xBnd[0].X = xBnd[1].X = X - xGradR * R;
    }

  return flagValid;
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
  _du = (double) (m - 1);
  _dv = (double) (n - 1);
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
    xColumns[0] = i6; xEntry[6] = 0;
    xColumns[1] = i5; xEntry[5] = 1;
    xColumns[2] = i4; xEntry[4] = 2;
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
  nDistinctSites = 5;

  // Now, depending on which border we are on, correct the index and weights of finite
  // differences in u and v directions
  if(i == 0)
    { xNeighbor[3] = iNode; wu = 1; nDistinctSites--; }
  else if(i == m - 1)
    { xNeighbor[1] = iNode; wu = 1; nDistinctSites--; }
  if(j == 0)
    { xNeighbor[4] = iNode; wv = 1; nDistinctSites--; }
  else if(j == n - 1)
    { xNeighbor[2] = iNode; wv = 1; nDistinctSites--; }

  // Check whether we are at a corner
  corner = (nDistinctSites == 3);

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

double FDInternalSite::ComputeJacobiIteration(double *Y)
{
   // Get the components
  double 
    F0 = Y[i0], F1 = Y[i1], F2 = Y[i2], 
    F3 = Y[i3], F4 = Y[i4], F5 = Y[i5], 
    F6 = Y[i6], F7 = Y[i7], F8 = Y[i8];
  
  // Compute the finite differences, without scaling factors
  double Fu = F1 - F5;
  double Fv = F3 - F7;
  double Fuu = F1 + F5;
  double Fvv = F3 + F7;
  double Fuv = F2 + F6 - F4 - F8;

  // Compute the generalized laplacian using precomputed coefficients
  
  return (Cuu * Fuu + Cuv * Fuv + Cvv * Fvv + Cu * Fu + Cv * Fv - rho) 
    / (2.0 * (Cuu + Cvv));
}

double FDBorderSite::ComputeEquation(double *Y)
{
  // Compute the finite differences, without scaling factors
  double Fu = Y[ xNeighbor[1] ] - Y[ xNeighbor[3] ];
  double Fv = Y[ xNeighbor[2] ] - Y[ xNeighbor[4] ];

  // Compute the generalized gradient magnitude using precomputed values
  return CuCu * Fu * Fu + CuCv * Fu * Fv + CvCv * Fv * Fv - C0 * Y[ xNeighbor[0] ];
}

double FDBorderSite::ComputeJacobiIteration(double *Y)
{
  // An array of weights for all combinations of parameters
  double A = CuCu, C = CvCv, B = 0.5 * CuCv;
  double Z[4][4] = {
    {  A,  B, -A, -B }, 
    {  B,  C, -B, -C },
    { -A, -B,  A,  B },
    { -B, -C,  B,  C }}; 

  // Create accumulators to hold coefficients of X0^2, X0 and const
  double a = 0.0, b = -C0, c = 0.0;

  // Loop over all the products
  for(unsigned int i = 0; i < 4; i++)
    for(unsigned int j = 0; j < 4; j++)
      {
      bool ti = (xNeighbor[i+1] == xNeighbor[0]);
      bool tj = (xNeighbor[j+1] == xNeighbor[0]);
      if(ti && tj)
        a += Z[i][j];
      else if(ti)
        b += Z[i][j] * Y[xNeighbor[j+1]];
      else if(tj)
        b += Z[i][j] * Y[xNeighbor[i+1]];
      else 
        c += Z[i][j] * Y[xNeighbor[i+1]] * Y[xNeighbor[j+1]];
      }

  // Solve the quadratic equation ax^2 + bx + c = 0
  double D = b * b - 4 * a * c;
  if(D < 0)
    return 0; // Solution does not exist

  // Get the two roots of the quadratic equation
  double x = ( - b - sqrt(D) ) * 0.5 / a; 

  // Select the root that results in the outward gradient
  return x;
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

// Compute the R, Ru and Rv given the solution - used to compute atoms
void FDInternalSite::ComputeRJet(double *Y, double &R, double &Ru, double &Rv)
{
  // Compute the derivatives in the U and V directions
  double Fu = 0.5 * _du * (Y[i1] - Y[i5]);
  double Fv = 0.5 * _dv * (Y[i3] - Y[i7]);
  
  // Convert to the derivatives of R
  R = sqrt( Y[i0] );
  Ru = Fu / (2.0 * R);
  Rv = Fv / (2.0 * R);     
}

// Compute the R, Ru and Rv given the solution - used to compute atoms
void FDBorderSite::ComputeRJet(double *Y, double &R, double &Ru, double &Rv)
{
  // Compute the derivatives in the U and V directions
  double Fu = _du * (Y[ xNeighbor[1] ] - Y[ xNeighbor[3] ]);
  double Fv = _dv * (Y[ xNeighbor[2] ] - Y[ xNeighbor[4] ]);
  
  // Convert to the derivatives of R
  R = sqrt( Y[xNeighbor[0]] );
  Ru = Fu / (2.0 * R);
  Rv = Fv / (2.0 * R);
}


/** Initialize internal site */
void FDInternalSite::SetGeometry(GeometryDescriptor *g, double rho)
{
  // Store the rho
  this->rho = rho;

  // Compute the finite difference premultipliers
  _du = m - 1;
  _dv = n - 1;
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
  _du = wu * (m - 1);
  _dv = wv * (n - 1);
  double _du2 = _du * _du;
  double _dv2 = _dv * _dv;
  double _duv = _du * _dv;
  
  // Here we need to compute the coefficients for the site
  CuCu = g->xContravariantTensor[0][0] * _du2;
  CuCv = 2.0 * g->xContravariantTensor[0][1] * _duv;
  CvCv = g->xContravariantTensor[1][1] * _dv2;
  C0 = 4.0;
}

void SparseMultiply(unsigned int n, unsigned int *rowIndex, unsigned int *colIndex, 
	double *values, double *x, double *y)
{
	for(unsigned int r = 1; r <=n ; r++)
    	{
        y[r-1] = 0;
    	for(unsigned int k = rowIndex[r-1]; k < rowIndex[r]; k++) 
      		{
      		// r and c are 1-based
      		int c = colIndex[k - 1] ;	
      		double v = values[k - 1];
      		y[r-1] += v * x[c - 1];
      		}
    	}	
}

void SparseLinearTest(unsigned int n, unsigned int *rowIndex, unsigned int *colIndex, 
	double *values, double *x, double *y,double *b)
{
	SparseMultiply(n, rowIndex, colIndex, values, x, y);
	cblas_daxpy(n, -1.0, b, 1, y, 1);
}

void DumpSparseMatrix(unsigned int n, unsigned int *rowIndex, unsigned int *colIndex, double *values)
{
  char ch = '{';
  cout << "SparseArray[";
  for(unsigned int r = 1; r <=n ; r++)
    {
    for(unsigned int c = rowIndex[r-1]; c < rowIndex[r]; c++) 
      {
      double v = values[c - 1];
      if(v < 0.000001 && v > -0.000001) v = 0;
      cout << ch << " {" << r << "," << colIndex[c-1] << "}->" << v;
      ch = ',';
      }
      cout << endl;
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

  // Initialize our three vectors
  b = new double[nSites];
  y = new double[nSites];
  eps = new double[nSites];
  zTest = new double[nSites];
  
  // Construct the evaluation grids
  double u = 0.0, v = 0.0;
  xGridU = new double[m]; xGridV = new double[n];
  for(size_t i = 0; i < m-1; i++, u+=du) xGridU[i] = u;
  for(size_t j = 0; j < n-1; j++, v+=dv) xGridV[j] = v;
  xGridU[m-1] = 1.0; xGridV[n-1] = 1.0;

  // Initialize the medial atom array
  xAtoms = new MedialAtom[nSites];

  // Initialize the PARDISO solver
  MTYPE = 11; // Nonsymmetric real
  memset(IPARM, 0, sizeof(int) * 64);
  pardisoinit_(PT,&MTYPE,IPARM);
  IPARM[2] = 1;

  // Compute the initial solution
  xInitSoln = new double[nSites];
  SetDefaultInitialGuess(0.001);

  // Initialize the atom grid to a Cartesian grid
  xGrid = new CartesianMedialAtomGrid(m, n);
}

void
MedialPDESolver
::SetDefaultInitialGuess(double xMagnitude)
{
  // Compute the initial solution
  for(unsigned int iSite = 0; iSite < nSites; iSite++)
    xInitSoln[iSite] = xMagnitude * xSites[iSite]->GetInitialValue();
}

void 
MedialPDESolver
::SetSolutionAsInitialGuess()
{
  for(unsigned int iSite = 0; iSite < nSites; iSite++)
    xInitSoln[iSite] = y[iSite];
}

void MedialPDESolver
::TestJacobi()
{
  // Compute the Jacobi iteration
  for(unsigned int k = 0; k < nSites; k++)
    {
    b[k] = xSites[k]->ComputeJacobiIteration(y);
    if(b[k] <= 0.0)
      {
      cout << "Bad Jacobi result at site " << k << endl;
      }
    b[k] -= y[k];
    }

  // Compute the max update value
  double bMax = fabs(b[cblas_idamax(nSites, b, 1)]);
  cout << "Jacobi Update: Max Differenece is " << bMax << endl;
}

void 
MedialPDESolver
::SetMedialSurface(IHyperSurface2D *xSurface)
{
  // Remember the surface
  this->xSurface = xSurface;

  // Tell the surface where it will be evaluated
  xSurface->SetEvaluationGrid(m, n, xGridU, xGridV);
}

void MedialPDESolver::InitializeSiteGeometry()
{
  // Initialize each site with the current surface properties
  for(size_t i = 0; i < m; i++) for(size_t j = 0; j < n; j++)
    {
    // Get the index of the site
    unsigned int iGrid = xGrid->GetAtomIndex(i, j);
    unsigned int iSite = xSiteIndex[i][j];

    // Access the medial atom underneath
    MedialAtom &xAtom = xAtoms[iGrid];

    // Set the atoms' domain coordinates
    xAtom.u = xGridU[i]; xAtom.v = xGridV[j];

    // Compute the surface jet and the laplacian
    xSurface->EvaluateAtGridIndex(i, j, 0, 0, xAtom.X.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 1, 0, xAtom.Xu.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 0, 1, xAtom.Xv.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 2, 0, xAtom.Xuu.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 1, 1, xAtom.Xuv.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 0, 2, xAtom.Xvv.data_block());

    // Compute the differential geometric tensors
    xAtom.ComputeDifferentialGeometry();

    // Compute the normal vector
    xAtom.ComputeNormalVector();

    // Compute the laplacian of R (TODO!)
    xAtom.xLapR = -0.25;

    // Compute the solution at this point
    xSites[iSite]->SetGeometry( &xAtom.G, xAtom.xLapR);
    }
}

void
MedialPDESolver
::ReconstructAtoms(double *ySolution)
{
  // Once the iterations are complete, reconstruct the atoms
  for(unsigned int i = 0; i < m; i++) for(unsigned int j = 0; j < n; j++)
    {
    // Map to a grid index and to our own index
    unsigned int iGrid = xGrid->GetAtomIndex(i, j);
    unsigned int iLocal = xSiteIndex[i][j];
    
    // The medial atom to update
    MedialAtom &xAtom = xAtoms[iGrid];

    // The case where phi is negative is undefined
    if(ySolution[iLocal] > 0)
      {
      // Compute the derivatives of R using finite differences
      xSites[iLocal]->ComputeRJet(ySolution, xAtom.R, xAtom.Ru, xAtom.Rv);

      // Compute the boundary properties of the medial point
      xAtom.ComputeBoundaryAtoms();
      }
    else
      {
      // What are we supposed to do?
      xAtom.R = xAtom.Ru = xAtom.Rv = 0.0;
      xAtom.flagValid = false;
      }
    }
}

void
MedialPDESolver
::Solve(double delta)
{
  unsigned int i, j, k;

  // Intialize the sites
  InitializeSiteGeometry();

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
      b[k] = -xSites[k]->ComputeEquation(y);
      }

    // Debug 
    // DumpSparseMatrix(nSites, xRowIndex, xColIndex, xSparseValues);

    // Solve the equation for X
    int MAXFCT = 1, MNUM = 1, PHASE = 11, N = nSites, NRHS = 1, MSGLVL = 0, ERROR = 0; 
    if(iIter == 0) 
      {
      pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, 
        xSparseValues, (int *) xRowIndex, (int *) xColIndex, 
        NULL, &NRHS, IPARM, &MSGLVL, 
        b, eps, &ERROR);
      }

    // Only perform the solution step
    PHASE=23;
    pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, 
      xSparseValues, (int *) xRowIndex, (int *) xColIndex, 
      NULL, &NRHS, IPARM, &MSGLVL, 
      b, eps, &ERROR);


    // Test the matrix result
    // SparseLinearTest(nSites, xRowIndex, xColIndex, xSparseValues, eps, zTest, b);

    // Get the largest error (eps)
    double epsMax = fabs(eps[cblas_idamax(nSites, eps, 1)]);
    double bMax = fabs(b[cblas_idamax(nSites, b, 1)]);
    // double zMax = fabs(zTest[cblas_idamax(nSites, zTest, 1)]);

    // Append the epsilon vector to the result
    cblas_daxpy(nSites, 1.0, eps, 1, y, 1);

    // Print the statistics
    
    // cout << "-----------" << endl;
    // cout << "Step " << iIter << ": " << endl;
    // cout << "  Largest Epsilon: " << epsMax << endl;
    // cout << "  Largest Eqn Error: " << bMax << endl;
    // cout << "  Largest Solver Error: " << zMax << endl;
   

    // Convergence is defined when epsilon is smaller than some threshold
    flagComplete = (epsMax < delta);
    }

  // Check for problems
  /* 
  MedialAtom *xAtoms = GetAtomArray();
  double gMin = 1.0e10;
  for(size_t i = 0; i < GetNumberOfAtoms(); i++)
    {
    if(xAtoms[i].G.g < gMin)
      gMin = xAtoms[i].G.g; 
    }
  cout << "Lowest g = " << gMin << endl;
  */
  
  // Reconstruct the medial atoms
  ReconstructAtoms(y);
}

void 
MedialPDESolver
::SolveByJacobiMethod(double delta)
{
  // Copy the initial solution to the current solution
  cblas_dcopy(nSites, xInitSoln, 1, y, 1);

  // Initialize the geometry
  InitializeSiteGeometry();
  
  // We are now ready to perform the Newton loop
  bool flagComplete = false;
  for(unsigned int iIter = 0; iIter < 30000 && !flagComplete; iIter++)
    {
    // Run the Jacobi iteration, put the result in b, difference in eps
    for(unsigned int k = 0; k < nSites; k++)
      {
      b[k] = xSites[k]->ComputeJacobiIteration(y);
      eps[k] = b[k] - y[k];
      zTest[k] = xSites[k]->ComputeEquation(y);
      if(b[k] <= 0) 
        cout << "Bad result " << b[k] << " at site " << k << endl;
      
      // Optional: Gauss-Seidell
      y[k] = b[k];
      }

    // Copy the result into y
    cblas_dcopy(nSites, b, 1, y, 1);

    // Check the maximum update difference - should be 0 if we have the true solution
    double epsMax = fabs(eps[cblas_idamax(nSites, eps, 1)]);
    double zMax = fabs(eps[cblas_idamax(nSites, zTest, 1)]);
    if(epsMax < delta)
      flagComplete = true;

    // Report the iteration progress
    cout << "-----------" << endl;
    cout << "Step " << iIter << ": " << endl;
    cout << "  Largest Epsilon: " << epsMax << endl;
    cout << "  Largest Error:   " << zMax << endl;
    }

  // Compute the atoms
  ReconstructAtoms(y);
}

double MedialPDESolver
::IntegrateBoundaryMeasure(EuclideanFunction *fun, double &area)
{
  // Integration over patches. First, compute the boundary measure at
  // each point on the medial surface.
  double *bndValues[2], *xAreas;
  bndValues[0]= new double[nSites];
  bndValues[1]= new double[nSites];

  // This vector stores the area assigned to each boundary point
  xAreas = new double[nSites];

  // Compute the boundary measure at each point
  for(unsigned int iSite = 0; iSite < nSites; iSite++)
    {
    // For boundary sites, only one side is used
    if(xSites[iSite]->IsBorderSite())
      {
      // Compute the crest measure
      bndValues[0][iSite] = bndValues[1][iSite] = 
        fun->Evaluate(xAtoms[iSite].xBnd[0].X);
      }
    else
      {
      // Compute the boundary measure
      bndValues[0][iSite] = fun->Evaluate(xAtoms[iSite].xBnd[0].X);
      bndValues[1][iSite] = fun->Evaluate(xAtoms[iSite].xBnd[1].X);
      }
    }

  // Integrate the surface area
  double xArea = 0.0, xMeasure = 0.0;
  for(unsigned int i = 0; i < m - 1; i++)
    for(unsigned int j = 0; j < n - 1; j++)
      for(unsigned int d = 0; d < 2; d++)
        {
        // Access the four medial points
        unsigned int i00 = xSiteIndex[i][j];
        unsigned int i01 = xSiteIndex[i][j+1];
        unsigned int i10 = xSiteIndex[i+1][j];
        unsigned int i11 = xSiteIndex[i+1][j+1];

        // Compute the areas of the two triangles involved
        double A1 = TriangleArea(
          xAtoms[i00].xBnd[d].X, xAtoms[i01].xBnd[d].X, xAtoms[i11].xBnd[d].X); 
        double A2 = TriangleArea(
          xAtoms[i00].xBnd[d].X, xAtoms[i11].xBnd[d].X, xAtoms[i10].xBnd[d].X); 
          
        // Integrate the measure
        xMeasure += A1 * (bndValues[d][i00] + bndValues[d][i01] + bndValues[d][i11]); 
        xMeasure += A2 * (bndValues[d][i00] + bndValues[d][i11] + bndValues[d][i10]); 
        xArea += (A1 + A2);
        }

  // Scale the measure by three
  xMeasure /= 3.0;
  // cout << "Area: " << xArea << " Measure: " << xMeasure << endl;

  // Clean up
  delete xAreas;
  delete bndValues[0];
  delete bndValues[1];

  // Return the result
  area = xArea;
  return xMeasure;
}

double MedialPDESolver
::IntegrateVolumeMeasure(EuclideanFunction *xFunction, unsigned int nSamples, double &volume)
{
  // Integrate some measure over the volume enclosed by the structure. A
  // triangle on the medial surface and a triangle on the boundary form a 
  // prism; we have to calculate its volume.
  unsigned int iSite, iSample, iSide, nSides, i, j;
  
  // First we sample the function along the profiles from the medial axis and
  // store this for future integration
  double **xSamples[2];
  xSamples[0] = new double *[nSites];
  xSamples[1] = new double *[nSites];

  for(iSite = 0; iSite < nSites; iSite++)
    {
    // Allocate the along-the-profile samples
    xSamples[0][iSite] = new double[nSamples];
    xSamples[1][iSite] = xSites[iSite]->IsBorderSite() 
      ? xSamples[0][iSite] : new double[nSamples];

    // Simpler processing at border sites
    nSides = xSites[iSite]->IsBorderSite() ? 1 : 2;
    for(iSide = 0; iSide < nSides; iSide++)
      {
      // Compute the delta vector
      SMLVec3d xSample = xAtoms[iSite].X;
      SMLVec3d xDelta = 
        (xAtoms[iSite].R / (nSamples - 1.0)) * xAtoms[iSite].xBnd[iSide].N;

      // Compute the delta vector
      for(iSample = 0; iSample < nSamples; iSample++, xSample += xDelta)
        {
        xSamples[iSide][iSite][iSample] = xFunction->Evaluate(xSample);
        }
      }
    }

  // Now, the more difficult part. For each section of the boundary / medial
  // grid, we must compute the volume and assign a sixth of this volume to 
  // each of the samples, thus generating a weight for every sample
  double xMeasure = 0.0, xVolume = 0.0;
  for(i = 0; i < m - 1; i++) for(j = 0; j < n - 1; j++)
    {
    // Access the four sites at this quad
    unsigned int i00 = xSiteIndex[i][j];
    unsigned int i01 = xSiteIndex[i][j+1];
    unsigned int i10 = xSiteIndex[i+1][j];
    unsigned int i11 = xSiteIndex[i+1][j+1];

    // Access the four medial atoms
    SMLVec3d X00 = xAtoms[i00].X;
    SMLVec3d X01 = xAtoms[i01].X;
    SMLVec3d X10 = xAtoms[i10].X;
    SMLVec3d X11 = xAtoms[i11].X;

    // Get the medial atom at this location
    for(iSide = 0; iSide < 2; iSide++)
      {
      // Get the four deltas
      SMLVec3d d00 = xAtoms[i00].xBnd[iSide].N;
      SMLVec3d d01 = xAtoms[i01].xBnd[iSide].N;
      SMLVec3d d10 = xAtoms[i10].xBnd[iSide].N;
      SMLVec3d d11 = xAtoms[i11].xBnd[iSide].N;

      // Compute the prism volumes
      for(iSample = 0; iSample < nSamples - 1; iSample++)
        {
        // Compute the volume of the first prism
        double t1 = iSample * 1.0f / (nSamples - 1); 
        double t2 = (1 + iSample) * 1.0f / (nSamples - 1); 

        double V1 = PrismVolume(
          X00 + d00 * t1, X01 + d01 * t1, X11 + d11 * t1,
          X00 + d00 * t2, X01 + d01 * t2, X11 + d11 * t2);
          
        // Compute the volume of the second prism
        double V2 = PrismVolume(
          X00 + d00 * t1, X11 + d11 * t1, X10 + d10 * t1,
          X00 + d00 * t2, X11 + d11 * t2, X10 + d10 * t2);

        // Integrate the measure
        xMeasure += V1 * (
          xSamples[iSide][i00][iSample] + xSamples[iSide][i00][iSample + 1] + 
          xSamples[iSide][i01][iSample] + xSamples[iSide][i01][iSample + 1] + 
          xSamples[iSide][i11][iSample] + xSamples[iSide][i11][iSample + 1] ); 

        xMeasure += V2 * (
          xSamples[iSide][i00][iSample] + xSamples[iSide][i00][iSample + 1] + 
          xSamples[iSide][i11][iSample] + xSamples[iSide][i11][iSample + 1] + 
          xSamples[iSide][i10][iSample] + xSamples[iSide][i10][iSample + 1] ); 

        xVolume += (V1 + V2);
        }
      }
    }

  // Clean up the samples
  for(iSite = 0; iSite < nSites; iSite++)
    {
    delete xSamples[0][iSite];
    if(!xSites[iSite]->IsBorderSite())
      delete xSamples[1][iSite];
    }
  delete xSamples[0];
  delete xSamples[1];

  // Finish
  volume = xVolume; 
  return xMeasure / 6.0;
}
