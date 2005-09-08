#include "MedialPDESolver.h"
#include <cmath>
#include <algorithm>
#include "smlmath.h"

// BLAS/PARDISO references
extern "C" {
  void pardisoinit_(size_t *, int *, int *);
  void pardiso_(size_t *, int *, int *, int *, int *, int *, double *, int *, int *, 
    int *, int *, int *, int *, double *, double *, int*);
}

using namespace std;

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
  // Compute the squares of the finite differences
  _du2 = _du * _du;
  _dv2 = _dv * _dv;
  _duv = _du * _dv;

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

  // Compute the finite difference premultipliers (they depend on the orientation)
  _du = wu * (m - 1);
  _dv = wv * (n - 1);
  _du2 = _du * _du;
  _dv2 = _dv * _dv;
  _duv = _du * _dv;

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

void FDInternalSite::ComputePhiJet(double *Y, double &F, double &Fu, double &Fv)
{
  // Compute the derivatives in the U and V directions
  F = Y[i0];
  Fu = 0.5 * _du * (Y[i1] - Y[i5]);
  Fv = 0.5 * _dv * (Y[i3] - Y[i7]);
}

// Compute the R, Ru and Rv given the solution - used to compute atoms
void FDBorderSite::ComputePhiJet(double *Y, double &F, double &Fu, double &Fv)
{
  // Compute the derivatives in the U and V directions
  F = Y[xNeighbor[0]]; 
  Fu = _du * (Y[ xNeighbor[1] ] - Y[ xNeighbor[3] ]);
  Fv = _dv * (Y[ xNeighbor[2] ] - Y[ xNeighbor[4] ]);
}

/** Initialize internal site */
void FDInternalSite::SetGeometry(GeometryDescriptor *g, double rho)
{
  // Store the rho
  this->rho = rho;

  // Compute the finite difference premultipliers
  
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

/** 
 * Compute the variational derivative with respect to N
 */
void
FDInternalSite::
ComputeVariationalDerivativeX(double *Y, double *A, double *b, 
  MedialAtom *xAtom, MedialAtom *dAtom)
{
  // Simply copy the values of the b's into A
  A[ xEntry[0] ] = b0; A[ xEntry[1] ] = b1; A[ xEntry[2] ] = b2;
  A[ xEntry[3] ] = b3; A[ xEntry[4] ] = b4; A[ xEntry[5] ] = b5;
  A[ xEntry[6] ] = b6; A[ xEntry[7] ] = b7; A[ xEntry[8] ] = b8;

  // The rest of this routine computes the right hand side

  // First get shorthands for the relevant vectors
  SMLVec3d &Nu = dAtom->Xu; SMLVec3d &Xu = xAtom->Xu;
  SMLVec3d &Nv = dAtom->Xv; SMLVec3d &Xv = xAtom->Xv;
  SMLVec3d &Nuu = dAtom->Xuu; SMLVec3d &Xuu = xAtom->Xuu;
  SMLVec3d &Nuv = dAtom->Xuv; SMLVec3d &Xuv = xAtom->Xuv;
  SMLVec3d &Nvv = dAtom->Xvv; SMLVec3d &Xvv = xAtom->Xvv;

  // Get shorthand for the differential geometric operators
  double (*G1)[2] = xAtom->G.xCovariantTensor;
  double (*G2)[2] = xAtom->G.xContravariantTensor;
  double (*K)[2][2] = xAtom->G.xChristoffelFirst;
  double (*K2)[2][2] = xAtom->G.xChristoffelSecond;
  double &g = xAtom->G.g;
  
  // Compute the derivatives of the contravariant tensor
  double NuXu = dot_product(Nu, Xu);
  double NvXu = dot_product(Nv, Xu);
  double NuXv = dot_product(Nu, Xv);
  double NvXv = dot_product(Nv, Xv);

  // The derivative of 'g'
  double dgdN = dAtom->G.g = 
    2 * ( NuXu * G1[1][1] + NvXv * G1[0][0] - (NuXv + NvXu) * G1[0][1] );
  double g2 = g * g;

  // The derivative of the contravariant tensor
  double (*z)[2] = dAtom->G.xContravariantTensor; 
  z[0][0] = (2 * NvXv * g - G1[1][1] * dgdN) / g2;
  z[1][1] = (2 * NuXu * g - G1[0][0] * dgdN) / g2;
  z[0][1] = z[1][0] = - ((NuXv + NvXu) * g - G1[1][0] * dgdN) / g2;

  // Compute the derivative of the Christoffel symbols of the second kind
  double NuXuu = dot_product(Nu, Xuu), NvXuu = dot_product(Nv, Xuu);
  double NuXuv = dot_product(Nu, Xuv), NvXuv = dot_product(Nv, Xuv);
  double NuXvv = dot_product(Nu, Xvv), NvXvv = dot_product(Nv, Xvv);
  double NuuXu = dot_product(Nuu, Xu), NuuXv = dot_product(Nuu, Xv);
  double NuvXu = dot_product(Nuv, Xu), NuvXv = dot_product(Nuv, Xv);
  double NvvXu = dot_product(Nvv, Xu), NvvXv = dot_product(Nvv, Xv);
  
  // The derivative of Christofel symbols of second kind  
  double (*Q)[2][2] = dAtom->G.xChristoffelSecond;
  Q[0][0][0] = 
    z[0][0] * K[0][0][0] + G2[0][0] * ( NuXuu + NuuXu ) +
    z[0][1] * K[0][0][1] + G2[0][1] * ( NvXuu + NuuXv );
  Q[0][0][1] = 
    z[1][0] * K[0][0][0] + G2[1][0] * ( NuXuu + NuuXu ) +
    z[1][1] * K[0][0][1] + G2[1][1] * ( NvXuu + NuuXv );
  
  Q[0][1][0] = Q[1][0][0] = 
    z[0][0] * K[1][0][0] + G2[0][0] * ( NuXuv + NuvXu ) +
    z[0][1] * K[1][0][1] + G2[0][1] * ( NvXuv + NuvXv );
  Q[0][1][1] = Q[1][0][1] = 
    z[1][0] * K[1][0][0] + G2[1][0] * ( NuXuv + NuvXu ) +
    z[1][1] * K[1][0][1] + G2[1][1] * ( NvXuv + NuvXv );

  Q[1][1][0] = 
    z[0][0] * K[1][1][0] + G2[0][0] * ( NuXvv + NvvXu ) +
    z[0][1] * K[1][1][1] + G2[0][1] * ( NvXvv + NvvXv );
  Q[1][1][1] = 
    z[1][0] * K[1][1][0] + G2[1][0] * ( NuXvv + NvvXu ) +
    z[1][1] * K[1][1][1] + G2[1][1] * ( NvXvv + NvvXv );

  // Compute the partials of phi
  double Fu = 0.5 * _du * ( Y[i1] - Y[i5] );
  double Fv = 0.5 * _dv * ( Y[i3] - Y[i7] );
  double Fuu = _du2 * ( Y[i1] + Y[i5] - Y[i0] - Y[i0] );
  double Fuv = 0.25 * _duv * ( Y[i2] + Y[i6] - Y[i4] - Y[i8] );
  double Fvv = _dv2 * ( Y[i3] + Y[i7] - Y[i0] - Y[i0] );

  // Compute the derivative (this should be the derivative)
  *b = -( 
    z[0][0] * (Fuu - K2[0][0][0] * Fu - K2[0][0][1] * Fv) - 
    G2[0][0] * (Q[0][0][0] * Fu + Q[0][0][1] * Fv) +
    z[1][1] * (Fvv - K2[1][1][0] * Fu - K2[1][1][1] * Fv) - 
    G2[1][1] * (Q[1][1][0] * Fu + Q[1][1][1] * Fv) +
    2.0 * ( 
      z[0][1] * (Fuv - K2[0][1][0] * Fu - K2[0][1][1] * Fv) - 
      G2[0][1] * (Q[0][1][0] * Fu + Q[0][1][1] * Fv) ));
}

void
FDInternalSite::
ComputeVariationalDerivativeRho(double *Y, double *A, double *b, 
  MedialAtom *xAtom, MedialAtom *dAtom)
{
  // Simply copy the values of the b's into A
  A[ xEntry[0] ] = b0; A[ xEntry[1] ] = b1; A[ xEntry[2] ] = b2;
  A[ xEntry[3] ] = b3; A[ xEntry[4] ] = b4; A[ xEntry[5] ] = b5;
  A[ xEntry[6] ] = b6; A[ xEntry[7] ] = b7; A[ xEntry[8] ] = b8;

  // The right hand side is just the derivative of rho
  *b = dAtom->xLapR;
}

void
FDBorderSite::ComputeVariationalDerivativeX(double *Y, double *A, double *b, 
  MedialAtom *xAtom, MedialAtom *dAtom)
{
  // First, let's compute grad phi
  double Fu = _du * (Y[ xNeighbor[1] ] - Y[ xNeighbor[3] ]);
  double Fv = _dv * (Y[ xNeighbor[2] ] - Y[ xNeighbor[4] ]);

  // Next, compute the weights for Hu and Hv and H
  double (*G1)[2] = xAtom->G.xCovariantTensor;
  double (*G2)[2] = xAtom->G.xContravariantTensor;
  double whu = 2.0 * _du * (G2[0][0] * Fu + G2[0][1] * Fv);
  double whv = 2.0 * _dv * (G2[1][0] * Fu + G2[1][1] * Fv);
  double wh = -4.0;

  // Clear all the entries to zero
  for(unsigned int j = 0; j < nDistinctSites; j++)
    A[j] = 0.0;
  A[ xEntry[0] ] += wh;
  A[ xEntry[1] ] += whu;
  A[ xEntry[2] ] += whv;
  A[ xEntry[3] ] -= whu;
  A[ xEntry[4] ] -= whv;

  // Compute the derivatives of the contravariant tensor
  SMLVec3d &Nu = dAtom->Xu; SMLVec3d &Nv = dAtom->Xv;
  SMLVec3d &Xu = xAtom->Xu; SMLVec3d &Xv = xAtom->Xv;

  double NuXu = dot_product(Nu, Xu);
  double NvXu = dot_product(Nv, Xu);
  double NuXv = dot_product(Nu, Xv);
  double NvXv = dot_product(Nv, Xv);
  
  // The derivative of 'g'
  double &g = xAtom->G.g;
  double dgdN = dAtom->G.g = 
    2 * ( NuXu * G1[1][1] + NvXv * G1[0][0] - (NuXv + NvXu) * G1[0][1] );
  double g2 = g * g;

  // The derivative of the contravariant tensor
  double (*z)[2] = dAtom->G.xContravariantTensor; 
  z[0][0] = (2 * NvXv * g - G1[1][1] * dgdN) / g2;
  z[1][1] = (2 * NuXu * g - G1[0][0] * dgdN) / g2;
  z[0][1] = z[1][0] = - ((NuXv + NvXu) * g - G1[1][0] * dgdN) / g2;

  // Compute the right hand side
  *b = - (z[0][0] * Fu * Fu + 2.0 * z[0][1] * Fu * Fv + z[1][1] * Fv * Fv);
}

void
FDBorderSite::ComputeVariationalDerivativeRho(double *Y, double *A, double *b, 
  MedialAtom *xAtom, MedialAtom *dAtom)
{
  // First, let's compute grad phi
  double Fu = _du * (Y[ xNeighbor[1] ] - Y[ xNeighbor[3] ]);
  double Fv = _dv * (Y[ xNeighbor[2] ] - Y[ xNeighbor[4] ]);

  // Next, compute the weights for Hu and Hv and H
  double (*G1)[2] = xAtom->G.xCovariantTensor;
  double (*G2)[2] = xAtom->G.xContravariantTensor;
  double whu = _du * (G2[0][0] * Fu + G2[0][1] * Fv);
  double whv = _dv * (G2[1][0] * Fu + G2[1][1] * Fv);
  double wh = -2.0;

  // Clear all the entries to zero
  for(unsigned int j = 0; j < nDistinctSites; j++)
    A[j] = 0.0;
  A[ xEntry[0] ] += wh;
  A[ xEntry[1] ] += whu;
  A[ xEntry[2] ] += whv;
  A[ xEntry[3] ] -= whu;
  A[ xEntry[4] ] -= whv;

  // The right hand side is just 0
  *b = 0.0;
}

/** Initialize the border site */
void FDBorderSite::SetGeometry(GeometryDescriptor *g, double)
{
  
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

double SparseLinearTest(unsigned int n, unsigned int *rowIndex, unsigned int *colIndex, 
	double *values, double *x, double *y,double *b)
{
  SparseMultiply(n, rowIndex, colIndex, values, x, y);
  double maxdiff = 0.0;
  for(unsigned int i = 0; i < n; i++)
    if(fabs(b[i] - y[i]) > maxdiff)
      maxdiff = fabs(b[i] - y[i]);
  return maxdiff;
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
  b.set_size(nSites);
  y.set_size(nSites);
  eps.set_size(nSites);
  zTest.set_size(nSites);
  
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
  xInitSoln.set_size(nSites);
  SetDefaultInitialGuess(0.001);

  // Initialize the atom grid to a Cartesian grid
  xGrid = new CartesianMedialAtomGrid(m, n);
}

void
MedialPDESolver
::SetDefaultInitialGuess(double xMagnitude)
{
  for(size_t k = 0; k < nSites; k++)
    xInitSoln[k] = xMagnitude * xSites[k]->GetInitialValue();
}

/**
void MedialPDESolver::ComputeInitialGuess()
{
  // In this method, we want to assign initial values so that grad R is as 
  // close to one as possible at the border sites
  double xMaxDist = 0.0;
  for(size_t i = 0; i < m; i++) for(size_t j = 0; j < n; j++)
    {
    // Get the site and atom locations
    unsigned int iGrid = xGrid->GetAtomIndex(i, j);
    unsigned int iSite = xSiteIndex[i][j];

    // If the site is internal, simply set the value to equal one
    if(xSites[iSite]->IsBorderSite())
      {
      unsigned int iOpp = i, jOpp = j;
      if(i == 0) iOpp = 1;
      if(i == m-1) iOpp = m-2;
      if(j == 0) jOpp = 1;
      if(j == n-1) jOpp = n-2;

      size_t iGridOpp = xGrid->GetAtomIndex(iOpp, jOpp);
      SMLVec3d X1 = xAtoms[iGridOpp].X;
      SMLVec3d X0 = xAtoms[iGrid].X;
      double xDist = (X0-X1).magnitude();
      xInitSoln[iSite] = - xDist;

      if(xMaxDist < xDist) xMaxDist = xDist;
      }
    else 
      xInitSoln[iSite] = 0.0;
    }

  // Correct for negative distances
  for(size_t k = 0; k < nSites; k++)
    {
    double r = xInitSoln[k] + xMaxDist + 0.01;
    xInitSoln[k] = r * r;
    }
}
*/

void 
MedialPDESolver
::SetSolutionAsInitialGuess()
{
  xInitSoln = y;
  flagReuseLastSolution = true;
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
    xSurface->EvaluateAtGridIndex(i, j, 0, 0, 0, 3, xAtom.X.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 1, 0, 0, 3, xAtom.Xu.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 0, 1, 0, 3, xAtom.Xv.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 2, 0, 0, 3, xAtom.Xuu.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 1, 1, 0, 3, xAtom.Xuv.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 0, 2, 0, 3, xAtom.Xvv.data_block());

    // Compute the differential geometric tensors
    xAtom.ComputeDifferentialGeometry();

    // Compute the normal vector
    xAtom.ComputeNormalVector();

    // Compute the laplacian of R 
    xSurface->EvaluateAtGridIndex(i, j, 0, 0, 3, 1, &xAtom.xLapR);
    if(xAtom.xLapR > 0 && xSites[iSite]->IsBorderSite()) 
      cout << xAtom.xLapR << " ! " << flush;

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
      xSites[iLocal]->ComputePhiJet(ySolution, xAtom.F, xAtom.Fu, xAtom.Fv);

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

double MedialPDESolver::SolveOnce(double delta)
{
  size_t i, j, k;
  double epsMax, bMax;
  
  // Initialize epsilon to zero
  eps.fill(0.0);
  
  // Copy the initial solution to the current solution
  y = xInitSoln;

  // We are now ready to perform the Newton loop
  bool flagComplete = false;
  for(unsigned int iIter = 0; iIter < 16 && !flagComplete; iIter++)
    {
    // Compute the A matrix and the b vector at each node
    for(k = 0; k < nSites; k++)
      {
      // Compute the non-zero values of A for this row
      xSites[k]->ComputeDerivative(y.data_block(), xSparseValues + xRowIndex[k] - 1, k+1);

      // Compute the value of b
      b[k] = -xSites[k]->ComputeEquation(y.data_block());
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
        b.data_block(), eps.data_block(), &ERROR);
      }

    // Only perform the solution step
    PHASE=23;
    pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, 
      xSparseValues, (int *) xRowIndex, (int *) xColIndex, 
      NULL, &NRHS, IPARM, &MSGLVL, 
      b.data_block(), eps.data_block(), &ERROR);

    // Get the largest error (eps)
    epsMax = eps.inf_norm();
    bMax = b.inf_norm();

    // Append the epsilon vector to the result
    y += eps;

    // Print the statistics
    /*
    cout << "-----------" << endl;
    cout << "Step " << iIter << ": " << endl;
    cout << "  Largest Epsilon: " << epsMax << endl;
    cout << "  Largest Eqn Error: " << bMax << endl; 
    
    // Test the matrix result
    SparseLinearTest(nSites, xRowIndex, xColIndex, xSparseValues, eps, zTest, b);
    double zMax = fabs(zTest[myidamax(nSites, zTest, 1)]);
    cout << "  Largest Solver Error: " << zMax << endl;
    */

    // Convergence is defined when epsilon is smaller than some threshold
    if(bMax < delta) 
      break;
    }

  return bMax;
}

void
MedialPDESolver
::Solve(double delta)
{
  unsigned int i, j, k;
  vnl_vector<double> xBestInit(nSites);

  // Intialize the sites
  InitializeSiteGeometry();

  // Try solving with the current initial guess
  double xBestGuess = SolveOnce(delta);
  if(xBestGuess > delta)
    {
    // Save the initial solution for all that it's worth
    xBestInit = xInitSoln;

    // Try using the following scale factors to get a decent solution
    double xTest = 1.0e-6, xScale = 2.0, nGuesses = 20;
    for(size_t iGuess = 0; iGuess < nGuesses; iGuess++)
      {
      
      // Try solving using the current guess
      SetDefaultInitialGuess(xTest);
      double xGuess = SolveOnce(delta);
      
      cout << "PDE Solver Initialization Failure; Initializing with " << xTest 
        << "; Result: " << xGuess << endl;
      
      if(xGuess < xBestGuess)
        {
        xBestGuess = xGuess;
        xBestInit = xInitSoln;
        }
      
      if(xGuess < delta)
        break;

      // Scale up to the next guess
      xTest *= xScale;
      }
    }
  
  // If we never succeeded, use the best initialization that we could find
  if(xBestGuess > delta)
    { 
    cerr << "Complete PDE Solver Failure!" << endl; 
    xInitSoln = xBestInit;
    }

  // Reconstruct the medial atoms
  ReconstructAtoms(y.data_block());
}

void ComputeMedialAtomBoundaryDerivative(
  MedialAtom *xAtom, MedialAtom *dAtom, bool isEdge)
{
  // Get the relevant elements of the atoms
  SMLVec3d &X = xAtom->X, &Xu = xAtom->Xu, &Xv = xAtom->Xv;
  SMLVec3d &N = dAtom->X, &Nu = dAtom->Xu, &Nv = dAtom->Xv;
  
  // Get the elements of the first fundamental form and its derivative
  double g11 = xAtom->G.xContravariantTensor[0][0];
  double g12 = xAtom->G.xContravariantTensor[0][1];
  double g22 = xAtom->G.xContravariantTensor[1][1];
  double z11 = dAtom->G.xContravariantTensor[0][0];
  double z12 = dAtom->G.xContravariantTensor[0][1];
  double z22 = dAtom->G.xContravariantTensor[1][1];

  // Get the g's
  double g = xAtom->G.g; double z = dAtom->G.g;

  // Get the partials of Phi and its variational derivative
  double F = xAtom->F, Fu = xAtom->Fu, Fv = xAtom->Fv;
  double H = dAtom->F, Hu = dAtom->Fu, Hv = dAtom->Fv;
  
  // This is the derivative of the normal vector
  dAtom->N = (vnl_cross_3d(Xu, Nv) + vnl_cross_3d(Nu, Xv)) / sqrt(g) - 
    ( 0.5 * z / g ) * xAtom->N;

  // We will compute several intermediate terms
  // This is the derivative of - grad phi
  SMLVec3d T2 = -(
    z11 * Fu * Xu + z12 * (Fu * Xv + Fv * Xu) + z22 * Fv * Xv +
    g11 * Hu * Xu + g12 * (Hu * Xv + Hv * Xu) + g22 * Hv * Xv +
    g11 * Fu * Nu + g12 * (Fu * Nv + Fv * Nu) + g22 * Fv * Nv );

  // Address the edge case first
  if(isEdge) 
    {
    dAtom->xBnd[1].X = dAtom->xBnd[0].X = N + 0.5 * T2;
    return;
    }

  // Compute the intermediate terms:
  // This is the coefficient of the normal before differentiation
  double T1 = sqrt(4.0 * F - 
    (g11 * Fu * Fu + 2.0 * g12 * Fu * Fv + g22 * Fv * Fv));

  // This is the ugly term
  double T3 = 4.0 * H - ( 
    z11 * Fu * Fu + 2.0 * z12 * Fu * Fv + z22 * Fv * Fv +
    2.0 * ( g11 * Hu * Fu + g12 * (Hu * Fv + Hv * Fu) + g22 * Hv * Fv ) );

  // The derivative of the out-out-tangent term
  SMLVec3d T4 = ( 0.5 * T3 / T1 ) * xAtom->N + T1 * dAtom->N;

  // Compute the Y atoms
  dAtom->xBnd[0].X = N + 0.5 * ( T2 - T4 );
  dAtom->xBnd[1].X = N + 0.5 * ( T2 + T4 );
}

void
MedialPDESolver
::ComputeVariationalDerivative(
  IHyperSurface2D *xVariation, bool flagRhoVariation, MedialAtom *dAtoms)
{
  size_t i, j;

  // First thing is to define the PDE's that will give us the variational
  // derivative of the function phi.
  xVariation->SetEvaluationGrid(m, n, xGridU, xGridV);

  // Compute the jet of nu at each site
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    {
    // Get the index of the site
    size_t iGrid = xGrid->GetAtomIndex(i, j);
    size_t iSite = xSiteIndex[i][j];

    // Access the medial atom underneath
    MedialAtom &dAtom = dAtoms[iGrid];

    // Set the atoms' domain coordinates
    dAtom.u = xGridU[i]; dAtom.v = xGridV[j];

    // Split depending on the type of variation
    if(flagRhoVariation)
      {
      // Compute the rho derivative with respect to the variation
      xVariation->EvaluateAtGridIndex(i, j, 0, 0, 3, 1, &dAtom.xLapR);

      // Prepare the matrices for the linear solver
      xSites[iSite]->ComputeVariationalDerivativeRho(
        y.data_block(), xSparseValues + xRowIndex[iSite] - 1, &b[iSite],
        &xAtoms[iGrid], &dAtom);
      }
    else
      {
      // Compute the derivative of the surface with respect to the variation
      xVariation->EvaluateAtGridIndex(i, j, 0, 0, 0, 3, dAtom.X.data_block());
      xVariation->EvaluateAtGridIndex(i, j, 1, 0, 0, 3, dAtom.Xu.data_block());
      xVariation->EvaluateAtGridIndex(i, j, 0, 1, 0, 3, dAtom.Xv.data_block());
      xVariation->EvaluateAtGridIndex(i, j, 2, 0, 0, 3, dAtom.Xuu.data_block());
      xVariation->EvaluateAtGridIndex(i, j, 1, 1, 0, 3, dAtom.Xuv.data_block());
      xVariation->EvaluateAtGridIndex(i, j, 0, 2, 0, 3, dAtom.Xvv.data_block());

      // Prepare the matrices for the linear solver
      xSites[iSite]->ComputeVariationalDerivativeX(
        y.data_block(), xSparseValues + xRowIndex[iSite] - 1, &b[iSite],
        &xAtoms[iGrid], &dAtom);
      }
    }

  // Solve the partial differential equation
  int MAXFCT = 1, MNUM = 1, PHASE = 13, N = nSites, NRHS = 1, MSGLVL = 0, ERROR = 0; 
  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, 
    xSparseValues, (int *) xRowIndex, (int *) xColIndex, 
    NULL, &NRHS, IPARM, &MSGLVL, 
    b.data_block(), eps.data_block(), &ERROR);

  // For each atom, compute the boundary derivatives
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    {
    // Get the index of the site
    size_t iGrid = xGrid->GetAtomIndex(i, j);
    size_t iSite = xSiteIndex[i][j];

    // Access the medial atom underneath
    MedialAtom &dAtom = dAtoms[iGrid];

    // Compute the gradient of phi for the new atom
    xSites[iSite]->ComputePhiJet(eps.data_block(), dAtom.F, dAtom.Fu, dAtom.Fv);

    // Compute the rest of the atom derivative
    ComputeMedialAtomBoundaryDerivative(
      &xAtoms[iGrid], &dAtom, xSites[iSite]->IsBorderSite());
    }
}

/** Remove this later */
double MedialPDESolver
::ComputeLaplaceBeltrami(unsigned int i, unsigned int j, vnl_vector<double> xField)
{
  // Get the index of the site
  unsigned int iGrid = xGrid->GetAtomIndex(i, j);
  unsigned int iSite = xSiteIndex[i][j];

  cout << "Computing Lappy at " << iSite << endl;

  return xSites[iSite]->ComputeEquation(xField.data_block());
  
  // Compute the Laplacian
  // return xSites[iSite]->ComputeLaplaceBeltrami(&xAtoms[iGrid], xField.data_block());
}
