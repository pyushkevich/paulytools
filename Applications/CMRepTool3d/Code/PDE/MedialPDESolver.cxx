#include "MedialPDESolver.h"
#include <cmath>
#include <algorithm>
#include "mspline.h"
#include "fastmath.h"

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
  return 0.1 * sqrt(x * x + y * y);
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
void FDInternalSite::SetGeometry(GDescriptor *g, double rho)
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
void FDBorderSite::SetGeometry(GDescriptor *g, double)
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

  // Compute the initial solution
  xInitSoln = new double[nSites];
  for(iSite = 0; iSite < nSites; iSite++)
    xInitSoln[iSite] = xSites[iSite]->GetInitialValue();

  // Initialize our three vectors
  b = new double[nSites];
  y = new double[nSites];
  eps = new double[nSites];
  zTest = new double[nSites];
  
  // Initialize the medial atom array
  xAtoms = new MedialPoint[nSites];
  xGeometryDescriptors = new GDescriptor[nSites];

  // Initialize the PARDISO solver
  MTYPE = 11; // Nonsymmetric real
  memset(IPARM, 0, sizeof(int) * 64);
  pardisoinit_(PT,&MTYPE,IPARM);
  IPARM[2] = 1;
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

      // Access the medial atom underneath
      MedialPoint *mp = xAtoms + iSite;
      mp->u = u; mp->v = v;

      // Compute the surface jet and the laplacian
      problem->ComputeJet2(mp);
      // problem->ComputeJet2(u, v, X, Xu, Xv, Xuu, Xuv, Xvv);

      // Compute the geometrical properties
      xGeometryDescriptors[iSite].SetJet(
        mp->F.data_block(), mp->Fu.data_block(), mp->Fv.data_block(), 
        mp->Fuu.data_block(), mp->Fuv.data_block(), mp->Fvv.data_block());

      // Compute the solution at this point
      xSites[iSite]->SetGeometry(
        &xGeometryDescriptors[iSite], problem->ComputeLaplacian(u,v));
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
      b[k] = -xSites[k]->ComputeEquation(y);
      }

    // Debug 
    // DumpSparseMatrix(nSites, xRowIndex, xColIndex, xSparseValues);

    // Solve the equation for X
    int MAXFCT = 1, MNUM = 1, PHASE = 11, N = nSites, NRHS = 1, MSGLVL = 1, ERROR = 0; 
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
    cout << "-----------" << endl;
    cout << "Step " << iIter << ": " << endl;
    cout << "  Largest Epsilon: " << epsMax << endl;
    cout << "  Largest Eqn Error: " << bMax << endl;
    // cout << "  Largest Solver Error: " << zMax << endl;

    // Convergence is defined when epsilon is smaller than some threshold
    flagComplete = (epsMax < delta);
    }

  // Once the iterations are complete, reconstruct the atoms
  for(k = 0; k < nSites; k++)
    {
    MedialPoint *mp = xAtoms + k;

    // The case where phi is negative is undefined
    if(y[k] > 0)
      {
      // Compute the grad r vector using finite differences
      double R, Ru, Rv;
      xSites[k]->ComputeRJet(y, R, Ru, Rv);
      mp->F[3] = R; mp->Fu[3] = Ru; mp->Fv[3] = Rv;

      // Compute the boundary properties of the medial point
      ComputeMedialAtom(mp, &xGeometryDescriptors[k]);
      }
    else
      {
      // What are we supposed to do?
      mp->F[3] = mp->Fu[3] = mp->Fv[3] = 0.0;
      }
    }
}

bool 
MedialPDESolver
::ComputeMedialAtom(MedialPoint *p, GDescriptor *gd)
{
  // The medial point should have the X 2-jet and the R 1-jet already set
  // We have to compute grad R and the boundary sites

  // Get a handle on the variables we will need
  SMLVec3f &X   = p->X();
  SMLVec3f &Xu  = p->Xu();
  SMLVec3f &Xv  = p->Xv();
  float &R   = p->R();
  float &Ru  = p->Ru();
  float &Rv  = p->Rv();

  // Elements of the first fundamental form
  float E = p->IE = gd->xCovariantTensor[0][0];
  float F = p->IF = gd->xCovariantTensor[0][1];
  float G = p->IG = gd->xCovariantTensor[1][1];
  
  // Terms going into gradR
  float CXu = Ru * gd->xContravariantTensor[0][0] 
    + Rv * gd->xContravariantTensor[0][1]; 
  float CXv = Ru * gd->xContravariantTensor[0][1] 
    + Rv * gd->xContravariantTensor[1][1];
    
  // Compute Grad R
  p->GradR[0] = Xv[0] * CXv + Xu[0] * CXu;
  p->GradR[1] = Xv[1] * CXv + Xu[1] * CXu;
  p->GradR[2] = Xv[2] * CXv + Xu[2] * CXu;

  // Compute squared length of GradR
  double xMagGradR2 = Ru * CXu + Rv * CXv;
  p->sinTheta2 = 1.0f - xMagGradR2;
  
  // Correct the floating point / solver error
  if(fabs(p->sinTheta2) < 1e-5)
    p->sinTheta2 = 0.0;
  
  // Compute the boundary sites
  if (p->sinTheta2 >= 0.0)
    {
    // Compute the derivatives of sine theta
    float sT = sqrtf(p->sinTheta2);

    // Compute the normal and tangent components of the boundaries
    SMLVec3f CN = p->N() * sT;

    // Compute the position and normals of boundary points
    p->bp[0].N = - p->GradR - CN;
    p->bp[1].N = - p->GradR + CN;
    p->bp[0].X = X + p->bp[0].N * R;
    p->bp[1].X = X + p->bp[1].N * R;

    // Return success
    return true;
    }

  // The point is outside of the boundary
  cout << "Boundary Failure: sin(theta)^2 = " << p->sinTheta2 << endl;
  return false;
}

float MedialPDESolver
::IntegrateBoundaryMeasure(BoundaryMeasure *bnd, float &area)
{
  // Integration over patches. First, compute the boundary measure at
  // each point on the medial surface.
  float *bndValues[2], *xAreas;
  bndValues[0]= new float[nSites];
  bndValues[1]= new float[nSites];

  // This vector stores the area assigned to each boundary point
  xAreas = new float[nSites];

  // Compute the boundary measure at each point
  for(unsigned int iSite = 0; iSite < nSites; iSite++)
    {
    // For boundary sites, only one side is used
    if(xSites[iSite]->IsBorderSite())
      {
      // Compute the crest measure
      bndValues[0][iSite] = bndValues[1][iSite] 
        = bnd->computeCrestBoundaryMeasure(xAtoms[iSite]);
      }
    else
      {
      // Compute the boundary measure
      bndValues[0][iSite] = bnd->computeBoundaryMeasure(xAtoms[iSite], 0);
      bndValues[1][iSite] = bnd->computeBoundaryMeasure(xAtoms[iSite], 1);
      }
    }

  // Integrate the surface area
  float xArea = 0.0, xMeasure = 0.0;
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
        float A1 = triangleArea( xAtoms[i00].F.data_block(), 
          xAtoms[i01].F.data_block(), xAtoms[i11].F.data_block());
        float A2 = triangleArea(xAtoms[i00].F.data_block(), 
          xAtoms[i11].F.data_block(), xAtoms[i10].F.data_block());

        // Integrate the measure
        xMeasure += A1 * (bndValues[d][i00] + bndValues[d][i01] + bndValues[d][i11]); 
        xMeasure += A2 * (bndValues[d][i00] + bndValues[d][i11] + bndValues[d][i10]); 
        xArea += (A1 + A2);
        }

  // Scale the measure by three
  xMeasure /= 3.0;

  // Clean up
  delete xAreas;
  delete bndValues[0];
  delete bndValues[1];
}
