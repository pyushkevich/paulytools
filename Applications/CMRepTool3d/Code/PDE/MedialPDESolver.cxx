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

/*
void SparseLinearTest(unsigned int n, unsigned int *rowIndex, unsigned int *colIndex, 
	double *values, double *x, double *y,double *b)
{
  SparseMultiply(n, rowIndex, colIndex, values, x, y);
  mydaxpy(n, -1.0, b, 1, y, 1);
}
*/

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
