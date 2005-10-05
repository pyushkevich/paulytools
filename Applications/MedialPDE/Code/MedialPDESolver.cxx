#include "MedialPDESolver.h"
#include <cmath>
#include <algorithm>
#include "smlmath.h"

#include "MedialPDEMasks.h"

using namespace std;

void SparseMultiply(
  unsigned int n, unsigned int *rowIndex, unsigned int *colIndex, 
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

double SparseLinearTest(unsigned int n, unsigned int *rowIndex, 
  unsigned int *colIndex, double *values, double *x, double *y,double *b)
{
  SparseMultiply(n, rowIndex, colIndex, values, x, y);
  double maxdiff = 0.0;
  for(unsigned int i = 0; i < n; i++)
    if(fabs(b[i] - y[i]) > maxdiff)
      maxdiff = fabs(b[i] - y[i]);
  return maxdiff;
}

void DumpSparseMatrix(unsigned int n, unsigned int *rowIndex, 
  unsigned int *colIndex, double *values)
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

MedialPDESolver
::MedialPDESolver(size_t nu, size_t nv, double xScale, size_t pu, size_t pv)
{
  // Compute the total number of grid points
  size_t tu = nu + 2 * pu;
  size_t tv = nv + 2 * pv;

  // Allocate the initialization grids
  Vec uuGrid(tu), vvGrid(tv);

  // Set the endpoints
  uuGrid[0] = vvGrid[0] = 0.0;
  uuGrid[tu-1] = vvGrid[tv-1] = 1.0;

  // Set the regularly spaced points
  double du = 1.0 / (nu - 1);
  double dv = 1.0 / (nv - 1);
  size_t iu, iv;

  double zu = du;
  for(iu = 1; iu <= nu-2; iu++, zu+=du)
    uuGrid[iu + pu] = zu;

  double zv = dv;
  for(iv = 1; iv <= nv-2; iv++, zv+=dv)
    vvGrid[iv + pv] = zv;

  // Set the irregularly spaced points
  for(iu = 0; iu < pu; iu++)
    {
    double duScaled = du * pow(xScale, iu+1.0); 
    uuGrid[pu - iu] = duScaled;
    uuGrid[(tu - 1) - (pu - iu)] = 1.0 - duScaled; 
    }
  
  for(iv = 0; iv < pv; iv++)
    {
    double dvScaled = dv * pow(xScale, iv+1.0); 
    vvGrid[pv - iv] = dvScaled;
    vvGrid[(tv - 1) - (pv - iv)] = 1.0 - dvScaled; 
    }

  // Apply random jitter
  /*
  for(iu = 1; iu <= nu-2; iu++)
    uuGrid[iu + pu] += 0.45 * (-1.0 + rand() * 2.0 / RAND_MAX) * du;

  for(iv = 1; iv <= nv-2; iv++)
    vvGrid[iv + pv] += 0.45 * (-1.0 + rand() * 2.0 / RAND_MAX) * dv;
  */

  Initialize(uuGrid, vvGrid);
}

MedialPDESolver
::MedialPDESolver(const Vec &uGrid, const Vec &vGrid)
{
  Initialize(uGrid, vGrid);
}

void 
MedialPDESolver
::Initialize(const Vec &uGrid, const Vec &vGrid)
{
  // Copy the grid vectors
  this->uGrid = uGrid; this->vGrid = vGrid;

  // Record the size of the grid
  m = uGrid.size(); n = vGrid.size();

  // Number of sites
  nSites = m * n;

  // Initialize the mask and site arrays
  xSites.resize(nSites, NULL);
  xMasks.resize(nSites, NULL);

  // Set the size of the size index
  xSiteIndex.set_size(m, n);

  // Set up the initial guess
  xDefaultInitSoln.set_size(m, n);

  // Initialize all the masks and sites
  size_t iSite = 0, i, j;
  for(i = 0; i < m; i++) for(j = 0; j < n; j++) 
    {
    // Flags indicating adjacence to a border
    bool uBorder = (i == 0 || i == m-1), vBorder = (j == 0 || j == n-1);

    // cout << i << ", " << j << " : ";

    // Check if we are at the border
    if(uBorder || vBorder)
      {
      if(uBorder && vBorder)
        {
        // Create a corner mask
        // cout << "Corner Mask" << endl;
        xMasks[iSite] = new CornerFDMask();
        }
      else
        {
        // Create a border mask
        // cout << "Border Mask" << endl;
        xMasks[iSite] = new BorderFDMask();
        if(vBorder)
          xMasks[iSite]->TransposeMask();
        }

      // In either case, flip the mask over
      xMasks[iSite]->FlipMask(i != 0, j != 0);
      }

    // We are at an internal site
    else
      {
      // Check the difference in grid spacing
      double ddu = fabs((uGrid[i+1] - uGrid[i]) - (uGrid[i] - uGrid[i-1]));
      double ddv = fabs((vGrid[j+1] - vGrid[j]) - (vGrid[j] - vGrid[j-1]));
      // Check if the grid is uniform
      if(ddu < 1.0e-13 && ddv < 1.0e-13)
        {
        // Create a uniform mask
        // cout << "Uniform Mask" << endl;
        xMasks[iSite] = new UniformFDMask();
        }
      else
        {
        // Create a non-uniform mask, and orient it correctly
        // cout << "Nonuniform Mask" << endl;
        xMasks[iSite] = new DoublyNonuniformFDMask();
        xMasks[iSite]->FlipMask( i > m/2, j > n/2 );
        }
      }

    // Set the location of the mask and specify the grid
    xMasks[iSite]->SetLocation(i, j);
    xMasks[iSite]->SortNodes();
    xMasks[iSite]->ComputeWeights(uGrid.data_block(), vGrid.data_block());

    // Compute the initial solution value as the distance from the nearest
    // edge
    // double uMin = i > m/2 ? uGrid[m-(i+1)] - uGrid[m-(i+2)] : uGrid[i+1] - uGrid[i];
    // double vMin = j > n/2 ? vGrid[n-(j+1)] - vGrid[n-(j+2)] : vGrid[j+1] - vGrid[j];
    
    // xDefaultInitSoln[i][j] = sqrt( uMin * uMin + vMin * vMin );
    xDefaultInitSoln[i][j] = 1.0;

    // xMasks[iSite]->PrintReport();

    // Make sure that every mask is in bounds
    for(size_t k = 0; k < xMasks[iSite]->Size(); k++)
      {
      int ku = xMasks[iSite]->GetOffsetU(k) + xMasks[iSite]->GetLocationU();
      int kv = xMasks[iSite]->GetOffsetV(k) + xMasks[iSite]->GetLocationV();
      if(ku < 0 || ku >= m || kv < 0 || kv >= n) 
        cout << "Site " << i << ", " << j << " is out of bounds!" << endl;
      }

    // Create a site wrapping around the mask
    if(uBorder || vBorder)
      xSites[iSite] = new FDBorderSite(xMasks[iSite]);
    else
      xSites[iSite] = new FDInternalSite(xMasks[iSite]);

    // Associate the raw offset with the u-v index
    xSiteIndex[i][j] = iSite++;
    }

  // Initialize the sparse matrix row index (use Fortran-based indexing)
  xRowIndex = new int[nSites + 1];
  xRowIndex[0] = 1;
  for(iSite = 0; iSite < nSites; iSite++)
    xRowIndex[iSite+1] = xRowIndex[iSite] + xMasks[iSite]->Size();
  
  // Initialize the column index
  nSparseEntries = xRowIndex[nSites] - 1;
  xColIndex = new int[nSparseEntries];
  
  // Initialize the sparse data
  xSparseValues = new double[nSparseEntries];
  memset(xSparseValues, 0, sizeof(double) * nSparseEntries);

  // Compute the column index and the sparse data
  for(iSite = 0; iSite < nSites; iSite++)
    {
    // Number of neighbors in this node
    size_t nEntries = xRowIndex[iSite+1] - xRowIndex[iSite];

    for(unsigned int iEntry = 0; iEntry < nEntries; iEntry++)
      {
      // Get the coordinates of the given neighbor node
      int iu = xMasks[iSite]->GetLocationU() + xMasks[iSite]->GetOffsetU(iEntry);
      int iv = xMasks[iSite]->GetLocationV() + xMasks[iSite]->GetOffsetV(iEntry);
      int iNeighbor = xSiteIndex[iu][iv];

      // Set the column entry in the sparse matrix
      xColIndex[ xRowIndex[iSite] + iEntry - 1 ] = iNeighbor + 1;
      }
    }

  // Initialize our three vectors
  b.set_size(m, n);
  y.set_size(m, n);
  eps.set_size(m, n);
  dy.set_size(m, n);
  zTest.set_size(m, n);

  // Initialize the medial atom array and atom grid
  xAtoms = new MedialAtom[nSites];
  xGrid = new CartesianMedialAtomGrid(m, n);

  // Compute the initial solution
  SetDefaultInitialGuess(1);
}

void
MedialPDESolver
::SetDefaultInitialGuess(double xMagnitude)
{
  // Set the initial guess to default values provided by the sites
  xInitSoln = xMagnitude * xDefaultInitSoln;
}

void 
MedialPDESolver
::SetSolutionAsInitialGuess()
{
  // Make the current solution the initial guess. This makes Newton converge
  // quickly when exploring near a point
  xInitSoln = y;
  flagReuseLastSolution = true;
}

void 
MedialPDESolver
::SetMedialSurface(IMutableHyperSurface2D *xSurface)
{
  // Remember the surface
  this->xSurface = xSurface;

  // Tell the surface where it will be evaluated
  xSurface->SetEvaluationGrid(uGrid, vGrid);
}

void MedialPDESolver::InitializeSiteGeometry()
{
  size_t i, j;
  
  /*
  // Generate a matrix of x, y, z values
  Mat vx(m, n), vy(m, n), vz(m, n);

  // Sample the surface along the current grid
	for(i = 0; i < m; i++) for(j = 0; j < n; j++)
	  {
	  double X[3];
	  xSurface->EvaluateAtGridIndex(i, j, 0, 0, 0, 3, X);
	  vx[i][j] = X[0]; vy[i][j] = X[1]; vz[i][j] = X[2];
	  }
  */

  // Initialize each site with the current surface properties
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    {
    // Get the index of the site
    unsigned int iGrid = xGrid->GetAtomIndex(i, j);
    unsigned int iSite = xSiteIndex[i][j];

    // Access the medial atom underneath
    MedialAtom &xAtom = xAtoms[iGrid];

    // Set the atoms' domain coordinates
    xAtom.u = uGrid[i]; xAtom.v = vGrid[j];
    double u = uGrid[i], v = vGrid[j];
    
    // Compute the surface jet and the laplacian    
    xSurface->EvaluateAtGridIndex(i, j, 0, 0, 0, 3, xAtom.X.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 1, 0, 0, 3, xAtom.Xu.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 0, 1, 0, 3, xAtom.Xv.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 2, 0, 0, 3, xAtom.Xuu.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 1, 1, 0, 3, xAtom.Xuv.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 0, 2, 0, 3, xAtom.Xvv.data_block());

		/*    
    xAtom.X[0] = xMasks[iSite]->ComputeTwoJet(vx, xAtom.Xu[0], xAtom.Xv[0], xAtom.Xuu[0], xAtom.Xuv[0], xAtom.Xvv[0]);
    xAtom.X[1] = xMasks[iSite]->ComputeTwoJet(vy, xAtom.Xu[1], xAtom.Xv[1], xAtom.Xuu[1], xAtom.Xuv[1], xAtom.Xvv[1]);
    xAtom.X[2] = xMasks[iSite]->ComputeTwoJet(vz, xAtom.Xu[2], xAtom.Xv[2], xAtom.Xuu[2], xAtom.Xuv[2], xAtom.Xvv[2]);
    */
/*
		double u = uGrid[i], v = vGrid[j];
    xAtom.X[0]   = 10 * u; 
    xAtom.X[1]   = 10 * v; 
    xAtom.X[2]   = sin(u * v) + cos(u*u + v*v);

    xAtom.Xu[0]  = 10;       xAtom.Xu[1]  = 0;        xAtom.Xu[2]  = v*cos(u*v) - 2*u*sin(u*u+v*v);
    xAtom.Xv[0]  = 0;        xAtom.Xv[1]  = 10;       xAtom.Xv[2]  = u*cos(u*v) - 2*v*sin(u*u+v*v);
    xAtom.Xuu[0] = 0;        xAtom.Xuu[1] = 0;        xAtom.Xuu[2] = -v*v*sin(u*v) - 2*sin(u*u+v*v) - 4*u*u*cos(u*u+v*v);
    xAtom.Xuv[0] = 0;        xAtom.Xuv[1] = 0;        xAtom.Xuv[2] = cos(u*v)-u*v*sin(u+v) - 4*u*v*cos(u*u+v*v);
    xAtom.Xvv[0] = 0;        xAtom.Xvv[1] = 0;        xAtom.Xvv[2] = -u*u*sin(u+v) - 2*sin(u*u+v*v) - 4*v*v*cos(u*u+v*v);
  */

    // Compute the differential geometric tensors
    xAtom.ComputeDifferentialGeometry();

    // Compute the normal vector
    xAtom.ComputeNormalVector();

    // Compute the laplacian of R 
    xSurface->EvaluateAtGridIndex(i, j, 0, 0, 3, 1, &xAtom.xLapR);
    
    // if(xAtom.xLapR > 0 && xSites[iSite]->IsBorderSite()) 
    //  cout << xAtom.xLapR << " ! " << flush;

    // Compute the solution at this point
    xSites[iSite]->SetGeometry( &xAtom.G, xAtom.xLapR);
    }
    
  /*
  // Define test function
	Mat phitest(m, n, 0.0);
	for(i = 0; i < m; i++) for(j = 0; j < n; j++)
		{
		double u = uGrid[i], v = vGrid[j];
		phitest(i, j) = sinh(u) + cosh(v * u);
		}
		
	Mat LB1(m, n, 0.0), LB2(m, n, 0.0);
	for(i = 1; i < m-1; i++) for(j = 1; j < n-1; j++)
		{
		LB1[i][j] = xSites[xSiteIndex[i][j]]->ComputeEquation(phitest) + xAtoms[xSiteIndex[i][j]].xLapR;
		LB2[i][j] = EstimateLBOperator(phitest, i, j);
		}

  cout << "SURFACE INFO" << endl;
  xSurface->PrintReport();

  MedialAtom B;
  B.X[0] = xMasks[xSiteIndex[8][8]]->ComputeTwoJet(vx, B.Xu[0], B.Xv[0], B.Xuu[0], B.Xuv[0], B.Xvv[0]);
  B.X[1] = xMasks[xSiteIndex[8][8]]->ComputeTwoJet(vy, B.Xu[1], B.Xv[1], B.Xuu[1], B.Xuv[1], B.Xvv[1]);
  B.X[2] = xMasks[xSiteIndex[8][8]]->ComputeTwoJet(vz, B.Xu[2], B.Xv[2], B.Xuu[2], B.Xuv[2], B.Xvv[2]);
		
  cout << "AT 0.04, 0.04 " << endl;
  MedialAtom &A = xAtoms[xSiteIndex[8][8]];
  cout << A.u << " " << A.v << endl;
  cout << "X : " << A.X << " " << B.X << endl;
  cout << "Xu : " << A.Xu << " " << B.Xu << endl;
  cout << "Xv : " << A.Xv << " " << B.Xv << endl;
  cout << "Xuu : " << A.Xuu << " " << B.Xuu << endl;
  cout << "Xuv : " << A.Xuv << " " << B.Xuv << endl;
  cout << "Xvv : " << A.Xvv << " " << B.Xvv << endl;
		
	cout << "ANALYTIC: " << LB1[8][8] << endl;
	cout << "NUMERIC: " << LB2[8][8] << endl;
  */
}

bool
MedialPDESolver
::ReconstructAtoms(const Mat &ySolution)
{
  // Keep track of invalid atoms
  bool flagAllAtomsValid = true;
  
  // Once the iterations are complete, reconstruct the atoms
  for(unsigned int i = 0; i < m; i++) for(unsigned int j = 0; j < n; j++)
    {
    // Map to a grid index and to our own index
    unsigned int iGrid = xGrid->GetAtomIndex(i, j);
    unsigned int iLocal = xSiteIndex[i][j];
    
    // The medial atom to update
    MedialAtom &xAtom = xAtoms[iGrid];

    // The case where phi is negative is undefined
    if( ySolution[i][j] > 0 )
      {
      // Compute the derivatives of R using finite differences
      double xTest = xSites[iLocal]->ComputeEquation(ySolution);
      xAtom.F = xMasks[iLocal]->ComputeOneJet(ySolution, xAtom.Fu, xAtom.Fv);

      // Compute the boundary properties of the medial point
      xAtom.ComputeBoundaryAtoms(xSites[iLocal]->IsBorderSite());
      }
    else
      {
      // What are we supposed to do?
      cout << "Negative F at " << xAtom.u << ", " << xAtom.v << endl;
      xAtom.R = xAtom.F = xAtom.Fu = xAtom.Fv = 0.0;
      xAtom.flagValid = false;
      }

    // Update the valid atoms flag
    flagAllAtomsValid &= xAtom.flagValid;
    }

  // Return value: whether all atoms are valid
  return flagAllAtomsValid;
}

double ArrayMinMax(double *array, size_t n, double &xMin, double &xMax)
{
	xMax = 0, xMin = 1e100;
  for(size_t q = 0; q < n; q++)
  	{
    if(xMax < fabs(array[q])) xMax = fabs(array[q]);
    if(xMin > fabs(array[q])) xMin = fabs(array[q]);
    }
}

/** This computes the right hand side of the site equations, given phi=x */
double MedialPDESolver::ComputeNewtonRHS(const Mat& x, Mat &b)
{
  size_t i, j;
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    b[i][j] = - xSites[xSiteIndex[i][j]]->ComputeEquation(x);
  return dot_product(b, b);
}

bool MedialPDESolver::SolveOnce(const Mat &xGuess, double delta)
{
  size_t i, j, k, iIter;
  double epsMax, bMax, bMagSqr;
  
  // Initialize epsilon to zero
  eps.fill(0.0);
  
  // Copy the initial solution to the current solution
  y = xGuess;

  // Compute the right hand side, put it into b
  bMagSqr = ComputeNewtonRHS(y, b);
  // cout << " " << bMagSqr << " ";

  // We are now ready to perform the Newton loop
  for(iIter = 0; iIter < 50; iIter++)
    {
    // Compute the Jacobian matrix
    for(size_t iSite = 0; iSite < nSites; iSite++)
      xSites[iSite]->
        ComputeDerivative(y, xSparseValues + xRowIndex[iSite] - 1, iSite+1);
      
    // Perform the symbolic factorization only for the first iteration
    if(iIter == 0)
      xPardiso.SymbolicFactorization(nSites, xRowIndex, xColIndex, xSparseValues);

    // Compute the Jacobian inverse
    xPardiso.NumericFactorization(xSparseValues);
    xPardiso.Solve(b.data_block(), eps.data_block());

    // A plus means solver step
    cout << "+";

    // Advance y to the new Newton position
    y += eps; 
       
    // Compute the solution at the full Newton step
    double bMagSqrTest = ComputeNewtonRHS(y, b);

    // Perform backtracking if necessary (this is boneheaded backtracking but
    // it appears to work). TODO: Implement lnsrch from NRC
    double lambda = 0.5;
    while(bMagSqrTest > bMagSqr && lambda > 1e-4)
      {
      // Go back along eps
      y -= lambda * eps; lambda *= 0.5;

      // Compute the right hand side again
      bMagSqrTest = ComputeNewtonRHS(y, b);

      // A "-" means backtrack step
      cout << "-";
      }

    // Store the new b magnitude value
    bMagSqr = bMagSqrTest;
    // cout << " " << bMagSqr << " ";

    // Get the largest error (eps)
    epsMax = eps.array_inf_norm();
    bMax = b.array_inf_norm();

    // Print the statistics
    // cout << "-----------" << endl;
    // cout << "Step " << iIter << ": " << endl;
    // cout << "  Largest Epsilon: " << epsMax << endl;
    // cout << "  Largest Eqn Error: " << bMax << endl; 

    /*
    if(iIter > 12)
      cout << "!";
    else
      cout << "+";
      */

      // cout << endl << "Step " << iIter << ":\t" << "eMax = " 
      //   << epsMax << "\t" << "bMax = " << bMax << endl;
    
    // Test the matrix result
    // SparseLinearTest(nSites, xRowIndex, xColIndex, xSparseValues, eps, zTest, b);
    // double zMax = fabs(zTest[myidamax(nSites, zTest, 1)]);
    // cout << "  Largest Solver Error: " << zMax << endl;
   
    // Convergence is defined when epsilon is smaller than some threshold
    if(bMax < delta || epsMax < delta) 
      break;
    }

  // Let the user know if we can't converge on the root
  if(bMax > delta && epsMax > delta)
    {
    cout << "  *** CONVERGENCE FAILURE *** ";
    cout << " epsMax = " << epsMax << "; bMax = " << bMax << endl;
    return false;
    }
  else
    {
    cout << endl;
    return true;
    }

  /* 
  // Show the value of phi at the four corners
  cout << "Corner Phi Values: " << 
    y[0][0] << ", " << y[0][n-1] << ", " << y[m-1][1] << ", " << y[m-1][n-1] << endl;

  cout << "Y piece: " << endl;
  cout << y.extract(6, 6) << endl;
  cout << "B piece: " << endl;
  cout << b.extract(6, 6) << endl;

  // Estimate the LBO for 1..5
  Mat LBO(m, n, 0.0);
  for(i = 1; i < m-1; i++) for(j = 1; j < n-1; j++)
    {
    LBO[i][j] = EstimateLBOperator(y, i, j);
    }

  cout << "LBO Estimate: " << endl;
  cout << LBO.extract(16, 16) << endl;
  

  xSites[xSiteIndex[1][1]]->PrintReport();
  xMasks[xSiteIndex[1][1]]->PrintReport();
  */
}

bool
MedialPDESolver
::Solve(const Mat &xGuess, double delta)
{
  // Intialize the sites
  InitializeSiteGeometry();

  // Just solve once - no initialization tricks!
  bool flagSolved = SolveOnce(xGuess, delta);
  
  // Reconstruct the medial atoms
  bool flagValid = ReconstructAtoms(y);

  return flagValid && flagSolved;
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

  // Record the grad-phi derivative
  dAtom->xGradPhi = - T2;

  // Address the edge case first
  if(isEdge) 
    {
    dAtom->xBnd[1].X = dAtom->xBnd[0].X = N + 0.5 * T2;
    dAtom->xBnd[0].N = dAtom->xBnd[1].N = 
      0.5 * ( (T2) / xAtom->R - xAtom->xBnd[0].N * (H / F) );
    dAtom->xGradRMagSqr = 0.0;
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
  dAtom->xBnd[0].N = 0.5 * ( (T2 - T4) / xAtom->R - xAtom->xBnd[0].N * (H / F) );
  dAtom->xBnd[1].N = 0.5 * ( (T2 + T4) / xAtom->R - xAtom->xBnd[1].N * (H / F) );

  // Compute the derivative of the gradient magnitude of R
  dAtom->xGradRMagSqr = 0.25 * 
    ( 2.0 * dot_product(dAtom->xGradPhi, xAtom->xGradPhi) * F
      - dot_product(xAtom->xGradPhi, xAtom->xGradPhi) * H ) / ( F * F );
}

void
MedialPDESolver
::ComputeVariationalDerivative(IHyperSurface2D *xVariation, MedialAtom *dAtoms)
{
  size_t i, j;

  // First thing is to define the PDE's that will give us the variational
  // derivative of the function phi.

  // Compute the jet of nu at each site
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    {
    // Get the index of the site
    size_t iGrid = xGrid->GetAtomIndex(i, j);
    size_t iSite = xSiteIndex[i][j];

    // Access the medial atom underneath
    MedialAtom &dAtom = dAtoms[iGrid];

    // Set the atoms' domain coordinates
    dAtom.u = uGrid[i]; dAtom.v = vGrid[j];

    // Evaluate the variation and its derivatives
    xVariation->EvaluateAtGridIndex(i, j, 0, 0, 0, 3, dAtom.X.data_block());
    xVariation->EvaluateAtGridIndex(i, j, 1, 0, 0, 3, dAtom.Xu.data_block());
    xVariation->EvaluateAtGridIndex(i, j, 0, 1, 0, 3, dAtom.Xv.data_block());
    xVariation->EvaluateAtGridIndex(i, j, 2, 0, 0, 3, dAtom.Xuu.data_block());
    xVariation->EvaluateAtGridIndex(i, j, 1, 1, 0, 3, dAtom.Xuv.data_block());
    xVariation->EvaluateAtGridIndex(i, j, 0, 2, 0, 3, dAtom.Xvv.data_block());
    xVariation->EvaluateAtGridIndex(i, j, 0, 0, 3, 1, &dAtom.xLapR);

    // Prepare the matrices for the linear solver
    xSites[iSite]->ComputeVariationalDerivative(
        y, xSparseValues + xRowIndex[iSite] - 1, &b[i][j],
        &xAtoms[iGrid], &dAtom);
    }

  // Solve the partial differential equation
  // TODO: Figure out if we can do factorization once for the whole gradient
  // computation. That could save a significant amount of time overall
  // xPardiso.SymbolicFactorization(nSites, xRowIndex, xColIndex, xSparseValues);
  xPardiso.NumericFactorization(xSparseValues);
  xPardiso.Solve(b.data_block(), dy.data_block());

  // For each atom, compute the boundary derivatives
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    {
    // Get the index of the site
    size_t iGrid = xGrid->GetAtomIndex(i, j);
    size_t iSite = xSiteIndex[i][j];

    // Access the medial atom underneath
    MedialAtom &dAtom = dAtoms[iGrid];

    // Compute the gradient of phi for the new atom
    dAtom.F = xMasks[iSite]->ComputeOneJet(dy, dAtom.Fu, dAtom.Fv);

    // Compute the rest of the atom derivative
    ComputeMedialAtomBoundaryDerivative(
      &xAtoms[iGrid], &dAtom, xSites[iSite]->IsBorderSite());
    }
}

double fnTestF01(double u, double v)
{
  return exp(u + v);
}

void fnTestF01Jet(double u, double v, double jet[])
{
  jet[0] = exp(u + v);
  jet[1] = jet[2] = jet[3] = jet[4] = jet[5] = jet[0];
}

void MedialPDESolver
::TestFiniteDifferenceConvergence()
{
  size_t i, j;
  
  // Create a field of values Phi
  Mat F(m, n);

  // Compute the values
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    F[i][j] = fnTestF01(uGrid[i], vGrid[j]);

  Vec xMaxDiff(6, 0.0);

  // Compute the partial derivatives using all the available masks
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    {
    // The jet vectors
    Vec j1(6), j2(6);
    
    // Compute the finite difference jet
    j1[0] = xMasks[xSiteIndex[i][j]]->ComputeTwoJet(F, j1[1], j1[2], j1[3], j1[4], j1[5]);
    
    // Compute the actual derivative
    fnTestF01Jet(uGrid[i], vGrid[j], j2.data_block());

    // Compute the differences in the jets
    Vec xDiff = (j1-j2).apply(fabs);

    int kn = xSites[xSiteIndex[i][j]]->IsBorderSite() ? 3 : 6;
    for(int k = 0; k < kn; k++)
      {
      if(xMaxDiff[k] < xDiff[k])
        {
        xMaxDiff[k] = xDiff[k];
        cout << "Max changed in " << k << " to " << xDiff[k] << " at " << i << "," << j << endl;
        }
      }
    }

  cout << "MedialPDESolver FD Test: " << xMaxDiff << endl;
}
