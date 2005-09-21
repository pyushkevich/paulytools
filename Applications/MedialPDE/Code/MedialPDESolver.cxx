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
        // cout << uGrid[i+1] - uGrid[i] << " <> " << uGrid[i] - uGrid[i-1] << endl;
        // cout << vGrid[j+1] - vGrid[j] << " <> " << vGrid[j] - vGrid[j-1] << endl;
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
    double uMin = i > m/2 ? uGrid[m-(i+1)] : uGrid[i];
    double vMin = j > n/2 ? vGrid[n-(j+1)] : vGrid[j];
    
    // xDefaultInitSoln[i][j] = 0.1 + uMin * uMin + vMin * vMin;
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
  
  // Generate a matrix of x, y, z values
  Mat vx(m, n), vy(m, n), vz(m, n);

  // Sample the surface along the current grid
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    {
    double X[3];
    xSurface->EvaluateAtGridIndex(i, j, 0, 0, 0, 3, X);
    vx[i][j] = X[0]; vy[i][j] = X[1]; vz[i][j] = X[2];
    }

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
    
    // Compute the surface jet and the laplacian
    /* 
    xSurface->EvaluateAtGridIndex(i, j, 0, 0, 0, 3, xAtom.X.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 1, 0, 0, 3, xAtom.Xu.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 0, 1, 0, 3, xAtom.Xv.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 2, 0, 0, 3, xAtom.Xuu.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 1, 1, 0, 3, xAtom.Xuv.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 0, 2, 0, 3, xAtom.Xvv.data_block());
    */
    
    xAtom.X[0] = xMasks[iSite]->ComputeTwoJet(vx, xAtom.Xu[0], xAtom.Xv[0], xAtom.Xuu[0], xAtom.Xuv[0], xAtom.Xvv[0]);
    xAtom.X[1] = xMasks[iSite]->ComputeTwoJet(vy, xAtom.Xu[1], xAtom.Xv[1], xAtom.Xuu[1], xAtom.Xuv[1], xAtom.Xvv[1]);
    xAtom.X[2] = xMasks[iSite]->ComputeTwoJet(vz, xAtom.Xu[2], xAtom.Xv[2], xAtom.Xuu[2], xAtom.Xuv[2], xAtom.Xvv[2]);

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
::ReconstructAtoms(const Mat &ySolution)
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
    if( ySolution[i][j] > 0 )
      {
      // Compute the derivatives of R using finite differences
      double xTest = xSites[iLocal]->ComputeEquation(ySolution);
      xAtom.F = xMasks[iLocal]->ComputeOneJet(ySolution, xAtom.Fu, xAtom.Fv);

      // Compute the boundary properties of the medial point
      xAtom.ComputeBoundaryAtoms();
      }
    else
      {
      // What are we supposed to do?
      cout << "Negative F at " << xAtom.u << ", " << xAtom.v << endl;
      xAtom.F = xAtom.Fu = xAtom.Fv = 0.0;
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
  for(unsigned int iIter = 0; iIter < 32 && !flagComplete; iIter++)
    {
    // Compute the A matrix and the b vector at each node
    for(i = 0; i < m; i++) for(j = 0; j < n; j++)
      {
      // Get the site index
      size_t iSite = xSiteIndex[i][j];
      
      // Compute the non-zero values of A for this row
      xSites[iSite]->
        ComputeDerivative(y, xSparseValues + xRowIndex[iSite] - 1, iSite+1);

      // Compute the value of b
      b[i][j] = -xSites[iSite]->ComputeEquation(y);
      }
      
    // Report the largest and smallest values in A and b
    Vec xSparse(xSparseValues, xRowIndex[nSites]-1);
    double xMax = 0, xMin = 1e100;
    for(size_t q = 0; q < xSparse.size(); q++)
      {
      if(xMax < fabs(xSparse[q])) xMax = fabs(xSparse[q]);
      if(xMin > fabs(xSparse[q])) xMin = fabs(xSparse[q]);
      }
    
    cout << "max(A) = " << xMax << "; min(A) = " << xMin << endl;

    // Perform the symbolic factorization only for the first iteration
    if(iIter == 0)
      xPardiso.SymbolicFactorization(nSites, xRowIndex, xColIndex, xSparseValues);

    // Peform numeric factorization and solution at each step
    xPardiso.NumericFactorization(xSparseValues);
    xPardiso.Solve(b.data_block(), eps.data_block());

    // Get the largest error (eps)
    epsMax = eps.array_inf_norm();
    bMax = b.array_inf_norm();

    // Append the epsilon vector to the result
    y += eps;

    // Print the statistics
    
    /* 
    cout << "-----------" << endl;
    cout << "Step " << iIter << ": " << endl;
    cout << "  Largest Epsilon: " << epsMax << endl;
    cout << "  Largest Eqn Error: " << bMax << endl; 
    */
    cout << "Step " << iIter << ":\t" << "eMax = " 
      << epsMax << "\t" << "bMax = " << bMax << endl;
    
    // Test the matrix result
    // SparseLinearTest(nSites, xRowIndex, xColIndex, xSparseValues, eps, zTest, b);
    // double zMax = fabs(zTest[myidamax(nSites, zTest, 1)]);
    // cout << "  Largest Solver Error: " << zMax << endl;
   

    // Convergence is defined when epsilon is smaller than some threshold
    if(bMax < delta) 
      break;
    }

  // Show the value of phi at the four corners
  cout << "Corner Phi Values: " << 
    y[0][0] << ", " << y[0][n-1] << ", " << y[m-1][1] << ", " << y[m-1][n-1] << endl;

  cout << "Y piece: " << endl;
  cout << y.extract(6, 6) << endl;
  cout << "B piece: " << endl;
  cout << b.extract(6, 6) << endl;

  return bMax;
}

void
MedialPDESolver
::Solve(double delta)
{
  /* 
  unsigned int i, j, k;
  Mat xBestInit(m, n);
  */
  
  // Intialize the sites
  InitializeSiteGeometry();

  // Just solve once - to initialization tricks!
  SolveOnce(delta);
  
  /* 

  // Try solving with the current initial guess
  double xBestGuess = SolveOnce(delta);
  if(false && xBestGuess > delta)
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
  */

  // Reconstruct the medial atoms
  ReconstructAtoms(y);
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

    // Split depending on the type of variation
    if(flagRhoVariation)
      {
      // Compute the rho derivative with respect to the variation
      xVariation->EvaluateAtGridIndex(i, j, 0, 0, 3, 1, &dAtom.xLapR);

      // Prepare the matrices for the linear solver
      xSites[iSite]->ComputeVariationalDerivativeRho(
        y, xSparseValues + xRowIndex[iSite] - 1, &b[i][j],
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
        y, xSparseValues + xRowIndex[iSite] - 1, &b[i][j],
        &xAtoms[iGrid], &dAtom);
      }
    }

  // Solve the partial differential equation
  // TODO: Figure out if we can do factorization once for the whole gradient
  // computation. That could save a significant amount of time overall
  xPardiso.SymbolicFactorization(nSites, xRowIndex, xColIndex, xSparseValues);
  xPardiso.NumericFactorization(xSparseValues);
  xPardiso.Solve(b.data_block(), eps.data_block());

  // For each atom, compute the boundary derivatives
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    {
    // Get the index of the site
    size_t iGrid = xGrid->GetAtomIndex(i, j);
    size_t iSite = xSiteIndex[i][j];

    // Access the medial atom underneath
    MedialAtom &dAtom = dAtoms[iGrid];

    // Compute the gradient of phi for the new atom
    dAtom.F = xMasks[iSite]->ComputeOneJet(eps, dAtom.Fu, dAtom.Fv);

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