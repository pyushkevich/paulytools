// Lapack definitions
extern "C" {
  void dgetrf(int *m,int *n,double *a,int *lda,int *ipiv,int *info);
  void dgetrs(char *trans,int *n,int *nrhs,double *a,int *lda,int *ipiv,double *b,int *ldb,int *info);
}

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::GenericBasisRepresentation2D(size_t ncu, size_t ncv)
: C(NComponents, ncu, ncv)
{
  Initialize(ncu, ncv);
}

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::Initialize(size_t ncu, size_t ncv)
{
  // initialize the coefficient array
  this->ncu = ncu; this->ncv = ncv;
  nc = ncu * ncv;
  ncRaw= nc * NComponents;

  // initialize the evaluation arrays
  uEval.resize(ncu); vEval.resize(ncv);
}


template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::EvaluateDerivative(double u, double v, size_t ou, size_t ov, double *x)
{
  // No grid is available - must evaluate each basis function
  size_t iu, iv, k;
  for(iu = 0; iu < ncu; iu++) uEval[iu] = fu.Evaluate(u, iu, ou);
  for(iv = 0; iv < ncv; iv++) vEval[iv] = fv.Evaluate(v, iv, ov);

  // Clear the output array
  for(k = 0; k < NComponents; k++) x[k] = 0.0;

  // Compute the cross product of the bases
  unsigned int iOffset = 0;
  for(iv = 0; iv < ncv; iv++) for(iu = 0; iu < ncu; iu++) 
    {
    // Compute the product of the basis functions
    double Fuv = uEval[iu] * vEval[iv];

    // Get the coefficient corresponding to iu and iv
    for(k = 0; k < NComponents; k++)
      x[k] += C[iOffset++] * Fuv;
    }
}

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::SetEvaluationGrid(size_t nu, size_t nv, double *uu, double *vv)
{
  // Initialize the evaluation grid to the specified dimensions
  uGrid.resize(nu * ncu * NOrder);
  vGrid.resize(nv * ncv * NOrder);
  uGridValues.resize(nu); 
  vGridValues.resize(nv);

  // Precompute the U-grid
  size_t iGrid = 0;
  for(size_t iu = 0; iu < nu; iu++)
    {
    uGridValues[iu] = uu[iu];
    for(size_t icu = 0; icu < ncu; icu++)
      for(size_t p = 0; p < NOrder; p++)
        uGrid[iGrid++] = fu.Evaluate(uGridValues[iu], icu, p);
    }
  
  // Precompute the V-grid
  iGrid = 0;
  for(size_t iv = 0; iv < nv; iv++)
    {
    vGridValues[iv] = vv[iv];
    for(size_t icv = 0; icv < ncv; icv++)
      for(size_t p = 0; p < NOrder; p++)
        vGrid[iGrid++] = fv.Evaluate(vGridValues[iv], icv, p);
    }
} 

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::EvaluateAtGridIndex(size_t iu, size_t iv, size_t ou, size_t ov, double *x)
{
  // Get a slice into the U grid. Indexing over this slice gets consecutive
  // basis functions of order ou at site iu
  // slice sGridU = 
  valarray<double> sGridU = uGrid[slice(NOrder * ncu * iu + ou, ncu, NOrder)]; 
  valarray<double> sGridV = vGrid[slice(NOrder * ncv * iv + ov, ncv, NOrder)];

  // Clear the output array
  size_t icu, icv, k;
  for(k = 0; k < NComponents; k++) x[k] = 0.0;

  // Compute the cross product of the bases
  size_t iOffset = 0;
  for(icv = 0; icv < ncv; icv++) for(icu = 0; icu < ncu; icu++) 
    {
    // Compute the product of the basis functions
    double Fuv = sGridU[icu] * sGridV[icv];

    // Get the coefficient corresponding to iu and iv
    for(k = 0; k < NComponents; k++)
      x[k] += C[iOffset++] * Fuv;
    }
}

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::FitToData(size_t n, size_t iComponent, double *uu, double *vv, double *xx)
{
  // The number of coefficients may not exceed the number of data points being
  // fitted (otherwise the fit will be nonsensical)
  size_t nu = ncu; size_t nv = ncv;
  while(nu * nv > n)
    { nu--; nv--; }

  // Get the number of unknown coefficients
  size_t iu, iv, i = 0, j, k;
  size_t nUnkowns = nu * nv;

  cout << "Fitting data with " << nUnkowns << " unknowns to " << n << " points." << endl;

  // Compute the basis for each point
  double **Z = new double*[nUnkowns];
  for(iu = 0, i = 0; iu < nu; iu++) for(iv = 0; iv < nv; iv++, i++)
    {
    Z[i] = new double[n];
    for(k = 0; k < n; k++) 
      Z[i][k] = fu.Evaluate(uu[k], iu, 0) * fv.Evaluate(vv[k], iv, 0); 
    }

  // Allocate the matrix A and vector b
  double *A = new double[nUnkowns * nUnkowns];
  double *b = new double[nUnkowns];

  // Set the elements of A and b
  unsigned int offset = 0;
  for(j = 0; j < nUnkowns; j++)
    {
    // Compute the b vector
    b[j] = 0;
    for(k = 0; k < n; k++) 
      b[j] += xx[k] * Z[j][k]; 

    // Compute the elements of the A matrix
    for(i = 0; i < nUnkowns; i++)
      {
      A[offset] = 0;
      for(k = 0; k < n; k++) 
        A[offset] += Z[i][k] * Z[j][k];
      ++offset;
      }
    }

  // Solve the system Ax = b (LU decomposition)
  int *iPivot = new int[nUnkowns], iInfo, nRows = nUnkowns;
  dgetrf( &nRows, &nRows, A, &nRows, iPivot, &iInfo);  
  if(iInfo < 0)
    { cerr << "Error calling dgetrf" << endl; return; }
  else if(iInfo > 0)
    { cerr << "dgetrf: Matrix is singular" << endl; return; }

  // Solve the system
  char cTrans = 'N';
  int nRhs = 1; 
  dgetrs(&cTrans, &nRows, &nRhs, A, &nRows, iPivot, b, &nRows, &iInfo);
  if(iInfo < 0)
    { cerr << "Error calling dgetrs" << endl; return; }

  // Solution has been placed into B. Map it to the coefficients
  i = 0;
  for(iu = 0; iu < nu; iu++) for(iv = 0; iv < nv; iv++)
    C(iComponent, iu, iv) = b[i++];

  // Clean up the resources
  delete b;
  delete A;
  for(i = 0; i < nUnkowns; i++)
    { delete Z[i]; }
}


template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::SaveToRegistry(Registry &R)
{
  // Store the main information
  R["Size.U"] << ncu;
  R["Size.V"] << ncv;
  R["Components"] << NComponents;

  // Store each coefficient
  for(size_t i = 0; i < ncu; i++) for(size_t j = 0; j < ncv; j++)
      for(size_t k = 0; k < NComponents; k++)
        R [ R.Key("Coefficient[%d][%d][%d]", i, j, k) ] << C(k, i, j);
}

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
bool
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::ReadFromRegistry(Registry &R)
{
  size_t mu = R["Size.U"][0];
  size_t mv = R["Size.V"][0];
  size_t mk = R["Components"][0];

  // Check the parameters
  if(mu == 0 || mv == 0 || mk != NComponents) return false;

  // Reinitialize the coefficient array
  C = Index3D(NComponents, mu, mv);
  Initialize(mu, mv);
  
  // Load the coefficients
  for(size_t i = 0; i < mu; i++) for(size_t j = 0; j < mv; j++)
    for(size_t k = 0; k < NComponents; k++)
      C(k, i, j) = R[ R.Key("Coefficient[%d][%d][%d]", i, j, k) ][0.0];
}

