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

template<class TSurface, class TLaplacian>
MedialPDESolver<TSurface,TLaplacian>
::MedialPDESolver(unsigned int m, unsigned int n)
{
  // Copy the site dimensions
  this->m = m; this->n = n;
  du = 1.0 / (m - 1); dv = 1.0 / (n - 1);

  nSites = m * n;
  xSites = new FDAbstractSite*[nSites];
  xSiteIndex = new unsigned int*[m];

  // Build a site index
  unsigned int iSite = 0, i, j;
  for(i = 0; i < m; i++) for(j = 0; j < n; j++) xSiteIndex[i][j] = iSite++;

  // Create the internal sites
  for(i = 1; i < m-1; i++)
    for(j = 1; j < n-1; j++)
      {
      xSites[xSiteIndex[i][j]] = new FDInternalSite(m, n, i, j);
      }

  // Create the border sites
  for(i = 1; i < m - 1; i++)
    {
    xSites[xSiteIndex[i][0]] = new FDBorderSite(m, n, i, 0);
    xSites[xSiteIndex[i][n-1]] = new FDBorderSite(m, n, i, n-1);
    }
  for(j = 1; j < n - 1; j++)
    {
    xSites[xSiteIndex[0][j]] = new FDBorderSite(m, n, 0, j);
    xSites[xSiteIndex[m-1][j]] = new FDBorderSite(m, n, m-1, j);
    }

  // Create the corner sites
  xSites[xSiteIndex[0][0]] = new FDCornerSite(m, n, 0, 0);
  xSites[xSiteIndex[0][n-1]] = new FDCornerSite(m, n, 0, n-1);
  xSites[xSiteIndex[m-1][0]] = new FDCornerSite(m, n, m-1, 0);
  xSites[xSiteIndex[m-1][n-1]] = new FDCornerSite(m, n, m-1, n-1);

  // Compute the initial solution
  xInitSoln = new double[nSites];
  for(iSite = 0; iSite < nSites; iSite++)
    xInitSoln = xSites[iSite]->GetInitialValue();
}

template<class TSurface, class TLaplacian>
void
MedialPDESolver<TSurface,TLaplacian>
::Solve(TSurface &surface, TLaplacian &lap)
{
  // The Jet of a surface
  double X[3], Xu[3], Xv[3], Xuu[3], Xuv[3], Xvv[3];
  double u = 0.0, v = 0.0;

  // Compute the jet at every point on the surface
  for(unsigned int i = 0; i < m; i++, u+=du)
    for(unsigned int j = 0; j < n; j++, v+=dv)
      {
      // Compute the surface jet
      surface.ComputeJet2(u,v,X,Xu,Xv,Xuu,Xuv,Xvv);

      // Compute the geometrical properties
      GeometryDescriptor gd(X,Xu,Xv,Xuu,Xuv,Xvv);

      // Compute the solution at this point
      xSites[xSiteIndex[i][j]]->ComputeEquation();
      }
  

}

/** Initialize the site with grid of size m and n */
FDInternalSite::FDInternalSite(unsigned int m, unsigned int n, unsigned int i, unsigned int j)
{
  // Set the indices of the eight neighbors in the array, go around a circle
  i0 = m * j + i;
  i1 = i0 + 1; i3 = i0 + m; i5 = i0 - 1; i7 = i0 - m;
  i2 = i3 + 1; i4 = i3 - 1; i6 = i7 - 1; i8 = i7 + 1;
}

/** Compute the equation corresponding to an internal site. */
double FDInternalSite::ComputeEquation(double *F)
{
  // Compute the finite differences, without scaling factors
  F0 = F[i0];
  Fu = F[i1] - F[i5];
  Fv = F[i3] - F[i7];
  Fuu = F[i1] + F[i5] - F0 - F0;
  Fvv = F[i3] - F[i7] - F0 - F0;
  Fuv = F[i2] + F[i6] - F[i4] - F[i8];

  // Compute the generalized laplacian using precomputed coefficients
  return Cuu * Fuu + Cuv * Fuv + Cvv * Fvv + Cu * Fu + Cv * Fv - rho;
}

/** Compute the derivatives for Newton's method */
void FDInternalSite::ComputeDerivative(double *X, double *A)
{
  // Simply copy the values of the b's into A
  A[i0] = b0; 
  A[i1] = b1; A[i2] = b2; A[i3] = b3; A[i4] = b4;
  A[i5] = b5; A[i6] = b6; A[i7] = b7; A[i8] = b8; 
}

/** Initialize internal site */
template<class TSurface, class TLaplacian>
void FDInternalSite::SetGeometry(GeometryDescriptor *g, double rho)
{
  // Store the rho
  this->rho = rho;

  // Compute the finite difference premultipliers
  double _du = 1.0 / (m-1);
  double _dv = 1.0 / (n-1);
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
  b0 = -2.0 * (Cuu + Cuv);
  b1 = Cuu + Cu;
  b2 = Cuv;
  b3 = Cvv + Cv;
  b4 = -Cuv;
  b5 = Cuu - Cu;
  b6 = Cuv;
  b7 = Cvv - Cv;
  b8 = -Cuv;
}
