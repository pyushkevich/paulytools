#include "MedialPDESolver.h"
#include <cmath>

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
  // we are at a border site
  i0 = iNode;
  i1 = i0 + stride_u; i2 = i0 + stride_v; 
  i3 = i0 - stride_u; i4 = i0 - stride_v;

  // Initialize the weights in u and v
  wu = 0.5; wv = 0.5;

  // Now, depending on which border we are on, correct the index and weights of finite
  // differences in u and v directions
  if(i == 0)
    { i3 = i0; wu = 1; }
  else if(i == m - 1)
    { i1 = i0; wu = 1; }
  if(j == 0)
    { i4 = i0; wv = 1; }
  else if(j == n - 1)
    { i2 = i0; wv = 1; }

  // Check whether we are at a corner
  corner = (wu == 1.0 && wv == 1.0);
  nDistinctSites = corner ? 3 : 4;

  // Determine the distinct sites (the non-corner situation)
  xDistinctSites[0] = i1; xEntry[1] = 0; 
  xDistinctSites[1] = i2; xEntry[2] = 1; 
  xDistinctSites[2] = i3; xEntry[3] = 2;
  xDistinctSites[3] = i4; xEntry[4] = 3;

  // Correct the above arrays in case of corners   
  if(i1 == i2) 
    { 
    xDistinctSites[xEntry[2]] = i4; 
    xEntry[4] = xEntry[2]; 
    xEntry[2] = xEntry[1]; 
    }
  else if(i2 == i3)
    { 
    xDistinctSites[xEntry[3]] = i4; 
    xEntry[4] = xEntry[3]; 
    xEntry[3] = xEntry[2]; 
    }
  else if(i3 == i4)
    { 
    xEntry[4] = xEntry[3]; 
    }
  else if(i4 == i1)
    { 
    xEntry[4] = xEntry[1]; 
    }

  // Find the entry corresponding to the center index
  for(unsigned int j = 0; j < nDistinctSites; j++)
    if(xDistinctSites[j] == i0)
      xEntry[0] = j;
}

/** Compute the equation corresponding to an internal site. */
double FDInternalSite::ComputeEquation(Vector *Y)
{
  // Get the components
  double 
    F0 = V_GetCmp(Y, i0 + 1), F1 = V_GetCmp(Y, i1 + 1), F2 = V_GetCmp(Y, i2 + 1), 
    F3 = V_GetCmp(Y, i3 + 1), F4 = V_GetCmp(Y, i4 + 1), F5 = V_GetCmp(Y, i5 + 1), 
    F6 = V_GetCmp(Y, i6 + 1), F7 = V_GetCmp(Y, i7 + 1), F8 = V_GetCmp(Y, i8 + 1);
  
  // Compute the finite differences, without scaling factors
  double Fu = F1 - F5;
  double Fv = F3 - F7;
  double Fuu = F1 + F5 - F0 - F0;
  double Fvv = F3 + F7 - F0 - F0;
  double Fuv = F2 + F6 - F4 - F8;

  // Compute the generalized laplacian using precomputed coefficients
  return Cuu * Fuu + Cuv * Fuv + Cvv * Fvv + Cu * Fu + Cv * Fv - rho;
}

double FDBorderSite::ComputeEquation(Vector *Y)
{
  // Compute the finite differences, without scaling factors
  double Fu = V_GetCmp(Y, i1 + 1) - V_GetCmp(Y, i3 + 1);
  double Fv = V_GetCmp(Y, i2 + 1) - V_GetCmp(Y, i4 + 1);

  // Compute the generalized gradient magnitude using precomputed values
  return CuCu * Fu * Fu + CuCv * Fu * Fv + CvCv * Fv * Fv - C0 * V_GetCmp(Y, i0 + 1);
}

/** Compute the derivatives for Newton's method */
void FDInternalSite::ComputeDerivative(Vector *, QMatrix *A, unsigned int iRow)
{
  // Simply copy the values of the b's into A
  Q_SetEntry(A, iRow, 0, i0 + 1, b0);
  Q_SetEntry(A, iRow, 1, i1 + 1, b1);
  Q_SetEntry(A, iRow, 2, i2 + 1, b2);
  Q_SetEntry(A, iRow, 3, i3 + 1, b3);
  Q_SetEntry(A, iRow, 4, i4 + 1, b4);
  Q_SetEntry(A, iRow, 5, i5 + 1, b5);
  Q_SetEntry(A, iRow, 6, i6 + 1, b6);
  Q_SetEntry(A, iRow, 7, i7 + 1, b7);
  Q_SetEntry(A, iRow, 8, i8 + 1, b8);
}

/** Compute the derivatives for Newton's method */
void FDBorderSite::ComputeDerivative(Vector *Y, QMatrix *A, unsigned int iRow)
{
  // This computes the derivative of the site's equation with respect to the
  // four variables involved in the equation.
  double Fu = V_GetCmp(Y, i1 + 1) - V_GetCmp(Y, i3 + 1);
  double Fv = V_GetCmp(Y, i2 + 1) - V_GetCmp(Y, i4 + 1);
  double d1 = (CuCu + CuCu) * Fu + CuCv * Fv;
  double d2 = CuCv * Fu + (CvCv + CvCv) * Fv;

  // Clear all the entries to zero
  for(unsigned int j = 0; j < nDistinctSites; j++)
    Q_SetEntry(A, iRow, j, xDistinctSites[j] + 1, 0.0);

  // Compute the derivatives
  Q_AddVal(A, iRow, xEntry[0], -C0 );
  Q_AddVal(A, iRow, xEntry[1], d1 );
  Q_AddVal(A, iRow, xEntry[2], d2 );
  Q_AddVal(A, iRow, xEntry[3], -d1 );
  Q_AddVal(A, iRow, xEntry[4], -d2 );
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

static char *dummy = "Newton";

MedialPDESolver
::MedialPDESolver(unsigned int m, unsigned int n)
{
  // Copy the site dimensions
  this->m = m; this->n = n;
  du = 1.0 / (m - 1); dv = 1.0 / (n - 1);

  nSites = m * n;
  xSites = new FDAbstractSite*[nSites];
  xSiteIndex = new unsigned int*[m];

  // Initialize the A Matrix and the b vector storage
  Q_Constr(&A,"Newton Mat",nSites,False,Rowws,Normal,True);
  V_Constr(&b,"Newton RHS",nSites,Normal,True);
  V_Constr(&y,"Newton Sol",nSites,Normal,True);
  V_Constr(&eps,"Newton Epsilon",nSites,Normal,True);

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

      // Initialize the row in the Q Matrix
      Q_SetLen(&A, iSite + 1, xSites[iSite]->GetNumberOfNetwonColumns());

      // Associate the raw offset with the u-v index
      xSiteIndex[i][j] = iSite++;
      }
    }

  // Compute the initial solution
  xInitSoln = new double[nSites];
  for(iSite = 0; iSite < nSites; iSite++)
    xInitSoln[iSite] = xSites[iSite]->GetInitialValue();
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

      // Copy the initial solution to the current solution
      V_SetCmp(&y, iSite+1, xInitSoln[iSite]);

      // Compute the surface jet and the laplacian
      problem->ComputeJet2(u,v,X,Xu,Xv,Xuu,Xuv,Xvv);

      // Compute the geometrical properties
      GeometryDescriptor gd(X,Xu,Xv,Xuu,Xuv,Xvv);

      // Compute the solution at this point
      xSites[iSite]->SetGeometry(&gd, problem->ComputeLaplacian(u,v));
      }
    }

  // Initialize epsilon to zero
  V_SetAllCmp(&eps, 0.0);

  // Set the accuracy of the solver
  SetRTCAccuracy(1e-8);

  // We are now ready to perform the Newton loop
  bool flagComplete = false;
  for(unsigned int iIter = 0; iIter < 16 && !flagComplete; iIter++)
    {
    // Compute the A matrix and the b vector
    for(k = 0; k < nSites; k++)
      {
      xSites[k]->ComputeDerivative(&y, &A, k+1);
      V_SetCmp(&b, k+1, - xSites[k]->ComputeEquation(&y) );
      }

    // Solve the equation for X
    BiCGIter(&A, &eps, &b, nSites, NULL, 1.0);

    /*
    cout << "A = "; DumpQMatrix(&A);
    cout << "b = "; DumpVector(&b);
    cout << "eps = "; DumpVector(&eps);
    */

    // Get the largest error (eps)
    double epsMax = MaxAbs_VV(&eps);
    double bMax = MaxAbs_VV(&b);

    // Append the epsilon vector to the result
    AddAsgn_VV(&y, &eps);
    
    // Print the statistics
    cout << "-----------" << endl;
    cout << "Step " << iIter << ": " << endl;
    cout << "  Lin. Solver Error: " << MaxAbs_VV( Sub_VV(Mul_QV(&A, &eps), &b) ) << endl;
    cout << "  Largest Epsilon: " << epsMax << endl;
    cout << "  Largest Eqn Error: " << bMax << endl;

    // Convergence is defined when epsilon is smaller than some threshold
    flagComplete = (epsMax < delta);
    }
}
