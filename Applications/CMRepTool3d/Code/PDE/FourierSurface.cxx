#include "FourierSurface.h"
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>

using namespace std;

// Lapack definitions
extern "C" {
  void dgetrf(int *m,int *n,double *a,int *lda,int *ipiv,int *info);
  void dgetrs(char *trans,int *n,int *nrhs,double *a,int *lda,int *ipiv,double *b,int *ldb,int *info);
}

FourierSurfaceOld
::FourierSurfaceOld(unsigned int nBasesU, unsigned int nBasesV, 
  unsigned int nComponents)
{
  // Store the parameters
  this->nBasesU = nBasesU;
  this->nBasesV = nBasesV;
  this->nComponents = nComponents;

  // Initialize the coefficient array to zero
  xCoeff = new double **[nBasesU];
  for(unsigned int iu = 0; iu < nBasesU; iu++)
    {
    xCoeff[iu] = new double *[nBasesV];
    for(unsigned int iv = 0; iv < nBasesV; iv++)
      {
      xCoeff[iu][iv] = new double[nComponents];
      memset(xCoeff[iu][iv], 0, sizeof(double) * nComponents);
      }
    }

  // Allocate the temp arrays
  xBasisTempU = new double[nBasesU];
  xBasisTempV = new double[nBasesV];
}

double FourierSurfaceOld
::Evaluate(double u, double v, 
  unsigned int iDerivativeU, unsigned int iDerivativeV, 
  unsigned int iComponent)
{
  /** Compute the basis values */
  for(unsigned int i = 0; i < nBasesU; i++) 
    xBasisTempU[i] = BasisFunction(i, iDerivativeU, u);

  for(unsigned int j = 0; j < nBasesV; j++) 
    xBasisTempV[j] = BasisFunction(j, iDerivativeV, v);
  
  /** Sum over the orders - time consuming! */
  double xValue = 0.0;
  for(unsigned int i = 0; i < nBasesU; i++) 
    for(unsigned int j = 0; j < nBasesV; j++) 
      xValue += xCoeff[i][j][iComponent] * xBasisTempU[i] * xBasisTempV[j];

  return xValue;      
}

double FourierSurfaceOld
::BasisFunction(unsigned int iOrder, unsigned int iDerivative, double uValue)
{
  // Compute the basis at the point
  double uu = 0.8 * uValue + 0.1;
  if(iDerivative == 0)
    return cos(M_PI * iOrder * uu);
  else if(iDerivative == 1)
    return - 0.8 * M_PI * iOrder * sin(M_PI * iOrder * uu);
  else 
    return - 0.8 * 0.8 * M_PI * M_PI * iOrder * iOrder * cos(M_PI * iOrder * uu);
}

void FourierSurfaceOld
::FitData(unsigned int iComponent, unsigned int nPoints, 
  double *xPointU, unsigned int iStrideU, 
  double *xPointV, unsigned int iStrideV,  
  double *xValues, unsigned int iStrideX,
  unsigned int nu, unsigned int nv)
{
  // Correct the number of bases
  nu = nu == 0 ? nBasesU : nu;
  nv = nv == 0 ? nBasesV : nv;

  // Get the number of unknown coefficients
  unsigned int iu, iv, i = 0, j, k;
  unsigned int nUnkowns = nu * nv;

  cout << "Fitting data with " << nUnkowns << " unknowns to " 
    << nPoints << " points." << endl;

  // Compute the basis for each point
  double **Z = new double*[nUnkowns];
  for(iu = 0, i = 0; iu < nu; iu++) for(iv = 0; iv < nv; iv++, i++)
    {
    Z[i] = new double[nPoints];
    
    double *pu = xPointU, *pv = xPointV;
    for(k = 0; k < nPoints; k++, pu += iStrideU, pv += iStrideV) 
      // Z[i][k] = BasisFunction(iu, 0, *pu) * BasisFunction(iv, 0, *pv); 
      Z[i][k] = BasisFunction(iu, 0, xPointU[k]) 
        * BasisFunction(iv, 0, xPointV[k]); 
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

    double *px = xValues;
    for(k = 0; k < nPoints; k++, px += iStrideX) 
      // b[j] += (*px) * Z[j][k];
      b[j] += xValues[k] * Z[j][k]; 

    // Compute the elements of the A matrix
    for(i = 0; i < nUnkowns; i++)
      {
      A[offset] = 0;
      for(k = 0; k < nPoints; k++) 
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

  // Clear the coefficients
  for(iu = 0; iu < nBasesU; iu++) for(iv = 0; iv < nBasesV; iv++)
    xCoeff[iu][iv][iComponent] = 0;

  // Solution has been placed into B. Map it to the coefficients
  i = 0;
  for(iu = 0; iu < nu; iu++) for(iv = 0; iv < nv; iv++)
    xCoeff[iu][iv][iComponent] = b[i++];

  // Clean up the resources
  delete b;
  delete A;
  for(i = 0; i < nUnkowns; i++)
    { delete Z[i]; }
}
