#include "FourierSurface.h"
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>

using namespace std;

// Lapack definitions
extern "C" {
  void sgetrf(int *m,int *n,float *a,int *lda,int *ipiv,int *info);
  void sgetrs(char *trans,int *n,int *nrhs,float *a,int *lda,int *ipiv,float *b,int *ldb,int *info);
}

FourierSurface
::FourierSurface(unsigned int nBasesU, unsigned int nBasesV, 
  unsigned int nComponents)
{
  // Store the parameters
  this->nBasesU = nBasesU;
  this->nBasesV = nBasesV;
  this->nComponents = nComponents;

  // Initialize the coefficient array to zero
  xCoeff = new float **[nBasesU];
  for(unsigned int iu = 0; iu < nBasesU; iu++)
    {
    xCoeff[iu] = new float *[nBasesV];
    for(unsigned int iv = 0; iv < nBasesV; iv++)
      {
      xCoeff[iu][iv] = new float[nComponents];
      memset(xCoeff[iu][iv], 0, sizeof(float) * nComponents);
      }
    }

  // Allocate the temp arrays
  xBasisTempU = new float[nBasesU];
  xBasisTempV = new float[nBasesV];
}

float FourierSurface
::Evaluate(float u, float v, 
  unsigned int iDerivativeU, unsigned int iDerivativeV, 
  unsigned int iComponent)
{
  /** Compute the basis values */
  for(unsigned int i = 0; i < nBasesU; i++) 
    xBasisTempU[i] = BasisFunction(i, iDerivativeU, u);

  for(unsigned int j = 0; j < nBasesV; j++) 
    xBasisTempV[j] = BasisFunction(j, iDerivativeV, v);
  
  /** Sum over the orders - time consuming! */
  float xValue = 0.0;
  for(unsigned int i = 0; i < nBasesU; i++) 
    for(unsigned int j = 0; j < nBasesV; j++) 
      xValue += xCoeff[i][j][iComponent] * xBasisTempU[i] * xBasisTempV[j];

  return xValue;      
}

float FourierSurface
::BasisFunction(unsigned int iOrder, unsigned int iDerivative, float uValue)
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

void FourierSurface
::FitData(unsigned int iComponent, unsigned int nPoints, 
  float *xPointU, unsigned int iStrideU, 
  float *xPointV, unsigned int iStrideV,  
  float *xValues, unsigned int iStrideX,
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
  float **Z = new float*[nUnkowns];
  for(iu = 0, i = 0; iu < nu; iu++) for(iv = 0; iv < nv; iv++, i++)
    {
    Z[i] = new float[nPoints];
    
    float *pu = xPointU, *pv = xPointV;
    for(k = 0; k < nPoints; k++, pu += iStrideU, pv += iStrideV) 
      // Z[i][k] = BasisFunction(iu, 0, *pu) * BasisFunction(iv, 0, *pv); 
      Z[i][k] = BasisFunction(iu, 0, xPointU[k]) 
        * BasisFunction(iv, 0, xPointV[k]); 
    }

  // Allocate the matrix A and vector b
  float *A = new float[nUnkowns * nUnkowns];
  float *b = new float[nUnkowns];

  // Set the elements of A and b
  unsigned int offset = 0;
  for(j = 0; j < nUnkowns; j++)
    {
    // Compute the b vector
    b[j] = 0;

    float *px = xValues;
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
  sgetrf( &nRows, &nRows, A, &nRows, iPivot, &iInfo);  
  if(iInfo < 0)
    { cerr << "Error calling sgetrf" << endl; return; }
  else if(iInfo > 0)
    { cerr << "sgetrf: Matrix is singular" << endl; return; }

  // Solve the system
  char cTrans = 'N';
  int nRhs = 1; 
  sgetrs(&cTrans, &nRows, &nRhs, A, &nRows, iPivot, b, &nRows, &iInfo);
  if(iInfo < 0)
    { cerr << "Error calling sgetrs" << endl; return; }

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
