#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>

#include "BasisFunctions2D.h"
#include "BasisFunctions2D.txx"
#include "CoefficientMask.h"

double CosineBasisFunction::Evaluate(double u, size_t k, size_t d)
{
  // Shift u so we don't have zero derivatives
  double b = 0.23525435234532, a = 1.0 - 2.0 * b;

  // Compute the basis at the point
  double T1 = M_PI * k * a, T2 = M_PI * k * (a * u + b);

  switch(d) 
    {
  case 0: return cos(T2);
  case 1: return - T1 * sin(T2);
  case 2: return - T1 * T1 * cos(T2);
  default: return 0.0;
    }
}

FourierSurface::FourierSurface(const FourierSurface &source)
  : FourierSurfaceBase(source.ncu, source.ncv)
{
  // Copy the raw coefficients
  this->SetCoefficientArray(source.GetCoefficientArray());
}

vnl_vector<double> FourierSurface::GetCenterOfRotation()
{
  vnl_vector<double> xCtr(4);
  xCtr[0] = C(0, 0, 0);
  xCtr[1] = C(1, 0, 0);
  xCtr[2] = C(2, 0, 0);
  xCtr[3] = C(3, 0, 0);
  return xCtr;
}

void FourierSurface::ApplyAffineTransform(const vnl_matrix<double> &A, 
  const vnl_vector<double> &b, const vnl_vector<double> &c)
{
  // Extract the upper 3x3 matrix from the input
  vnl_matrix<double> A3 = A.extract(3, 3);
  vnl_vector<double> B3 = b.extract(3);
  vnl_vector<double> C3 = c.extract(3);

  // If the matrix is 4 x 4, the last component is the rho scaling
  double xRhoScale = (A.rows() > 3) ? A[3][3] : 1.0;
    
  // Sample points from the medial surface and apply rotation
  for(size_t i = 0; i < ncu; i++) for(size_t j = 0; j < ncv; j++)
    {
    SMLVec3d y, x( C(0, i, j), C(1, i, j), C(2, i, j));

    if(i == 0 && j == 0)
      {
      y = A3 * (x - C3) + B3 + C3;
      }
    else
      {
      y = A3 * x;
      }

    // Update the coefficients
    C(0, i, j) = y[0]; C(1, i, j) = y[1]; C(2, i, j) = y[2]; 
    
    // Scale the rho coefficient
    C(3, i, j) *= xRhoScale;
    }
}

/** Get a coefficient mask of a given coarseness in each component.
 * Coarseness of zero should mean minimal number of coefficients, and
 * coarseness of one is the maximal numner of coefficients */
IMedialCoefficientMask * FourierSurface::NewCoarseToFineCoefficientMask(double *xCoarseness)
{
  // Get the number of coefficients in each direction that has been requested
  vector<size_t> iSelect;
  for(size_t iComp = 0; iComp < 4; iComp ++)
    {
    size_t nu = (size_t) (xCoarseness[iComp] * ncu);
    size_t nv = (size_t) (xCoarseness[iComp] * ncv);
    for(size_t i = 0; i < nu; i++) for(size_t j = 0; j < nv; j++)
      iSelect.push_back(GetCoefficientIndex(i,j,iComp));
    
    cout << "Component " << iComp << ", Taking " << nu << " x " << nv << " coeffs. " << endl;
    }



  // Create the mask
  return new SelectionMedialCoefficientMask(this, iSelect);
}


IHyperSurface2D *
FourierSurface
::GetVariationSurface(const double *xCoeff)
{
  // Create a new surface with the same parameters as this one
  FourierSurface *xVariation = new FourierSurface(ncu, ncv);

  // Set the components of the variation to the components passed in; 
  // we can do this because the function computed here is linear in the 
  // coefficients
  xVariation->SetCoefficientArray(xCoeff);

  // Copy the evaluation grid
  VectorType uu, vv;
  GetEvaluationGrid(uu, vv);
  xVariation->SetEvaluationGrid(uu, vv);

  // Return this variation  
  return xVariation;
}

void
FourierSurface
::ReleaseVariationSurface(IHyperSurface2D *xSurface)
{
  delete xSurface;
}







template class GenericBasisRepresentation2D< 4, 3,        
         CosineBasisFunction, CosineBasisFunction>;
