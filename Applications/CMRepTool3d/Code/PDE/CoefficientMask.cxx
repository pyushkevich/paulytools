#include "CoefficientMask.h"
#include "BasisFunctions2D.h"

SelectionMedialCoefficientMask::SelectionMedialCoefficientMask(
  IBasisRepresentation2D *source, const vector<size_t> &mask)
: xMask(mask), 
  xRawCoefficients(mask.size()), xMaskCoefficients(mask.size())
{
  xSource = source;
  nCoeff = mask.size();
}

/** Get the coefficient array */
double *SelectionMedialCoefficientMask::GetCoefficientArray()
{
  double *xRaw = xSource->GetRawCoefficientArray();
  for(size_t i=0; i < nCoeff; i++)
    { xMaskCoefficients[i] = xRaw[ xMask[i] ]; }
  return &(xMaskCoefficients[0]);
}

/** Set the coefficient array */
void SelectionMedialCoefficientMask::SetCoefficientArray(double *xData)
{
  double *xRaw = xSource->GetRawCoefficientArray();
  for(size_t i=0; i < nCoeff; i++)
    {
    xRaw[ xMask[i] ] = xData[i];
    }
}

double SelectionMedialCoefficientMask::GetCoefficient(size_t i)
{
  return xSource->GetRawCoefficientArray()[xMask[i]];
}

void SelectionMedialCoefficientMask::SetCoefficient(size_t i, double x)
{ 
  xSource->SetRawCoefficient( xMask[i], x ); 
}
  

PassThroughCoefficientMask::PassThroughCoefficientMask(
  IBasisRepresentation2D *source) 
: xSource(source) 
{
}

double *PassThroughCoefficientMask::GetCoefficientArray()
{ 
  return xSource->GetRawCoefficientArray(); 
}

void PassThroughCoefficientMask::SetCoefficientArray(double *xData)
{ 
  xSource->SetRawCoefficientArray(xData); 
}

size_t PassThroughCoefficientMask::GetNumberOfCoefficients()
{ 
  return xSource->GetNumberOfRawCoefficients(); 
}

double PassThroughCoefficientMask::GetCoefficient(size_t i)
{ 
  return xSource->GetRawCoefficient(i); 
}

void PassThroughCoefficientMask::SetCoefficient(size_t i, double x)
{ 
  xSource->SetRawCoefficient(i, x); 
}

AffineTransformCoefficientMask
::AffineTransformCoefficientMask(IBasisRepresentation2D *surface)
: xOriginalCoefficients(
    surface->GetRawCoefficientArray(), surface->GetNumberOfRawCoefficients()), 
  xSurface(surface), 
  nDim(surface->GetNumberOfDimensions()),
  A(nDim, nDim), b(nDim), xData(nDim * nDim + nDim)
{ 
  // Initialize the affine transform
  A.set_identity(); b.fill(0);
  c = surface->GetCenterOfRotation();

  // Copy it to the data vector
  A.copy_out(xData.data_block());

  iIndexB = nDim * nDim; 
  nCoeff = iIndexB + nDim;
  b.copy_out(xData.data_block() + iIndexB);
}

void AffineTransformCoefficientMask::SetCoefficientArray(double *inData)
{ 
  xData.copy_in(inData);
  A.copy_in(inData);
  b.copy_in(inData + iIndexB);

  xSurface->SetRawCoefficientArray(xOriginalCoefficients.data_block());
  xSurface->ApplyAffineTransform(A, b, c);
}

