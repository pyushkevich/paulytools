#include "CoefficientMask.h"
#include "BasisFunctions2D.h"

/******************************************************************************
 * SelectionMedialCoefficientMask 
 *****************************************************************************/

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

IHyperSurface2D *
SelectionMedialCoefficientMask
::GetComponentSurface(size_t iCoefficient)
{
  return xSource->GetComponentSurface( xMask[iCoefficient] );
}

void
SelectionMedialCoefficientMask
::ReleaseComponentSurface(IHyperSurface2D *xSurface)
{
  xSource->ReleaseComponentSurface(xSurface);
}

/******************************************************************************
 * PassThroughCoefficientMask 
 *****************************************************************************/

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

IHyperSurface2D *
PassThroughCoefficientMask
::GetComponentSurface(size_t iCoefficient)
{
  return xSource->GetComponentSurface(iCoefficient);
}

void
PassThroughCoefficientMask
::ReleaseComponentSurface(IHyperSurface2D *xSurface)
{
  xSource->ReleaseComponentSurface(xSurface);
}

/******************************************************************************
 * AffineComponentSurface 
 *****************************************************************************/

class AffineComponentSurface : public IHyperSurface2D
{
public:
  // Typed definitions
  typedef vnl_matrix<double> MatrixType;
  typedef vnl_vector<double> VectorType;

  // Constructor: pass in matrix derivatives
  AffineComponentSurface(
    IHyperSurface2D *xSource, size_t nDim,
    const MatrixType &A, const VectorType &b, const VectorType &c)
    {
    this->A = A;
    this->b = b;
    this->c = c;
    this->xSource = xSource;
    this->nDim = nDim;
    }

  // Apply the affine transform (this is in the same way as in FourierSurface)
  void ApplyAffineTransform(VectorType &X, VectorType &Y, bool flagDeriv)
    {
    if(flagDeriv)
      { X[0] -= c[0]; X[1] -= c[1]; X[2] -= c[2]; }
    
    Y[0] = A[0][0] * X[0] + A[0][1] * X[1] + A[0][2] * X[2];
    Y[1] = A[1][0] * X[0] + A[1][1] * X[1] + A[1][2] * X[2];
    Y[2] = A[2][0] * X[0] + A[2][1] * X[1] + A[2][2] * X[2];
    Y[3] = X[3] * A[3][3];

    if(flagDeriv)
      { Y[0] += b[0]; Y[1] += b[1]; Y[2] += b[2]; }
    }

  /** Evaluate the function at a particular u, v coordinate. The return value
   * is a vector of doubles */
  void Evaluate(double u, double v, double *x)
    {
    VectorType X(nDim), Y(nDim);
    xSource->Evaluate(u, v, X.data_block());
    ApplyAffineTransform(X, Y, true);
    Y.copy_out(x);
    }
 
  /** Evaluate a partial derivative of the function at the u, v coordinate.
   * The first component is c0, number of components nc */
  void EvaluateDerivative(
    double u, double v, size_t ou, size_t ov, size_t c0, size_t nc, double *x)
    {
    VectorType X(nDim), Y(nDim);

    // Compute the surface point
    xSource->EvaluateDerivative(u, v, ou, ov, 0, nDim, X.data_block());
    
    // Apply the affine transform
    ApplyAffineTransform(X, Y, ou + ov == 0);

    // Copy the values into x
    for(size_t i = 0; i < nc; i++)
      x[i] =  Y[i + c0];
    }

  /** Evaluate the function at a grid point */
  void EvaluateAtGridIndex(
    size_t iu, size_t iv, size_t ou, size_t ov, size_t c0, size_t nc, double *x)
    {
    VectorType X(nDim), Y(nDim);

    // Compute the surface point
    xSource->EvaluateAtGridIndex(iu, iv, ou, ov, 0, nDim, X.data_block());
    
    // Apply the affine transform
    ApplyAffineTransform(X, Y, ou + ov == 0);

    // Copy the values into x
    for(size_t i = 0; i < nc; i++)
      x[i] =  Y[i + c0];
    }

  /** Get the evaluation grid parameters */
  void GetEvaluationGrid(VectorType &uu, VectorType &vv) const
    { return xSource->GetEvaluationGrid(uu, vv); }

  /** Get the number of dimensions of the surface */
  size_t GetNumberOfDimensions()
    { return xSource->GetNumberOfDimensions(); }

private:
  size_t nDim;
  IHyperSurface2D *xSource;
  MatrixType A;
  VectorType b, c;
};

/******************************************************************************
 * AffineTransformCoefficientMask 
 *****************************************************************************/

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

IHyperSurface2D *
AffineTransformCoefficientMask
::GetComponentSurface(size_t iCoefficient)
{
  // Create the derivative coefficient vector
  VectorType xIndex(nCoeff, 0.0);
  xIndex[iCoefficient] = 1.0;
  
  // Compute the matrices that are to be applied
  MatrixType dA(nDim, nDim);
  VectorType db(nDim);

  dA.copy_in(xIndex.data_block());
  db.copy_in(xIndex.data_block() + iIndexB);

  // Create the special surface
  return new AffineComponentSurface(xSurface, nDim, dA, db, c);
}

void
AffineTransformCoefficientMask
::ReleaseComponentSurface(IHyperSurface2D *xSurface)
{
  delete xSurface;
}

