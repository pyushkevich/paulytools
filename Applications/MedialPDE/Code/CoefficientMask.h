#ifndef __CoefficientMask_h_
#define __CoefficientMask_h_

#include <vector>
#include <smlmath.h>
#include "BasisFunctions2D.h"

using namespace std;

class IMedialCoefficientMask : public ICoefficientSettable
{
public:
};

class SelectionMedialCoefficientMask : public IMedialCoefficientMask
{
public:
  SelectionMedialCoefficientMask(
    ICoefficientSettable *source, const vector<size_t> &mask);

  SelectionMedialCoefficientMask(
    ICoefficientSettable *source, size_t iMaskSize, const size_t *iMaskIndex);

  double *GetCoefficientArray();
  void SetCoefficientArray(const double *xData);
  size_t GetNumberOfCoefficients() const
    { return nCoeff; }
  double GetCoefficient(size_t i) const; 
  void SetCoefficient(size_t i, double x);

  // Get a surface corresponding to a single component
  IHyperSurface2D *GetComponentSurface(size_t iCoefficient);
  void ReleaseComponentSurface(IHyperSurface2D *xSurface);

  // Get a surface corresponding to some variation
  IHyperSurface2D *GetVariationSurface(const double *xData);
  void ReleaseVariationSurface(IHyperSurface2D *xSurface);

private:
  ICoefficientSettable *xSource;
  vnl_vector<double> xRawCoefficients, xMaskCoefficients;
  vector<size_t> xMask;
  size_t nCoeff;
};

class PassThroughCoefficientMask : public virtual IMedialCoefficientMask
{
public:
  PassThroughCoefficientMask(ICoefficientSettable *source);
  double *GetCoefficientArray();
  void SetCoefficientArray(const double *xData);
  size_t GetNumberOfCoefficients() const;
  double GetCoefficient(size_t i) const;
  void SetCoefficient(size_t i, double x);

  // Get a surface corresponding to a single component
  IHyperSurface2D *GetComponentSurface(size_t iCoefficient);
  void ReleaseComponentSurface(IHyperSurface2D *xSurface);

  // Get a surface corresponding to some variation
  IHyperSurface2D *GetVariationSurface(const double *xData);
  void ReleaseVariationSurface(IHyperSurface2D *xSurface);
  
private:
  ICoefficientSettable *xSource;
};

class AffineTransformCoefficientMask : virtual public IMedialCoefficientMask
{
public:
  AffineTransformCoefficientMask(IBasisRepresentation2D *surface);

  double *GetCoefficientArray()
    { return xData.data_block(); }
  
  void SetCoefficientArray(const double *inData);
  
  size_t GetNumberOfCoefficients() const
    { return nCoeff; }

  double GetCoefficient(size_t i) const
    { return xData[i]; }

  void SetCoefficient(size_t i, double x)
    { xData[i] = x; SetCoefficientArray(xData.data_block()); }

  // Get a surface corresponding to a single component
  IHyperSurface2D *GetComponentSurface(size_t iCoefficient);
  void ReleaseComponentSurface(IHyperSurface2D *xSurface);

  // Get a surface corresponding to some variation
  IHyperSurface2D *GetVariationSurface(const double *xData);
  void ReleaseVariationSurface(IHyperSurface2D *xSurface);

private:
  typedef vnl_matrix<double> MatrixType;
  typedef vnl_vector<double> VectorType;

  
  IBasisRepresentation2D *xSurface;
  size_t nDim, iIndexB, nCoeff;
  
  VectorType xOriginalCoefficients;
  MatrixType A;
  VectorType b, c;
  VectorType xData;
  double xRhoScale;
};


/** 
 * This is an affine transform in 3x1D, it uses the standard 12 rather than
 * 20 coefficients. It's basically a selection mask wrapped around an affine
 * mask
 */
class AffineTransform3DCoefficientMask :
  public SelectionMedialCoefficientMask
{
public:
  AffineTransform3DCoefficientMask(IBasisRepresentation2D *surface) : 
    SelectionMedialCoefficientMask(
      xAffineMask = new AffineTransformCoefficientMask(surface), 12, xIndexArray) 
  { }

  ~AffineTransform3DCoefficientMask()
    { delete xAffineMask; }
    
private:
  AffineTransformCoefficientMask *xAffineMask;
  static const size_t xIndexArray[];
};



#endif //__CoefficientMask_h_

