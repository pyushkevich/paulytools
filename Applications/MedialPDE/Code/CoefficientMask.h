#ifndef __CoefficientMask_h_
#define __CoefficientMask_h_

#include <vector>
#include <smlmath.h>

class IBasisRepresentation2D;

using namespace std;

class IMedialCoefficientMask
{
public:
  virtual double *GetCoefficientArray() = 0;
  virtual void SetCoefficientArray(double *xData) = 0;
  virtual size_t GetNumberOfCoefficients() = 0;
  virtual double GetCoefficient(size_t i) = 0;
  virtual void SetCoefficient(size_t i, double x) = 0;
};

class SelectionMedialCoefficientMask : public IMedialCoefficientMask
{
public:
  SelectionMedialCoefficientMask(
    IBasisRepresentation2D *source, const vector<size_t> &mask);

  double *GetCoefficientArray();
  void SetCoefficientArray(double *xData);
  size_t GetNumberOfCoefficients()
    { return nCoeff; }
  double GetCoefficient(size_t i);
  void SetCoefficient(size_t i, double x);

private:
  IBasisRepresentation2D *xSource;
  vnl_vector<double> xRawCoefficients, xMaskCoefficients;
  vector<size_t> xMask;
  size_t nCoeff;
};

class PassThroughCoefficientMask : public virtual IMedialCoefficientMask
{
public:
  PassThroughCoefficientMask(IBasisRepresentation2D *source);
  double *GetCoefficientArray();
  void SetCoefficientArray(double *xData);
  size_t GetNumberOfCoefficients();
  double GetCoefficient(size_t i);
  void SetCoefficient(size_t i, double x);
private:
  IBasisRepresentation2D *xSource;
};

class AffineTransformCoefficientMask : virtual public IMedialCoefficientMask
{
public:
  AffineTransformCoefficientMask(IBasisRepresentation2D *surface);

  double *GetCoefficientArray()
    { return xData.data_block(); }
  
  void SetCoefficientArray(double *inData);
  
  size_t GetNumberOfCoefficients()
    { return nCoeff; }

  double GetCoefficient(size_t i)
    { return xData[i]; }

  void SetCoefficient(size_t i, double x)
    { xData[i] = x; SetCoefficientArray(xData.data_block()); }

private:
  IBasisRepresentation2D *xSurface;
  size_t nDim, iIndexB, nCoeff;
  
  vnl_vector<double> xOriginalCoefficients;
  vnl_matrix<double> A;
  vnl_vector<double> b, c;
  vnl_vector<double> xData;
  double xRhoScale;
};
#endif //__CoefficientMask_h_
