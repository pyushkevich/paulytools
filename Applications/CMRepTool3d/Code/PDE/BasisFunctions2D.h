#ifndef __BasisFunctions2D_h_
#define __BasisFunctions2D_h_

#include <valarray>
#include <smlmath.h>
#include "Registry.h"

using namespace std;

class IHyperSurface2D
{
public:
  /** Evaluate the function at a particular u, v coordinate. The return value
   * is a vector of doubles */
  virtual void Evaluate(double u, double v, double *x) = 0;
 
  /** Evaluate a partial derivative of the function at the u, v coordinate.
   * The first component is c0, number of components nc */
  virtual void EvaluateDerivative(
    double u, double v, size_t ou, size_t ov, size_t c0, size_t nc, double *x) = 0;

  /** Specify a grid of values along which the function will need to be
   * evaluated repeatedly. */
  virtual void SetEvaluationGrid(size_t nu, size_t nv, double *uu, double *vv) = 0;

  /** Evaluate the function at a grid point */
  virtual void EvaluateAtGridIndex(
    size_t iu, size_t iv, size_t ou, size_t ov, size_t c0, size_t nc, double *x) = 0;
};

/** This is an interface that should be implemented by any 2D basis function
 * implementation */
class IBasisRepresentation2D : virtual public IHyperSurface2D
{
public:

  /** Get the coefficients as an array of real values */
  virtual size_t GetNumberOfRawCoefficients() = 0;
  virtual double *GetRawCoefficientArray() = 0;
  virtual double GetRawCoefficient(size_t iCoeff) = 0;

  /** Set the coefficient array */
  virtual void SetRawCoefficientArray(const double *array) = 0;
  virtual void SetRawCoefficient(size_t iCoeff, double xValue) = 0;
  
  /** Fit the i-th component of the coefficients to some data points */
  virtual void FitToData(size_t n, size_t i, double *uu, double *vv, double *xx) = 0;

  /** Save the coefficients to a registry */
  virtual void SaveToRegistry(Registry &R) = 0;

  /** Read the coefficients from a registry */
  virtual bool ReadFromRegistry(Registry &R) = 0;

  /** Apply affine transform to the surface */
  virtual void ApplyAffineTransform(
    const vnl_matrix<double> &A, const vnl_vector<double> &b, 
    const vnl_vector<double> &c) = 0;
};
  
/**
 * A 3D index into an array. As the array is traversed, the x increases
 * fastest, z slowest
 */
class Index3D : public valarray<double> 
{
public:
  Index3D(size_t nx, size_t ny, size_t nz) 
    : valarray<double> (nx * ny * nz)
    { stride_y = nx; stride_z = nx * ny; }

  double &operator() (size_t ix, size_t iy, size_t iz)
    { return (*this)[iz * stride_z + iy * stride_y + ix]; }

  double *GetPointer()
    { return &((*this)[0]); }

private:
  size_t stride_z, stride_y;
};

/** 
 * A generic implementation that uses some arbitrary basis function class.
 * This representation is modeled after Fourier basis functions, but it could
 * be arbitrary. The basic idea is that the 2D basis functions are separable
 * into 1D components, i.e., F_k(u,v) = f_i(u) f_j(v)
 */
template< size_t NComponents, size_t NOrder, 
  typename BasisFunctionU, typename BasisFunctionV >
class GenericBasisRepresentation2D : virtual public IBasisRepresentation2D
{
public:

  /** Construct the representation with a specific number of coefficients in u
   * and v components */
  GenericBasisRepresentation2D(size_t ncu, size_t ncv);
  
  /** Evaluate the function at a particular u, v coordinate. The return value
   * is a vector of doubles */
  void Evaluate(double u, double v, double *x)
    { return EvaluateDerivative(u, v, 0, 0, 0, NComponents, x); } 
 
  /** Evaluate a partial derivative of the function at the u, v coordinate */
  void EvaluateDerivative(double u, double v, size_t ou, size_t ov, size_t c0, size_t nc, double *x); 

  /** Specify a grid of values along which the function will need to be
   * evaluated repeatedly. Call PrecomputeGrid() before evaluating at grid and
   * after changing any of the coefficients */
  void SetEvaluationGrid(size_t nu, size_t nv, double *uu, double *vv); 

  /** Evaluate the function at a grid point */
  void EvaluateAtGridIndex(size_t iu, size_t iv, size_t ou, size_t ov, size_t c0, size_t nc, double *x); 

  /** Get number of raw coefficients */
  size_t GetNumberOfRawCoefficients()
    { return C.size(); }

  /** Get the coefficients as an array of real values */
  double *GetRawCoefficientArray()
    { return C.GetPointer(); }

  /** Set the array of coefficients */
  void SetRawCoefficientArray(const double *data)
    { for(size_t i = 0; i < C.size(); i++) C[i] = data[i]; }

  // Direct coefficient access
  double GetRawCoefficient(size_t iCoeff) 
    { return C[iCoeff]; }
  
  void SetRawCoefficient(size_t iCoeff, double xValue) 
    { C[iCoeff] = xValue; }

  /** Fit the i-th component of the coefficients to some data points */
  void FitToData(size_t n, size_t i, double *uu, double *vv, double *xx); 

  /** Save the coefficients to a registry */
  void SaveToRegistry(Registry &R);

  /** Read the coefficients from a registry */
  bool ReadFromRegistry(Registry &R);

protected:
  // The raw array of coefficients (without pointers)
  Index3D C;

  // The evaluation grids along u and v. 
  valarray<double> uGrid, vGrid, uGridValues, vGridValues, uEval, vEval;

  // Number of coefficients, raw coefficients
  size_t ncu, ncv, nc, ncRaw;

  // The basis functions
  BasisFunctionU fu; 
  BasisFunctionV fv;

  // Get the slice corresponding to coefficient iu, iv
  slice GetCoefficientIndex(size_t iu, size_t iv)
    { return slice(((ncu * iv) + iu) * NComponents, NComponents, 1); }

  size_t GetCoefficientIndex(size_t iu, size_t iv, size_t k)
    { return ((ncu * iv) + iu) * NComponents + k; }

  void Initialize(size_t ncu, size_t ncv);
};

/** The cosine basis function */
class CosineBasisFunction
{
public:
  double Evaluate(double u, size_t k, size_t d);
};

/** Declare the Fourier surface class */
typedef 
GenericBasisRepresentation2D <4, 3, CosineBasisFunction, CosineBasisFunction>
FourierSurfaceBase;

class FourierSurface : public FourierSurfaceBase
{
public:
  FourierSurface(size_t ncu, size_t ncv) : FourierSurfaceBase(ncu, ncv) {}

  // Special way to apply the affine transform
  void ApplyAffineTransform( const vnl_matrix<double> &A, 
    const vnl_vector<double> &b, const vnl_vector<double> &c);
};

#endif

