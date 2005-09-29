#ifndef __BasisFunctions2D_h_
#define __BasisFunctions2D_h_

#include <valarray>
#include <smlmath.h>
#include "Registry.h"

using namespace std;

class IMedialCoefficientMask;

class IHyperSurface2D
{
public:
  // Typed definitions
  typedef vnl_matrix<double> MatrixType;
  typedef vnl_vector<double> VectorType;

  /** Evaluate the function at a particular u, v coordinate. The return value
   * is a vector of doubles */
  virtual void Evaluate(double u, double v, double *x) = 0;
 
  /** Evaluate a partial derivative of the function at the u, v coordinate.
   * The first component is c0, number of components nc */
  virtual void EvaluateDerivative(
    double u, double v, size_t ou, size_t ov, size_t c0, size_t nc, double *x) = 0;

  /** Evaluate the function at a grid point */
  virtual void EvaluateAtGridIndex(
    size_t iu, size_t iv, size_t ou, size_t ov, size_t c0, size_t nc, double *x) = 0;

  /** Get the evaluation grid parameters */
  virtual void GetEvaluationGrid(VectorType &uu, VectorType &vv) const = 0;

  /** Get the number of dimensions of the surface */
  virtual size_t GetNumberOfDimensions() = 0;

  /** Print a report about this surface */
  virtual void PrintReport() {};
};

class IMutableHyperSurface2D : virtual public IHyperSurface2D
{
public:
  /** Specify a grid of values along which the function will need to be
   * evaluated repeatedly. */
  virtual void SetEvaluationGrid(const VectorType &uu, const VectorType &vv) = 0;

  /** Apply affine transform to the hypersurface */
  virtual void ApplyAffineTransform(
    const MatrixType &A, const VectorType &b, const VectorType &c) = 0;

  /** Find the center of the surface (for rotations) */
  virtual VectorType GetCenterOfRotation() = 0;
};

/**
 * This is an interface that contains functions for setting and getting
 * coefficients and corresponding component surfaces
 */
class ICoefficientSettable
{
public:
  virtual double *GetCoefficientArray() = 0;
  virtual void SetCoefficientArray(const double *xData) = 0;
  virtual size_t GetNumberOfCoefficients() const = 0 ;
  virtual double GetCoefficient(size_t i) const = 0;
  virtual void SetCoefficient(size_t i, double x) = 0;

  // Get a surface corresponding to a single component
  virtual IHyperSurface2D *GetComponentSurface(size_t iCoefficient) = 0;
  virtual void ReleaseComponentSurface(IHyperSurface2D *xSurface) = 0;
};

/** This is an interface that should be implemented by any 2D basis function
 * implementation */
class IBasisRepresentation2D : 
  virtual public IMutableHyperSurface2D, 
  virtual public ICoefficientSettable
{
public:

  /** Fit the i-th component of the coefficients to some data points */
  virtual void FitToData(size_t n, size_t i, double *uu, double *vv, double *xx) = 0;

  /** Save the coefficients to a registry */
  virtual void SaveToRegistry(Registry &R) = 0;

  /** Read the coefficients from a registry */
  virtual bool ReadFromRegistry(Registry &R) = 0;
  
  /** Get a coefficient mask of a given coarseness in each component.
   * Coarseness of zero should mean minimal number of coefficients, and
   * coarseness of one is the maximal numner of coefficients */
  virtual IMedialCoefficientMask *
    NewCoarseToFineCoefficientMask(double *xCoarseness) = 0;
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

  void resize(size_t nx, size_t ny, size_t nz)
    {
    stride_y = nx; 
    stride_z = nx * ny; 
    valarray<double>::resize(nx * ny * nz);
    (*this) *= 0.0; 
    }
    

  double &operator() (size_t ix, size_t iy, size_t iz) 
    { return (*this)[iz * stride_z + iy * stride_y + ix]; }

  const double &operator() (size_t ix, size_t iy, size_t iz) const
    { return (*this)[iz * stride_z + iy * stride_y + ix]; }

  double *GetPointer()
    { return &((*this)[0]); }

  const double *GetPointer() const
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

  /** Class that allows you to return a single basis function as a surface;
   * used in variational calculus computations */
  class SingleFunctionAdapter : public IHyperSurface2D
  {
  public:
    void Evaluate(double u, double v, double *x);
    void EvaluateDerivative(
      double u, double v, size_t ou, size_t ov, size_t c0, size_t nc, double *x);

    void GetEvaluationGrid(VectorType &uu, VectorType &vv) const;
    void EvaluateAtGridIndex(
      size_t iu, size_t iv, size_t ou, size_t ov, size_t c0, size_t nc, double *x);
    size_t GetNumberOfDimensions();

  private:
    GenericBasisRepresentation2D *xParent;
    size_t icu, icv, iComp;

    friend class GenericBasisRepresentation2D;
  };
      
  /** Construct the representation with a specific number of coefficients in u
   * and v components */
  GenericBasisRepresentation2D(size_t ncu, size_t ncv);

  void SetNumberOfCoefficients(size_t ncu, size_t ncv);
  
  /** Evaluate the function at a particular u, v coordinate. The return value
   * is a vector of doubles */
  void Evaluate(double u, double v, double *x)
    { return EvaluateDerivative(u, v, 0, 0, 0, NComponents, x); } 
 
  /** Evaluate a partial derivative of the function at the u, v coordinate */
  void EvaluateDerivative(double u, double v, size_t ou, size_t ov, size_t c0, size_t nc, double *x); 

  /** Specify a grid of values along which the function will need to be
   * evaluated repeatedly. Call PrecomputeGrid() before evaluating at grid and
   * after changing any of the coefficients */
  void SetEvaluationGrid(const VectorType &uu, const VectorType &vv);

  /** Get the vectors that define the evaluation grid */
  void GetEvaluationGrid(VectorType &uu, VectorType &vv) const;

  /** Evaluate the function at a grid point */
  void EvaluateAtGridIndex(
    size_t iu, size_t iv, size_t ou, size_t ov, size_t c0, size_t nc, double *x); 

  /** Get the number of components/dimensions */
  size_t GetNumberOfDimensions() 
    { return NComponents; }

  /** Get number of raw coefficients */
  size_t GetNumberOfCoefficients() const 
    { return C.size(); }

  /** Get the coefficients as an array of real values */
  double *GetCoefficientArray() 
    { return C.GetPointer(); }
  
  /** Get the coefficients as an array of real values */
  const double *GetCoefficientArray() const
    { return C.GetPointer(); }

  /** Set the array of coefficients */
  void SetCoefficientArray(const double *data)
    { for(size_t i = 0; i < C.size(); i++) C[i] = data[i]; }

  // Direct coefficient access
  double GetCoefficient(size_t iCoeff) const
    { return C[iCoeff]; }
  
  void SetCoefficient(size_t iCoeff, double xValue) 
    { C[iCoeff] = xValue; }

  void SetAllCoefficients(double xValue)
    {  
    for(size_t i=0; i < C.size(); i++)
      C[i] = xValue;
    }

  // Set a coefficient by index
  void SetCoefficient(size_t iu, size_t iv, size_t iComp, double value)
    { C(iComp,iu,iv) = value; }

  double GetCoefficient(size_t iu, size_t iv, size_t iComp) const
    { return C(iComp,iu,iv); }

  /** Fit the i-th component of the coefficients to some data points */
  void FitToData(size_t n, size_t i, double *uu, double *vv, double *xx); 

  /** Save the coefficients to a registry */
  void SaveToRegistry(Registry &R);

  /** Read the coefficients from a registry */
  bool ReadFromRegistry(Registry &R);

  /** Get the number of coefficients in U and V */
  void GetNumberOfCoefficientsUV(size_t &ncu, size_t &ncv)
    { ncu = this->ncu; ncv = this->ncv; }

  /** Debug report */
  void PrintReport();

  /** Get one of the components as a surface */
  IHyperSurface2D *GetComponentSurface(size_t iCoefficient);

  /** Free the memory associated with the component surface */
  void ReleaseComponentSurface(IHyperSurface2D *xSurface);

  /** Get the component by u, v, comp index */
  IHyperSurface2D *GetComponentSurface(size_t icu, size_t icv, size_t iComp);

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

  void GetRawCoefficientIndices(size_t iRaw, size_t &icu, size_t &icv, size_t &iComp)
    {
    iComp = iRaw % NComponents;
    icu = (iRaw / NComponents) % ncu;
    icv = iRaw / (NComponents * ncu);
    }

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
  // Constructor
  FourierSurface(size_t ncu, size_t ncv) : FourierSurfaceBase(ncu, ncv) {}

  // Copy constructor
  FourierSurface(const FourierSurface &source); 

  /** Find the center of the surface (for rotations) */
  vnl_vector<double> GetCenterOfRotation();

  // Special way to apply the affine transform
  void ApplyAffineTransform( const vnl_matrix<double> &A, 
    const vnl_vector<double> &b, const vnl_vector<double> &c);
  
  /** Get a coefficient mask of a given coarseness in each component.
   * Coarseness of zero should mean minimal number of coefficients, and
   * coarseness of one is the maximal numner of coefficients */
  IMedialCoefficientMask *NewCoarseToFineCoefficientMask(double *xCoarseness);
};

#endif

