#ifndef _itkWindowedSincInterpolateImageFunction_h
#define _itkWindowedSincInterpolateImageFunction_h

#include "itkConstNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkInterpolateImageFunction.h"

namespace itk
{

namespace Function {

/** 
 * \class CosineWindowFunction
 * \brief Window function for sinc interpolation.
 * \f[ w(x) = cos ( \frac{\pi x}{2 m} ) \f]
 * \sa WindowedSincInterpolateImageFunction 
 */
template< unsigned int VRadius, 
  class TInput=double, class TOutput=double>
class CosineWindowFunction
{
public:
  inline TOutput operator()( const TInput & A ) const
    { return (TOutput) cos( A * m_Factor ); }
private:
  /** Equal to \f$ \frac{\pi}{2 m} \f$ */
  static const double m_Factor;
}; 

/** 
 * \class HammingWindowFunction
 * \brief Window function for sinc interpolation.
 * \f[ w(x) = 0.54 + 0.46 cos ( \frac{\pi x}{m} ) \f]
 * \sa WindowedSincInterpolateImageFunction 
 */
template< unsigned int VRadius, 
  class TInput=double, class TOutput=double>
class HammingWindowFunction
{
public:
  inline TOutput operator()( const TInput & A ) const
   { return (TOutput) 0.54 + 0.46 * cos( A * m_Factor ); }
private:
  /** Equal to \f$ \frac{\pi}{m} \f$ */
  static const double m_Factor;
}; 

/** 
 * \class WelchWindowFunction
 * \brief Window function for sinc interpolation.
 * \f[ w(x) = 1 - ( \frac{x^2}{m^2} ) \f]
 * \sa WindowedSincInterpolateImageFunction 
 */
template< unsigned int VRadius, 
  class TInput=double, class TOutput=double>
class WelchWindowFunction
{
public:
  inline TOutput operator()( const TInput & A ) const
   { return (TOutput) (1.0 - A * m_Factor * A); }
private:
  /** Equal to \f$ \frac{1}{m^2} \f$ */
  static const double m_Factor;
}; 

/** 
 * \class LancosWindowFunction
 * \brief Window function for sinc interpolation.
 * \f[ w(x) = \textrm{sinc} ( \frac{x}{m} ) \f]
 * Note: Paper referenced in WindowedSincInterpolateImageFunction gives
 * an incorrect definition of this window function. 
 * \sa WindowedSincInterpolateImageFunction 
 */
template< unsigned int VRadius, 
  class TInput=double, class TOutput=double>
class LancosWindowFunction
{
public:
  inline TOutput operator()( const TInput & A ) const
    {
    if(A == 0.0) return (TOutput) 1.0;
    double z = m_Factor * A;
    return (TOutput) ( sin(z) / z ); 
    }
private:
  /** Equal to \f$ \frac{\pi}{m} \f$ */
  static const double m_Factor;
}; 

/** 
 * \class BlackmanWindowFunction
 * \brief Window function for sinc interpolation.
 * \f[ w(x) = 0.42 + 0.5 cos(\frac{\pi x}{m}) + 0.08 cos(\frac{2 \pi x}{m}) \f]
 * \sa WindowedSincInterpolateImageFunction 
 */
template< unsigned int VRadius, 
  class TInput=double, class TOutput=double>
class BlackmanWindowFunction
{
public:
  inline TOutput operator()( const TInput & A ) const
    {
    return (TOutput)
      (0.42 + 0.5 * cos(A * m_Factor1) + 0.08 * cos(A * m_Factor2));
    }
private:
  /** Equal to \f$ \frac{\pi}{m} \f$ */
  static const double m_Factor1;
  
  /** Equal to \f$ \frac{2 \pi}{m} \f$  */
  static const double m_Factor2;
}; 

} // namespace Function

/** 
  \class WindowedSincInterpolateImageFunction
  \brief Use the windowed sinc function to interpolate
  \author Paul A. Yushkevich
  
  \par THEORY

  This function is intended to provide an interpolation function that 
  has minimum aliasing artifacts, in contrast to linear interpolation.
  According to sampling theory, the infinite-support sinc filter, 
  whose Fourier transform is the box filter, is optimal for resampling 
  a function. In practice, the infinite support sinc filter is 
  approximated using a limited support 'windowed' sinc filter. 

  \par
  This function is based on the following publication:

  \par
  Erik H. W. Meijering, Wiro J. Niessen, Josien P. W. Pluim,
  Max A. Viergever: Quantitative Comparison of Sinc-Approximating
  Kernels for Medical Image Interpolation. MICCAI 1999, pp. 210-217

  \par
  In this work, several 'windows' are estimated. In two dimensions, the 
  interpolation at a position (x,y) is given by the following
  expression: 

  \par
  \f[
    I(x,y) = 
      \sum_{i = \lfloor x \rfloor + 1 - m}^{\lfloor x \rfloor + m} 
      \sum_{j = \lfloor y \rfloor + 1 - m}^{\lfloor y \rfloor + m}
      I_{i,j} K(x-i) K(y-j),
  \f]

  \par
  where m is the 'radius' of the window, (3,4 are reasonable numbers),
  and K(t) is the kernel function, composed of the sinc function and
  one of several possible window functions:

  \par
  \f[
    K(t) = w(t) \textrm{sinc}(t) = w(t) \frac{\sin(\pi t)}{\pi t}
  \f]

  \par
  Several window functions are provided here in the itk::Function
  namespace. The conclusions of the referenced paper suggest to use the
  Welch, Cosine, Kaiser, and Lancos windows for m = 4,5. These are based
  on error in rotating medical images w.r.t. the linear interpolation
  method. In some cases the results achieve a 20-fold improvement in
  accuracy.

  \par USING THIS FILTER

  Use this filter the way you would use any ImageInterpolationFunction,
  so for instance, you can plug it into the ResampleImageFilter class.
  In order to initialize the filter you must choose several template
  parameters. 

  \par
  The first (TInputImage) is the image type, that's standard. 

  \par
  The second (VRadius) is the radius of the kernel, i.e., the 
  \f$ m \f$ from the formula above. 

  \par
  The third (TWindowFunction) is the window function object, which you 
  can choose from about five different functions defined in this 
  header. The default is the Hamming window, which is commonly used
  but not optimal according to the cited paper. 

  \par
  The fourth (TBoundaryCondition) is the boundary condition class used
  to determine the values of pixels that fall off the image boundary.
  This class has the same meaning here as in the NeighborhoodItetator
  classes. 

  \par
  The fifth (TCoordRep) is again standard for interpolating functions,
  and should be float or double.

  \par CAVEATS

  There are a few improvements that an enthusiasting ITK developer could
  make to this filter. First, there is a performance issue associated
  with creating an internal copy of the m_Neighborhood object in the
  method EvaluateAtContnuousIndex. This has to do with the fact that the
  latter method is const, and I wanted to use the SetLocation method of
  m_Neighborhood. One solution would be to not use neighborhood
  iterators, another to make the method not const in the interpolation
  function hierarchy.

  \par
  The second issue is with the way that the kernel is applied. The
  computational expense comes from two sources: computing the kernel
  weights K(t) and multiplying the pixels in the window by the kernel
  weights. The first is done more or less efficiently in \f$ 2 m d \f$
  operations (where d is the dimensionality of the image). The second
  can be done better. Presently, each pixel \f$ I(i,j,k) \f$ is 
  multiplied by the weights \f$ K(x-i), K(y-j), K(z-k) \f$ and added to
  the running total. This results in \f$ d (2m)^d \f$ multiplication
  operations. However, by keeping intermediate sums, it would be
  possible to do the operation in \f$ O ( (2m)^d ) \f$ operations. This
  would require some creative coding. In addition, in the case when one
  of the coordinates is integer, the computation could be reduced by an
  order of magnitude.
 
  \sa LinearInterpolateImageFunction ResampleImageFilter
  \sa Function::HammingWindowFunction 
  \sa Function::CosineWindowFunction 
  \sa Function::WelchWindowFunction
  \sa Function::LancosWindowFunction 
  \sa Function::BlackmanWindowFunction
  \ingroup ImageFunctions ImageInterpolators
 */
template <
  class TInputImage, 
  unsigned int VRadius, 
  class TWindowFunction = Function::HammingWindowFunction<VRadius>,
  class TBoundaryCondition = ConstantBoundaryCondition<TInputImage>,
  class TCoordRep=float >
class ITK_EXPORT WindowedSincInterpolateImageFunction : 
  public InterpolateImageFunction<TInputImage, TCoordRep>
{
public:
  /** Standard class typedefs. */
  typedef WindowedSincInterpolateImageFunction Self;
  typedef InterpolateImageFunction<TInputImage,TCoordRep> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(WindowedSincInterpolateImageFunction, 
    InterpolateImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** OutputType typedef support. */
  typedef typename Superclass::OutputType OutputType;

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType InputImageType;

  /** RealType typedef support. */
  typedef typename Superclass::RealType RealType;

  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

  /** Index typedef support. */
  typedef typename Superclass::IndexType IndexType;
  
  /** Image type definition */
  typedef TInputImage ImageType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  virtual void SetInputImage(const ImageType *image);

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the interpolated image intensity at a 
   * specified point position.  Bounds checking is based on the 
   * type of the TBoundaryCondition specified. 
   */
  virtual OutputType EvaluateAtContinuousIndex( 
    const ContinuousIndexType & index ) const;

protected:
  WindowedSincInterpolateImageFunction();
  virtual ~WindowedSincInterpolateImageFunction();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  WindowedSincInterpolateImageFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  // Internal typedefs
  typedef ConstNeighborhoodIterator<
    ImageType, TBoundaryCondition> IteratorType;

  // Constant to store twice the radius
  static const unsigned int m_WindowSize;

  /** Neighborhood iterator used to access the pixels in the window */
  IteratorType m_Neighborhood;

  /** The function object, used to compute window */
  TWindowFunction m_WindowFunction;
  
  /** The offset array, used to keep a list of relevant
   * offsets in the neihborhoodIterator */
  unsigned int *m_OffsetTable;

  /** Size of the offset table */
  unsigned int m_OffsetTableSize;

  /** Index into the weights array for each offset */
  unsigned int **m_WeightOffsetTable;

  /** The sinc function */
  inline double Sinc(double x) const
    { 
    double px = vnl_math::pi * x;
    return (x == 0.0) ? 1.0 : sin(px) / px;
    }
};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkWindowedSincInterpolateImageFunction.txx"
#endif

#endif // _itkWindowedSincInterpolateImageFunction_h
