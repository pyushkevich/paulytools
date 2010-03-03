/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSymmetricImageToImageMetric_h
#define __itkSymmetricImageToImageMetric_h

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"

#ifdef ITK_USE_OPTIMIZED_REGISTRATION_METHODS
#include "itkOptImageToImageMetric.h"
#else

#include "itkImageBase.h"
#include "itkTransform.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkInterpolateImageFunction.h"
#include "itkSingleValuedCostFunction.h"
#include "itkResampleImageFilter.h"
#include "itkExceptionObject.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkSpatialObject.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/vnl_det.h"
#include "vnl/vnl_inverse.h"


namespace itk
{
  
/** \class SymmetricImageToImageMetric
 * \brief Computes similarity between regions of two images.
 *
 * This Class is templated over the type of the two input images.
 * It expects a Transform and an Interpolator to be plugged in.
 * This particular class is the base class for a hierarchy of 
 * similarity metrics.
 *
 * This class computes a value that measures the similarity 
 * between the Fixed image and the transformed Moving image.
 * The Interpolator is used to compute intensity values on 
 * non-grid positions resulting from mapping points through 
 * the Transform.
 * 
 *
 * \ingroup RegistrationMetrics
 *
 */

template <class TFixedImage,  class TMovingImage> 
class ITK_EXPORT SymmetricImageToImageMetric : public SingleValuedCostFunction 
{
public:
  /** Standard class typedefs. */
  typedef SymmetricImageToImageMetric           Self;
  typedef SingleValuedCostFunction     Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Type used for representing point components  */
  typedef typename Superclass::ParametersValueType CoordinateRepresentationType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);


  /** Run-time type information (and related methods). */
  itkTypeMacro(SymmetricImageToImageMetric, SingleValuedCostFunction);

  /**  Type of the moving Image. */
  typedef TMovingImage                               MovingImageType;
  typedef typename TMovingImage::PixelType           MovingImagePixelType;
  typedef typename MovingImageType::Pointer     MovingImagePointer;

  /**  Type of the fixed Image. */
  typedef TFixedImage                                FixedImageType;
  typedef typename FixedImageType::Pointer      FixedImagePointer;
  typedef typename FixedImageType::RegionType        FixedImageRegionType;

  /**  Type of the halfway Image. */
  typedef TFixedImage                                HalfwayImageType;
  typedef typename HalfwayImageType::Pointer      HalfwayImagePointer;
  typedef typename HalfwayImageType::RegionType        HalfwayImageRegionType;

  /** Matrix type for RAS transform */
  typedef vnl_matrix_fixed<double, TFixedImage::ImageDimension + 1, TFixedImage::ImageDimension + 1> MatrixType;

  /** Constants for the image dimensions */
  itkStaticConstMacro(MovingImageDimension, 
                      unsigned int,
                      TMovingImage::ImageDimension);
  itkStaticConstMacro(FixedImageDimension, 
                     unsigned int,
                      TFixedImage::ImageDimension);
  
  /**  Type of the Transform Base class */
  //typedef itk::MatrixOffsetTransformBase<CoordinateRepresentationType, 
   //                 itkGetStaticConstMacro(MovingImageDimension),
    //                itkGetStaticConstMacro(FixedImageDimension)> 
     //                                                TransformType;
  typedef itk::MatrixOffsetTransformBase<CoordinateRepresentationType, 
                    MovingImageDimension,
                    FixedImageDimension> 
                                                     TransformType;
  

  typedef typename TransformType::Pointer            TransformPointer;
  typedef typename TransformType::InputPointType     InputPointType;
  typedef typename TransformType::OutputPointType    OutputPointType;
  typedef typename TransformType::ParametersType     TransformParametersType;
  typedef typename TransformType::JacobianType       TransformJacobianType;

  /**  Type of the Interpolator Base class */
  typedef InterpolateImageFunction<
    MovingImageType,
    CoordinateRepresentationType > InterpolatorType;


  /** Gaussian filter to compute the gradient of the Moving Image */
  typedef typename NumericTraits<MovingImagePixelType>::RealType 
                                                     RealType;
  typedef CovariantVector<RealType,
                          //itkGetStaticConstMacro(MovingImageDimension)>
                          MovingImageDimension>
                                                     GradientPixelType;
  typedef Image<GradientPixelType,
                //itkGetStaticConstMacro(MovingImageDimension)> 
                MovingImageDimension> 
                                                     GradientImageType;
  typedef SmartPointer<GradientImageType>            GradientImagePointer;
  typedef GradientRecursiveGaussianImageFilter< MovingImageType,
                                                GradientImageType >
                                                     GradientImageFilterType;  
  typedef typename GradientImageFilterType::Pointer  GradientImageFilterPointer;


  typedef typename InterpolatorType::Pointer         InterpolatorPointer;


  /**  Type for the mask of the fixed image. Only pixels that are "inside"
       this mask will be considered for the computation of the metric */
  //typedef SpatialObject< itkGetStaticConstMacro(FixedImageDimension) >
  typedef SpatialObject< FixedImageDimension >
                                                     FixedImageMaskType;
  typedef typename FixedImageMaskType::Pointer       FixedImageMaskPointer;
  typedef typename FixedImageMaskType::ConstPointer  FixedImageMaskConstPointer;


  /**  Type for the mask of the moving image. Only pixels that are "inside"
       this mask will be considered for the computation of the metric */
 // typedef SpatialObject< itkGetStaticConstMacro(MovingImageDimension) >
  typedef SpatialObject< MovingImageDimension >
                                                     MovingImageMaskType;
  typedef typename MovingImageMaskType::Pointer      MovingImageMaskPointer;
  typedef typename MovingImageMaskType::ConstPointer MovingImageMaskConstPointer;


  /**  Type of the measure. */
  typedef typename Superclass::MeasureType                    MeasureType;

  /**  Type of the derivative. */
  typedef typename Superclass::DerivativeType                 DerivativeType;

  /**  Type of the parameters. */
  typedef typename Superclass::ParametersType                 ParametersType;

  /** Connect the Fixed Image.  */
  itkSetObjectMacro( FixedImage, FixedImageType );

  /** Get the Fixed Image. */
  itkGetObjectMacro( FixedImage, FixedImageType );

  /** Connect the Moving Image.  */
  itkSetObjectMacro( MovingImage, MovingImageType );

  /** Get the Moving Image. */
  itkGetObjectMacro( MovingImage, MovingImageType );

  /** Connect the Fixed Image.  */
  itkSetObjectMacro( HalfwayImage, HalfwayImageType );

  /** Get the Fixed Image. */
  itkGetObjectMacro( HalfwayImage, HalfwayImageType );

  /** Connect the Transform. */
  itkSetObjectMacro( Transform, TransformType );

  /** Get a pointer to the Transform.  */
  itkGetObjectMacro( Transform, TransformType );
 
  /** Connect the Interpolator. */
  itkSetObjectMacro( Interpolator, InterpolatorType );

  /** Get a pointer to the Interpolator.  */
  itkGetConstObjectMacro( Interpolator, InterpolatorType );

  /** Get the number of pixels considered in the computation. */
  itkGetConstReferenceMacro( NumberOfPixelsCounted, unsigned long );

  /** Set the region over which the metric will be computed */
  itkSetMacro( FixedImageRegion, FixedImageRegionType );

  /** Get the region over which the metric will be computed */
  itkGetConstReferenceMacro( FixedImageRegion, FixedImageRegionType );
 
  /** Set/Get the moving image mask. */
  itkSetObjectMacro( MovingImageMask, MovingImageMaskType );
#ifdef ITK_LEGACY_REMOVE
  itkSetConstObjectMacro( MovingImageMask, MovingImageMaskType );
#else
  virtual void SetMovingImageMask( const MovingImageMaskType* mask )
    { this->SetMovingImageMask(const_cast<MovingImageMaskType*>(mask)); }
#endif
  itkGetConstObjectMacro( MovingImageMask, MovingImageMaskType );

  /** Set/Get the fixed image mask. */
  itkSetObjectMacro( FixedImageMask, FixedImageMaskType );
#ifdef ITK_LEGACY_REMOVE
  itkSetConstObjectMacro( FixedImageMask, FixedImageMaskType );
#else
  virtual void SetFixedImageMask( const FixedImageMaskType* mask )
    { this->SetFixedImageMask(const_cast<FixedImageMaskType*>(mask)); }
#endif
  itkGetConstObjectMacro( FixedImageMask, FixedImageMaskType );

  /** Set/Get gradient computation. */
  itkSetMacro( ComputeGradient, bool);
  itkGetConstReferenceMacro( ComputeGradient, bool);
  itkBooleanMacro(ComputeGradient);

  /** Computes the gradient image and assigns it to m_GradientImage */
  virtual void ComputeGradient();

  /** Get Gradient Image. */
  itkGetConstObjectMacro( GradientImage, GradientImageType );

  // Stuff regarding the asymmetric metric
    /**  Type of the metric. */
  typedef itk::ImageToImageMetric<FixedImageType, MovingImageType>  AsymmetricMetricType;
  typedef typename AsymmetricMetricType::Pointer                AsymmetricMetricPointer;

  /** Set/Get the Metric. */
  itkSetObjectMacro( AsymMetric, AsymmetricMetricType );
  itkGetObjectMacro( AsymMetric, AsymmetricMetricType );
  typedef typename AsymmetricMetricType::ParametersType AsymParametersType;

  /** Set the parameters defining the Transform. */
  void SetTransformParameters( const ParametersType & parameters ) const;

  /** Get the parameters defining the Transform. */
  ParametersType GetTransformParameters()
  { return m_Transform->GetParameters(); }

  /** Return the number of parameters required by the Transform */
  unsigned int GetNumberOfParameters(void) const 
    { return m_Transform->GetNumberOfParameters(); }

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  void Initialize(void) throw ( ExceptionObject );

  // Stuff related to the symmetric metric
  /**  Get the value for single valued optimizers. */

  MeasureType GetValue( const TransformParametersType & parameters ) const;

  //void GetHalfTransform( void );
  void GetHalfTransform( MatrixType  Z, TransformPointer atran ) const;
  void CreateHalfwayImageSpace(void);


  /** This method returns the derivative of the cost function corresponding
    * to the specified parameters.   */
  virtual void GetDerivative( const ParametersType & parameters,
                              DerivativeType & derivative ) const
    { derivative = 0;}

  /** This method returns the value and derivative of the cost function corresponding
    * to the specified parameters    */
  virtual void GetValueAndDerivative( const ParametersType & parameters,
                                      MeasureType & value,
                                      DerivativeType & derivative ) const
    {
    value = this->GetValue( parameters );
    this->GetDerivative( parameters, derivative );
    }




protected:
  SymmetricImageToImageMetric();
  virtual ~SymmetricImageToImageMetric();
  void PrintSelf(std::ostream& os, Indent indent) const;

  void PrintMatrix(vnl_matrix<double> mat, const char *fmt, const char *prefix);

  mutable unsigned long       m_NumberOfPixelsCounted;

  FixedImagePointer      m_FixedImage;
  MovingImagePointer     m_MovingImage;

  HalfwayImagePointer         m_HalfwayImage;

  mutable TransformPointer    m_Transform;
  InterpolatorPointer         m_Interpolator;

  bool                        m_ComputeGradient;
  GradientImagePointer        m_GradientImage;

#ifdef ITK_LEGACY_REMOVE
  FixedImageMaskConstPointer  m_FixedImageMask;
  MovingImageMaskConstPointer m_MovingImageMask;
#else
  mutable FixedImageMaskPointer   m_FixedImageMask;
  mutable MovingImageMaskPointer  m_MovingImageMask;
#endif


  // Stuff regarding the asymmetric metric
  mutable AsymmetricMetricPointer                    m_AsymMetric;

private:
  SymmetricImageToImageMetric(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  FixedImageRegionType        m_FixedImageRegion;  

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSymmetricImageToImageMetric.txx"
#endif

#endif

#endif
