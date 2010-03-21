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
#ifndef __itkSymmetricImageToImageMetric_txx
#define __itkSymmetricImageToImageMetric_txx

// First, make sure that we include the configuration file.
// This line may be removed once the ThreadSafeTransform gets
// integrated into ITK.
#include "itkConfigure.h"

// Second, redirect to the optimized version if necessary
#ifdef ITK_USE_OPTIMIZED_REGISTRATION_METHODS
#include "itkOptImageToImageMetric.txx"
#else

#include "itkSymmetricImageToImageMetric.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"


namespace itk
{

/**
 * Constructor
 */
template <class TFixedImage, class TMovingImage> 
SymmetricImageToImageMetric<TFixedImage,TMovingImage>
::SymmetricImageToImageMetric()
{
  m_FixedImage    = 0; // has to be provided by the user.
  m_MovingImage   = 0; // has to be provided by the user.
  m_HalfwayImage  = 0; // may to be provided by the user.
  m_Transform     = 0; // has to be provided by the user.
  m_FixedTransform  = 0; // internally computed.
  m_MovingTransform = 0; // internally computed.
  m_Interpolator  = 0; // has to be provided by the user.
  m_FixedInterpolator  = 0; // has to be provided by the user.
  m_GradientImage = 0; // will receive the output of the filter;
  m_UseSymmetric    = true; // use symmetric metric
  m_UseSlaveMetric  = true; // use slave metric
  m_ComputeGradient = true; // metric computes gradient by default
  m_NumberOfPixelsCounted = 0; // initialize to zero
  m_GradientImage = NULL; // computed at initialization
  m_AsymMetric = NULL; 
  m_SubtractMean = false;
}

/**
 * Destructor
 */
template <class TFixedImage, class TMovingImage> 
SymmetricImageToImageMetric<TFixedImage,TMovingImage>
::~SymmetricImageToImageMetric()
{

}

// * Get the cost function value at the given parameters
template <class TFixedImage, class TMovingImage> 
typename SymmetricImageToImageMetric <TFixedImage,TMovingImage>
::MeasureType
SymmetricImageToImageMetric <TFixedImage,TMovingImage>
::GetValueInternalSymmetric( const TransformParametersType & parameters ) const
{

  MeasureType measure = NumericTraits< MeasureType >::Zero;
  if (!strcmp( m_AsymMetric->GetNameOfClass(), "MeanSquaresImageToImageMetric"))
    {
    //std::cerr << "Internal method: " << parameters << std::endl;
    typedef  itk::ImageRegionConstIteratorWithIndex< TFixedImage > FixedIteratorType;
    typedef  itk::ImageRegionConstIteratorWithIndex< HalfwayImageType > HalfwayIteratorType;

    HalfwayIteratorType ti( m_HalfwayImage , m_HalfwayImage->GetLargestPossibleRegion() );

    typename HalfwayImageType::IndexType index;


    this->m_NumberOfPixelsCounted = 0;


    while(!ti.IsAtEnd())
      {

      index = ti.GetIndex();

      InputPointType inputPoint;
      m_HalfwayImage->TransformIndexToPhysicalPoint( index, inputPoint );

      OutputPointType transformedFixedPoint = this->m_FixedTransform->TransformPoint( inputPoint );

      if( this->m_FixedImageMask && !this->m_FixedImageMask->IsInside( transformedFixedPoint ) )
        {
        ++ti;
        continue;
        }

      OutputPointType transformedMovingPoint = this->m_MovingTransform->TransformPoint( inputPoint );

      if( this->m_MovingImageMask && !this->m_MovingImageMask->IsInside( transformedMovingPoint ) )
        {
        ++ti;
        continue;
        }

      if( this->m_Interpolator->IsInsideBuffer( transformedMovingPoint ) &&
          this->m_FixedInterpolator->IsInsideBuffer( transformedFixedPoint ))
        {
        const RealType movingValue  = this->m_Interpolator->Evaluate( transformedMovingPoint );
        const RealType fixedValue   = this->m_FixedInterpolator->Evaluate( transformedFixedPoint );
        this->m_NumberOfPixelsCounted++;
        const RealType diff = movingValue - fixedValue;
        measure += diff * diff;
        }

      ++ti;
      }

    if( !this->m_NumberOfPixelsCounted )
      {
      //itkExceptionMacro(<<"All the points mapped to outside of the moving image");
      std::cerr << "WARNING: All the points mapped to outside of the moving image" << std::endl;
      }
    else
      {
      measure /= this->m_NumberOfPixelsCounted;
      }
    }
  else if (!strcmp( m_AsymMetric->GetNameOfClass(), "NormalizedCorrelationImageToImageMetric"))
    {

    typedef  itk::ImageRegionConstIteratorWithIndex< HalfwayImageType > HalfwayIteratorType;


    HalfwayIteratorType ti( m_HalfwayImage , m_HalfwayImage->GetLargestPossibleRegion() );

    typename HalfwayImageType::IndexType index;


    this->m_NumberOfPixelsCounted = 0;


    typedef  typename NumericTraits< MeasureType >::AccumulateType AccumulateType;

    AccumulateType sff = NumericTraits< AccumulateType >::Zero;
    AccumulateType smm = NumericTraits< AccumulateType >::Zero;
    AccumulateType sfm = NumericTraits< AccumulateType >::Zero;
    AccumulateType sf  = NumericTraits< AccumulateType >::Zero;
    AccumulateType sm  = NumericTraits< AccumulateType >::Zero;

    while(!ti.IsAtEnd())
      {

      index = ti.GetIndex();

      InputPointType inputPoint;
      m_HalfwayImage->TransformIndexToPhysicalPoint( index, inputPoint );

      OutputPointType transformedFixedPoint = this->m_FixedTransform->TransformPoint( inputPoint );

      if( this->m_FixedImageMask && !this->m_FixedImageMask->IsInside( transformedFixedPoint ) )
        {
        ++ti;
        continue;
        }

      OutputPointType transformedMovingPoint = this->m_MovingTransform->TransformPoint( inputPoint );

      if( this->m_MovingImageMask && !this->m_MovingImageMask->IsInside( transformedMovingPoint ) )
        {
        ++ti;
        continue;
        }

      if( this->m_Interpolator->IsInsideBuffer( transformedMovingPoint ) &&
          this->m_FixedInterpolator->IsInsideBuffer( transformedFixedPoint ))
        {
        const RealType movingValue  = this->m_Interpolator->Evaluate( transformedMovingPoint );
        const RealType fixedValue   = this->m_FixedInterpolator->Evaluate( transformedFixedPoint );

        sff += fixedValue  * fixedValue;
        smm += movingValue * movingValue;
        sfm += fixedValue  * movingValue;
        if ( this->m_SubtractMean )
          {
          sf += fixedValue;
          sm += movingValue;
          }
        this->m_NumberOfPixelsCounted++;
        }

      ++ti;
      }

    if ( this->m_SubtractMean && this->m_NumberOfPixelsCounted > 0 )
      {
      sff -= ( sf * sf / this->m_NumberOfPixelsCounted );
      smm -= ( sm * sm / this->m_NumberOfPixelsCounted );
      sfm -= ( sf * sm / this->m_NumberOfPixelsCounted );
      }

    const RealType denom = -1.0 * vcl_sqrt(sff * smm );

    if( this->m_NumberOfPixelsCounted > 0 && denom != 0.0)
      {
      measure = sfm / denom;
      }
    else
      {
      measure = NumericTraits< MeasureType >::Zero;
      }
 
    }
  return measure;

}


// * Get the cost function value at the given parameters
template <class TFixedImage, class TMovingImage> 
typename SymmetricImageToImageMetric <TFixedImage,TMovingImage>
::MeasureType
SymmetricImageToImageMetric <TFixedImage,TMovingImage>
::GetValueInternal( const TransformParametersType & parameters ) const
{

  MeasureType measure = NumericTraits< MeasureType >::Zero;
  if (!strcmp( m_AsymMetric->GetNameOfClass(), "MeanSquaresImageToImageMetric"))
    {
    typedef  itk::ImageRegionConstIteratorWithIndex< TFixedImage > FixedIteratorType;

    FixedIteratorType ti( m_FixedImage , this->GetFixedImageRegion() );

    typename TFixedImage::IndexType index;


    this->m_NumberOfPixelsCounted = 0;


    while(!ti.IsAtEnd())
      {

      index = ti.GetIndex();

      InputPointType inputPoint;
      m_FixedImage->TransformIndexToPhysicalPoint( index, inputPoint );

      if( this->m_FixedImageMask && !this->m_FixedImageMask->IsInside( inputPoint ) )
        {
        ++ti;
        continue;
        }

      OutputPointType transformedPoint = this->m_Transform->TransformPoint( inputPoint );

      if( this->m_MovingImageMask && !this->m_MovingImageMask->IsInside( transformedPoint ) )
        {
        ++ti;
        continue;
        }

      if( this->m_Interpolator->IsInsideBuffer( transformedPoint ) )
        {
        const RealType movingValue  = this->m_Interpolator->Evaluate( transformedPoint );
        const RealType fixedValue   = ti.Get();
        this->m_NumberOfPixelsCounted++;
        const RealType diff = movingValue - fixedValue;
        measure += diff * diff;
        }

      ++ti;
      }

    if( !this->m_NumberOfPixelsCounted )
      {
      itkExceptionMacro(<<"All the points mapped to outside of the moving image");
      }
    else
      {
      measure /= this->m_NumberOfPixelsCounted;
      }
    }
  if (!strcmp( m_AsymMetric->GetNameOfClass(), "NormalizedCorrelationImageToImageMetric"))
    {
    typedef  itk::ImageRegionConstIteratorWithIndex< TFixedImage > FixedIteratorType;

    FixedIteratorType ti( m_FixedImage, this->GetFixedImageRegion() );

    typename TFixedImage::IndexType index;


    this->m_NumberOfPixelsCounted = 0;


    typedef  typename NumericTraits< MeasureType >::AccumulateType AccumulateType;

    AccumulateType sff = NumericTraits< AccumulateType >::Zero;
    AccumulateType smm = NumericTraits< AccumulateType >::Zero;
    AccumulateType sfm = NumericTraits< AccumulateType >::Zero;
    AccumulateType sf  = NumericTraits< AccumulateType >::Zero;
    AccumulateType sm  = NumericTraits< AccumulateType >::Zero;

    while(!ti.IsAtEnd())
      {

      index = ti.GetIndex();

      InputPointType inputPoint;
      m_FixedImage->TransformIndexToPhysicalPoint( index, inputPoint );

      if( this->m_FixedImageMask && !this->m_FixedImageMask->IsInside( inputPoint ) )
        {
        ++ti;
        continue;
        }

      OutputPointType transformedPoint = this->m_Transform->TransformPoint( inputPoint );

      if( this->m_MovingImageMask && !this->m_MovingImageMask->IsInside( transformedPoint ) )
        {
        ++ti;
        continue;
        }

      if( this->m_Interpolator->IsInsideBuffer( transformedPoint ) )
        {
        const RealType movingValue  = this->m_Interpolator->Evaluate( transformedPoint );
        const RealType fixedValue   = ti.Get();
        sff += fixedValue  * fixedValue;
        smm += movingValue * movingValue;
        sfm += fixedValue  * movingValue;
        if ( this->m_SubtractMean )
          {
          sf += fixedValue;
          sm += movingValue;
          }
        this->m_NumberOfPixelsCounted++;
        }

      ++ti;
      }
    if ( this->m_SubtractMean && this->m_NumberOfPixelsCounted > 0 )
      {
      sff -= ( sf * sf / this->m_NumberOfPixelsCounted );
      smm -= ( sm * sm / this->m_NumberOfPixelsCounted );
      sfm -= ( sf * sm / this->m_NumberOfPixelsCounted );
      }

    const RealType denom = -1.0 * vcl_sqrt(sff * smm );

    if( this->m_NumberOfPixelsCounted > 0 && denom != 0.0)
      {
      measure = sfm / denom;
      }
    else
      {
      measure = NumericTraits< MeasureType >::Zero;
      }

    }

  return measure;

}

// * Get the cost function value at the given parameters
template <class TFixedImage, class TMovingImage> 
typename SymmetricImageToImageMetric <TFixedImage,TMovingImage>
::MeasureType
SymmetricImageToImageMetric <TFixedImage,TMovingImage>
::GetValue( const TransformParametersType & parameters ) const
{
  itkDebugMacro("Computing Cost Function value at " <<  parameters);

  if(!m_AsymMetric)
    {
    ExceptionObject ex;
    ex.SetLocation(__FILE__);
    ex.SetDescription("The cost function must be set prior to calling GetValue");
    throw ex;
    }

  this->SetTransformParameters( parameters );
  MeasureType metric;

  if ( !m_UseSymmetric )
    {

    if ( m_UseSlaveMetric )
      {
      //m_AsymMetric->SetTransformParameters( parameters );
      m_AsymMetric->SetTransform( m_Transform );
      m_AsymMetric->SetFixedImage( m_FixedImage );
      m_AsymMetric->SetMovingImage( m_MovingImage );
      m_AsymMetric->SetFixedImageRegion(
        m_FixedImage->GetLargestPossibleRegion() );
      if ( m_FixedImageMask )
        m_AsymMetric->SetFixedImageMask(
          m_FixedImageMask );
      if ( m_MovingImageMask )
        m_AsymMetric->SetMovingImageMask(
          m_MovingImageMask );
      m_AsymMetric->SetInterpolator( m_Interpolator );
      // TODO Not all metrics require initialization -- take this out
      m_AsymMetric->Initialize();
      metric = m_AsymMetric->GetValue( parameters );
      }
    else
      {
      metric =  this->GetValueInternal( parameters );    
      }
    
    
    }
  else
    {
  
    // This is the callback function, so everytime we are called,
    // we need to construct the transform and then split it in half


    itk::Matrix<double, FixedImageDimension + 1, FixedImageDimension +1 > amat;
    itk::Matrix<double, FixedImageDimension, FixedImageDimension> tmat = m_Transform->GetMatrix();
    itk::Vector<double, FixedImageDimension> toff = m_Transform->GetOffset();
    Flip_LPS_to_RAS( amat, tmat, toff);


    // Peform Denman-Beavers iteration to compute half and half inverse transforms
    MatrixType Y, Z;
    Y = amat.GetVnlMatrix();
    Z.set_identity();

    for(size_t i = 0; i < 16; i++)
      {
      MatrixType Ynext = 0.5 * (Y + vnl_inverse(Z));
      MatrixType Znext = 0.5 * (Z + vnl_inverse(Y));
      Y = Ynext;
      Z = Znext;
      }

    MatrixType Yinv = vnl_inverse( Y );



    m_FixedTransform = TransformType::New();
    m_MovingTransform = TransformType::New();
  
    this->GetHalfTransform( Y, m_MovingTransform);
    this->GetHalfTransform( Yinv, m_FixedTransform);
  
  
    itkDebugMacro(  "Callback for parameters: " << parameters );
    itkDebugMacro( "Fixed half transform parameters are: " << m_FixedTransform->GetParameters() );
    itkDebugMacro( "Moving half transform parameters are: " << m_MovingTransform->GetParameters() );

    if ( m_UseSlaveMetric )
      {
      typedef typename itk::ResampleImageFilter< TMovingImage, HalfwayImageType >    ResampleMovingFilterType;
      typedef typename itk::ResampleImageFilter< TFixedImage, HalfwayImageType >    ResampleFixedFilterType;
      typename ResampleMovingFilterType::Pointer resamplemoving = ResampleMovingFilterType::New();
      typename ResampleFixedFilterType::Pointer resamplefixed = ResampleFixedFilterType::New();

      resamplemoving->SetInput( m_MovingImage );
      resamplefixed->SetInput( m_FixedImage );
      resamplefixed->SetTransform( m_FixedTransform );
      resamplemoving->SetTransform( m_MovingTransform );
      resamplefixed->SetDefaultPixelValue( 0.0 );
      resamplemoving->SetDefaultPixelValue( 0.0 );

      resamplefixed->UseReferenceImageOn();
      resamplefixed->SetReferenceImage(m_HalfwayImage);
      resamplemoving->UseReferenceImageOn();
      resamplemoving->SetReferenceImage(m_HalfwayImage);

      FixedImagePointer fixedhalfim = resamplefixed->GetOutput();  
      MovingImagePointer movinghalfim = resamplemoving->GetOutput();  

      resamplefixed->Update();
      resamplemoving->Update();

      //TransformPointer idtran = TransformType::New();
      //idtran->SetIdentity();
      //m_AsymMetric->SetTransform (idtran);
      m_AsymMetric->SetFixedImage( fixedhalfim );
      m_AsymMetric->SetMovingImage( movinghalfim );
      m_AsymMetric->SetFixedImageRegion(
         fixedhalfim->GetBufferedRegion() );
      m_AsymMetric->SetInterpolator( m_Interpolator );
      TransformPointer idtran = TransformType::New();
      idtran->SetIdentity();
      m_AsymMetric->SetTransform (idtran);

  
      // TODO Not all metrics require initialization -- take this out
      m_AsymMetric->Initialize();
      metric = m_AsymMetric->GetValue( parameters );
      }
    else
      {
      metric = this->GetValueInternalSymmetric( parameters );
      }
  }
  if (this->GetDebug())
    std::cout << "metric " << metric << " parameters: " << parameters << std::endl; 
    //TransformParametersType finalparam(6) ;
    //finalparam[0] = 0.0289415;
    //finalparam[1] = 0.00178908;
    //finalparam[2] = -0.0275856;
    //finalparam[3] = -5.94439;
    //finalparam[4] = -14.5248;
    //finalparam[5] = 6.53772;
    //if ( finalparam[0]  parameters[0] && finalparam[1] == parameters[1] && finalparam[2] == parameters[2] && finalparam[3] == parameters[3] && finalparam[4] == parameters[4] && finalparam[5] == parameters[5] )
/*
    if (fabs(metric -  -1.16884) < 0.0001)
    { 
    // Debug -- what is the resampled image here ?
std::cout << "This is the optimal parameter" << std::endl;
    typedef typename itk::ResampleImageFilter< TMovingImage, TFixedImage >    ResampleFilterType;
    typename ResampleFilterType::Pointer resample = ResampleFilterType::New();
    itk::Matrix<double, FixedImageDimension + 1, FixedImageDimension +1 > amat;
    itk::Matrix<double, FixedImageDimension, FixedImageDimension> tmat = m_Transform->GetMatrix();
    itk::Vector<double, FixedImageDimension> toff = m_Transform->GetOffset();
    TransformPointer tran = TransformType::New();
    tran->SetMatrix( tmat );
    tran->SetOffset( toff );
    resample->SetInput( m_MovingImage );
    resample->SetTransform( tran );
    resample->SetDefaultPixelValue( 0.0 );

    resample->UseReferenceImageOn();
    resample->SetReferenceImage(m_FixedImage);

    resample->Update();
    MovingImagePointer movingtranim = resample->GetOutput();  
    typedef typename itk::ImageFileWriter< TFixedImage > WriterType;
    typename WriterType::Pointer      writer =  WriterType::New();

    writer->SetFileName( "resampled.nii.gz" );
    writer->SetInput( resample->GetOutput() );
    writer->Update();
    }
*/
  return metric;
}

template <class TFixedImage, class TMovingImage>
void
SymmetricImageToImageMetric<TFixedImage,TMovingImage>
::GetHalfTransform( MatrixType  m, TransformPointer atran ) const
{

  itk::Matrix<double,FixedImageDimension + 1,FixedImageDimension + 1> matrix;
  itk::Matrix<double,FixedImageDimension,FixedImageDimension> amat;
  itk::Vector<double, FixedImageDimension> aoff;

  matrix.GetVnlMatrix().update( m , 0, 0);  
  itkDebugMacro( "Halfinv right before flipping: " );
  if (this->GetDebug())
    this->PrintMatrix( matrix.GetVnlMatrix(), "%12.5f ", "    ");
  Flip_RAS_to_LPS(matrix, amat, aoff);
  atran->SetMatrix(amat);
  atran->SetOffset(aoff);


}

/**
 * Set the parameters that define a unique transform
 */
template <class TFixedImage, class TMovingImage> 
void
SymmetricImageToImageMetric<TFixedImage,TMovingImage>
::SetTransformParameters( const ParametersType & parameters ) const
{
  if( !m_Transform )
    {
    itkExceptionMacro(<<"Transform has not been assigned");
    }
  m_Transform->SetParameters( parameters );
}




/**
 * Initialize
 */
template <class TFixedImage, class TMovingImage> 
void
SymmetricImageToImageMetric<TFixedImage,TMovingImage>
::Initialize(void) throw ( ExceptionObject )
{

  if( !m_Transform )
    {
    itkExceptionMacro(<<"Transform is not present");
    }

  if( !m_AsymMetric )
    {
    itkExceptionMacro(<<"Metric is not present");
    }

  if( !m_Interpolator )
    {
    itkExceptionMacro(<<"Interpolator is not present");
    }

  if( !m_MovingImage )
    {
    itkExceptionMacro(<<"MovingImage is not present");
    }

  if( !m_FixedImage )
    {
    itkExceptionMacro(<<"FixedImage is not present");
    }

  if( m_FixedImageRegion.GetNumberOfPixels() == 0 )
    {
    itkExceptionMacro(<<"FixedImageRegion is empty");
    }

  // If the image is provided by a source, update the source.
  if( m_MovingImage->GetSource() )
    {
    m_MovingImage->GetSource()->Update();
    }

  // If the image is provided by a source, update the source.
  if( m_FixedImage->GetSource() )
    {
    m_FixedImage->GetSource()->Update();
    }

  // Make sure the FixedImageRegion is within the FixedImage buffered region
  if ( !m_FixedImageRegion.Crop( m_FixedImage->GetBufferedRegion() ) )
    {
    itkExceptionMacro(
      <<"FixedImageRegion does not overlap the fixed image buffered region" );
    }

  m_Interpolator->SetInputImage( m_MovingImage );

  if( !m_HalfwayImage )
    {
    //if (m_UseSymmetric)
      {
      this->CreateHalfwayImageSpace();
      itkDebugMacro( "Halfway image space created" );
      }
    }

  // Initialize the slave metric
  TransformPointer idtran = TransformType::New();
  idtran->SetIdentity();
  if (this->GetDebug())
    m_AsymMetric->DebugOn();
  if (m_UseSymmetric)
    {
    m_AsymMetric->SetTransform (idtran);
    if ( !m_UseSlaveMetric )
      if( !m_FixedInterpolator )
      {
      itkExceptionMacro(<<"Interpolator is not present");
      }
      m_FixedInterpolator->SetInputImage( m_FixedImage ); 
    }
  else
    {
    m_AsymMetric->SetTransform( this->GetTransform() );
    }
  m_AsymMetric->SetFixedImage( m_FixedImage );
  m_AsymMetric->SetMovingImage( m_MovingImage );
  m_AsymMetric->SetFixedImageRegion(
       m_FixedImage->GetBufferedRegion() );

  if ( m_InternalFixedImageMask )
    {
    itkDebugMacro( "Setting fixed image mask from internal mask" );
    if ( m_FixedImageMask )
      m_AsymMetric->SetFixedImageMask( m_FixedImageMask );
    else
      {
      typename InternalFixedImageBinaryMaskType::Pointer fixedimageMask = InternalFixedImageBinaryMaskType::New();
      typename InternalFixedImageMaskType::RegionType fixedRegion;
      fixedRegion.SetSize( m_InternalFixedImageMask->GetLargestPossibleRegion().GetSize() );
      fixedimageMask->SetRegions( fixedRegion );
      fixedimageMask->Allocate();
      fixedimageMask->SetOrigin(  m_InternalFixedImageMask->GetOrigin() );
      fixedimageMask->SetSpacing( m_InternalFixedImageMask->GetSpacing() );
      fixedimageMask->SetDirection( m_InternalFixedImageMask->GetDirection() );
      fixedimageMask->FillBuffer(1);
      typedef  typename itk::ImageRegionIteratorWithIndex<InternalFixedImageMaskType> InternalFixedImageMaskIteratorType;
      typedef  typename itk::ImageRegionIteratorWithIndex<InternalFixedImageBinaryMaskType> InternalFixedImageBinaryMaskIteratorType;
      InternalFixedImageMaskIteratorType itsource( m_InternalFixedImageMask, m_InternalFixedImageMask->GetLargestPossibleRegion() );
      InternalFixedImageBinaryMaskIteratorType ittarget( fixedimageMask, fixedimageMask->GetLargestPossibleRegion() );
      while(!itsource.IsAtEnd())
        {
        if (itsource.Get() < 0.5)
            ittarget.Set( 0 );
        ++itsource;
        ++ittarget;
        }
      typename itk::ImageMaskSpatialObject< FixedImageDimension >::Pointer  fixedMask = itk::ImageMaskSpatialObject< FixedImageDimension >::New();
      fixedMask->SetImage( fixedimageMask );
      m_AsymMetric->SetFixedImageMask(
         fixedMask );
      }
    }

  if ( m_InternalMovingImageMask )
    {
    itkDebugMacro( "Setting moving image mask from internal mask" );
    if ( m_MovingImageMask )
      m_AsymMetric->SetMovingImageMask( m_MovingImageMask );
    else
      {
      typename InternalMovingImageBinaryMaskType::Pointer movingimageMask = InternalMovingImageBinaryMaskType::New();
      typename InternalMovingImageMaskType::RegionType movingRegion;
      movingRegion.SetSize( m_InternalMovingImageMask->GetLargestPossibleRegion().GetSize() );
      movingimageMask->SetRegions( movingRegion );
      movingimageMask->Allocate();
      movingimageMask->SetOrigin(  m_InternalMovingImageMask->GetOrigin() );
      movingimageMask->SetSpacing( m_InternalMovingImageMask->GetSpacing() );
      movingimageMask->SetDirection( m_InternalMovingImageMask->GetDirection() );
      movingimageMask->FillBuffer(1);
      typedef  typename itk::ImageRegionIteratorWithIndex<InternalMovingImageMaskType> InternalMovingImageMaskIteratorType;
      typedef  typename itk::ImageRegionIteratorWithIndex<InternalMovingImageBinaryMaskType> InternalMovingImageBinaryMaskIteratorType;
      InternalMovingImageMaskIteratorType itsource( m_InternalMovingImageMask, m_InternalMovingImageMask->GetLargestPossibleRegion() );
      InternalMovingImageBinaryMaskIteratorType ittarget( movingimageMask, movingimageMask->GetLargestPossibleRegion() );
      while(!itsource.IsAtEnd())
        {
        if (itsource.Get() < 0.5)
            ittarget.Set( 0 );
        ++itsource;
        ++ittarget;
        }
      typename itk::ImageMaskSpatialObject< MovingImageDimension >::Pointer  movingMask = itk::ImageMaskSpatialObject< MovingImageDimension >::New();
      movingMask->SetImage( movingimageMask );
      m_AsymMetric->SetMovingImageMask(
         movingMask );
      }
    }



      
  
  m_AsymMetric->SetInterpolator( m_Interpolator );
  m_AsymMetric->Initialize();

 
  if ( m_ComputeGradient )
    {
    this->ComputeGradient();
    }

  // If there are any observers on the metric, call them to give the
  // user code a chance to set parameters on the metric
  this->InvokeEvent( InitializeEvent() );
}

template <class TFixedImage, class TMovingImage> 
void
SymmetricImageToImageMetric<TFixedImage,TMovingImage>
::CreateHalfwayImageSpace(void)
{
//  c3d_affine_tool -sform $TP1 -sform $TP0 -inv -mult -sqrt -sform $TP0 -mult -o $WDIR/hwspace.mat
  MatrixType mfixed  = this->m_FixedImage->GetVoxelSpaceToRASPhysicalSpaceMatrix().GetVnlMatrix();
  MatrixType mmoving = this->m_MovingImage->GetVoxelSpaceToRASPhysicalSpaceMatrix().GetVnlMatrix();
  
  MatrixType mcomb = mmoving * vnl_inverse(mfixed);
  // Peform Denman-Beavers iteration
  MatrixType Z, Y = mcomb;
  Z.set_identity();

  for(size_t i = 0; i < 16; i++)
    {
    MatrixType Ynext = 0.5 * (Y + vnl_inverse(Z));
    MatrixType Znext = 0.5 * (Z + vnl_inverse(Y));
    Y = Ynext;
    Z = Znext;
    }

  MatrixType mhalf = Y * mfixed;


  // Create the image
  FixedImagePointer hwspaceimg = FixedImageType::New();
  hwspaceimg->SetBufferedRegion(this->m_FixedImage->GetBufferedRegion());
  // Set the voxel size
  hwspaceimg->SetSpacing(m_FixedImage->GetSpacing());
  hwspaceimg->Allocate();
  hwspaceimg->FillBuffer( 0.0 );

  // Set the matrix
  hwspaceimg->SetVoxelSpaceToRASPhysicalSpaceMatrix( mhalf );

  // Now apply half transform to fixed image and resample to this space, thus creating halfway reference image
  typedef typename itk::ResampleImageFilter< FixedImageType, HalfwayImageType >    ResampleFixedFilterType;
  typename ResampleFixedFilterType::Pointer resamplefixed = ResampleFixedFilterType::New();

  resamplefixed->SetInput( m_FixedImage );



  itk::Matrix<double, FixedImageDimension + 1, FixedImageDimension +1 > amat;
  itk::Matrix<double, FixedImageDimension, FixedImageDimension> tmat = m_Transform->GetMatrix();
  itk::Vector<double, FixedImageDimension> toff = m_Transform->GetOffset();
  Flip_LPS_to_RAS( amat, tmat, toff);

  itkDebugMacro( "Full transform: " );
  if (this->GetDebug())
    PrintMatrix( amat.GetVnlMatrix(), "%12.5f ", "    ");

  // Peform Denman-Beavers iteration to compute half and half inverse transforms
  Y = amat.GetVnlMatrix();
  Z.set_identity();

  for(size_t i = 0; i < 16; i++)
    {
    MatrixType Ynext = 0.5 * (Y + vnl_inverse(Z));
    MatrixType Znext = 0.5 * (Z + vnl_inverse(Y));
    Y = Ynext;
    Z = Znext;
    }

  MatrixType Yinv = vnl_inverse( Y );

  TransformPointer fixedtran = TransformType::New();

  this->GetHalfTransform( Yinv, fixedtran);
  itkDebugMacro( "Halfinv transform: " );
  if (this->GetDebug())
    PrintMatrix( Yinv, "%12.5f ", "    ");

  resamplefixed->SetTransform( fixedtran );
  resamplefixed->SetDefaultPixelValue( 0.0 );

  resamplefixed->UseReferenceImageOn();
  resamplefixed->SetReferenceImage( hwspaceimg);

  HalfwayImagePointer hwimg = resamplefixed->GetOutput();

  resamplefixed->Update();

  this->SetHalfwayImage( hwimg );

  itkDebugMacro( "Fixed image sform is " );
  if (this->GetDebug())
    PrintMatrix( mfixed, "%12.5f ", "    ");
  itkDebugMacro( "Moving image sform is " );
  if (this->GetDebug())
    PrintMatrix( mmoving, "%12.5f ", "    ");
  itkDebugMacro( "Halfway image sform is " );
  if (this->GetDebug())
    PrintMatrix( m_HalfwayImage->GetVoxelSpaceToRASPhysicalSpaceMatrix().GetVnlMatrix(), "%12.5f ", "    ");
  

}

/*
 * Compute the gradient image and assign it to m_GradientImage.
 */
template <class TFixedImage, class TMovingImage> 
void
SymmetricImageToImageMetric<TFixedImage,TMovingImage>
::ComputeGradient() 
{
  GradientImageFilterPointer gradientFilter
    = GradientImageFilterType::New();

  gradientFilter->SetInput( m_MovingImage );

  const typename MovingImageType::SpacingType&
    spacing = m_MovingImage->GetSpacing();
  double maximumSpacing=0.0;
  for(unsigned int i=0; i<MovingImageDimension; i++)
    {
    if( spacing[i] > maximumSpacing )
      {
      maximumSpacing = spacing[i];
      }
    }
  gradientFilter->SetSigma( maximumSpacing );
  gradientFilter->SetNormalizeAcrossScale( true );

#ifdef ITK_USE_ORIENTED_IMAGE_DIRECTION
  gradientFilter->SetUseImageDirection( true );
#endif
  
  gradientFilter->Update();
  
  m_GradientImage = gradientFilter->GetOutput();
}

template <class TFixedImage, class TMovingImage> 
void
SymmetricImageToImageMetric<TFixedImage,TMovingImage>
::PrintMatrix(vnl_matrix<double> mat, const char *fmt, const char *prefix) const
{
  // Print each row and column of the matrix
  char buffer[256];
  for(size_t i = 0; i < mat.rows(); i++)
    {
    std::cout << prefix;
    for(size_t j = 0; j < mat.columns(); j++)
      {
      sprintf(buffer, fmt, mat(i,j));
      std::cout << buffer;
      }
    std::cout << std::endl;
    }
}

template <class TFixedImage, class TMovingImage> 
void
SymmetricImageToImageMetric<TFixedImage,TMovingImage>
::Flip_LPS_to_RAS(  itk::Matrix<double,FixedImageDimension+1,FixedImageDimension+1> &matrix, 
                    itk::Matrix<double,FixedImageDimension,FixedImageDimension> &amat, 
                    itk::Vector<double, FixedImageDimension> &aoff) const
{   

    // Convert lps to ras
    vnl_vector<double> v_ras_to_lps(FixedImageDimension, 1.0);
    v_ras_to_lps[0] = v_ras_to_lps[1] = -1.0;
    vnl_diag_matrix<double> m_ras_to_lps(v_ras_to_lps);

    vnl_matrix<double> amatvnl = amat.GetVnlMatrix();
    amatvnl = m_ras_to_lps * amatvnl * m_ras_to_lps;
    vnl_vector_fixed<double, FixedImageDimension > aoffs ;
    vnl_vector_fixed<double, FixedImageDimension + 1> aoffl ;
    aoffs = m_ras_to_lps * aoff.GetVnlVector();
    aoffl.fill(1.0);
    for (size_t i=0; i<FixedImageDimension; i++)
      aoffl(i) = aoffs(i);    

    matrix.GetVnlMatrix().set_identity();
    matrix.GetVnlMatrix().update( amatvnl, 0, 0); 
    matrix.GetVnlMatrix().set_column(FixedImageDimension, aoffl);


}


template <class TFixedImage, class TMovingImage> 
void
SymmetricImageToImageMetric<TFixedImage,TMovingImage>
::Flip_RAS_to_LPS(  itk::Matrix<double,FixedImageDimension+1,FixedImageDimension+1> &matrix, 
                    itk::Matrix<double,FixedImageDimension,FixedImageDimension> &amat, 
                    itk::Vector<double, FixedImageDimension> &aoff) const
{

    amat.GetVnlMatrix().update(
      matrix.GetVnlMatrix().extract(FixedImageDimension, FixedImageDimension));
    aoff.GetVnlVector().update(
      matrix.GetVnlMatrix().get_column(FixedImageDimension).extract(FixedImageDimension));


    // External matrices are assumed to be RAS to RAS, so we must convert to LPS to LPS
    vnl_vector<double> v_lps_to_ras(FixedImageDimension, 1.0);
    v_lps_to_ras[0] = v_lps_to_ras[1] = -1.0;
    vnl_diag_matrix<double> m_lps_to_ras(v_lps_to_ras);
    vnl_matrix<double> mold = amat.GetVnlMatrix();
    amat.GetVnlMatrix().update(m_lps_to_ras * mold * m_lps_to_ras);
    aoff.GetVnlVector().update(m_lps_to_ras * aoff.GetVnlVector());


}



/**
 * PrintSelf
 */
template <class TFixedImage, class TMovingImage> 
void
SymmetricImageToImageMetric<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "ComputeGradient: "
     << static_cast<typename NumericTraits<bool>::PrintType>(m_ComputeGradient)
     << std::endl;
  os << indent << "Moving   Image: " << m_MovingImage.GetPointer()  << std::endl;
  os << indent << "Fixed    Image: " << m_FixedImage.GetPointer()   << std::endl;
  os << indent << "Halfway  Image: " << m_FixedImage.GetPointer()   << std::endl;
  os << indent << "Gradient Image: " << m_HalfwayImage.GetPointer() 
     << std::endl;
  os << indent << "Transform:    " << m_Transform.GetPointer()    << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
  os << indent << "Fixed Interpolator: " << m_FixedInterpolator.GetPointer() << std::endl;
  os << indent << "FixedImageRegion: " << m_FixedImageRegion << std::endl;
  os << indent << "Moving Image Mask: " << m_MovingImageMask.GetPointer() 
     << std::endl;
  os << indent << "Fixed Image Mask: " << m_FixedImageMask.GetPointer() 
     << std::endl;
  os << indent << "Moving Image Internal Mask: " << m_InternalMovingImageMask.GetPointer() 
     << std::endl;
  os << indent << "Fixed Image Internal Mask: " << m_InternalFixedImageMask.GetPointer() 
     << std::endl;
  os << indent << "Number of Pixels Counted: " << m_NumberOfPixelsCounted 
     << std::endl;

}


} // end namespace itk

#endif

#endif
