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


namespace itk
{

/**
 * Connect a Cost Function
 */
void
SymmetricImageToImageMetric
::SetCostFunction( CostFunctionType * costFunction )
{
  if( m_CostFunction.GetPointer() == costFunction )
    {
    return;
    }

  itkDebugMacro("setting CostFunction  to " <<  costFunction);

  m_CostFunction = costFunction;

  if(!m_ScalesInitialized)
    {
    const unsigned int numberOfParameters =
      m_CostFunction->GetNumberOfParameters();

    ScalesType scales( numberOfParameters );
    scales.Fill( 1.0f );
    SetScales( scales );
    m_ScalesInitialized = true;
    }

  this->Modified();
}

/**
 * Get the cost function value at the given parameters
 */
SymmetricImageToImageMetric::MeasureType
SymmetricImageToImageMetric
::GetValue( const ParametersType & parameters ) const
{
  itkDebugMacro("Computing CostFunction value at " <<  parameters);

  if(!m_CostFunction)
    {
    ExceptionObject ex;
    ex.SetLocation(__FILE__);
    ex.SetDescription("The costfunction must be set prior to calling GetValue");
    throw ex;
    }
/**************** code for symmetrization goes here *********************/

  return this->GetCostFunction()->GetValue(parameters);
}


/**
 * Constructor
 */
template <class TFixedImage, class TMovingImage> 
SymmetricImageToImageMetric<TFixedImage,TMovingImage>
::SymmetricImageToImageMetric()
{
  m_FixedImage    = 0; // has to be provided by the user.
  m_MovingImage   = 0; // has to be provided by the user.
  m_Transform     = 0; // has to be provided by the user.
  m_Interpolator  = 0; // has to be provided by the user.
  m_GradientImage = 0; // will receive the output of the filter;
  m_ComputeGradient = true; // metric computes gradient by default
  m_NumberOfPixelsCounted = 0; // initialize to zero
  m_GradientImage = NULL; // computed at initialization
}

/**
 * Destructor
 */
template <class TFixedImage, class TMovingImage> 
SymmetricImageToImageMetric<TFixedImage,TMovingImage>
::~SymmetricImageToImageMetric()
{

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
 
  if ( m_ComputeGradient )
    {
    this->ComputeGradient();
    }

  // If there are any observers on the metric, call them to give the
  // user code a chance to set parameters on the metric
  this->InvokeEvent( InitializeEvent() );
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
  os << indent << "Moving Image: " << m_MovingImage.GetPointer()  << std::endl;
  os << indent << "Fixed  Image: " << m_FixedImage.GetPointer()   << std::endl;
  os << indent << "Gradient Image: " << m_GradientImage.GetPointer() 
     << std::endl;
  os << indent << "Transform:    " << m_Transform.GetPointer()    << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
  os << indent << "FixedImageRegion: " << m_FixedImageRegion << std::endl;
  os << indent << "Moving Image Mask: " << m_MovingImageMask.GetPointer() 
     << std::endl;
  os << indent << "Fixed Image Mask: " << m_FixedImageMask.GetPointer() 
     << std::endl;
  os << indent << "Number of Pixels Counted: " << m_NumberOfPixelsCounted 
     << std::endl;

}


} // end namespace itk

#endif

#endif
