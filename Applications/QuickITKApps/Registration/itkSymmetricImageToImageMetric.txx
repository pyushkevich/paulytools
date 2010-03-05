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
 * Constructor
 */
template <class TFixedImage, class TMovingImage> 
SymmetricImageToImageMetric<TFixedImage,TMovingImage>
::SymmetricImageToImageMetric()
{
  m_FixedImage    = 0; // has to be provided by the user.
  m_MovingImage   = 0; // has to be provided by the user.
  m_HalfwayImage  = 0; // has to be provided by the user.
  m_Transform     = 0; // has to be provided by the user.
  m_Interpolator  = 0; // has to be provided by the user.
  m_GradientImage = 0; // will receive the output of the filter;
  m_UseSymmetric     = true; // use symmetric metric
  m_ComputeGradient = true; // metric computes gradient by default
  m_NumberOfPixelsCounted = 0; // initialize to zero
  m_GradientImage = NULL; // computed at initialization
  m_AsymMetric = NULL; 
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
::GetValue( const TransformParametersType & parameters ) const
{
  itkDebugMacro("Computing CostFunction value at " <<  parameters);

  if(!m_AsymMetric)
    {
    ExceptionObject ex;
    ex.SetLocation(__FILE__);
    ex.SetDescription("The costfunction must be set prior to calling GetValue");
    throw ex;
    }

  this->SetTransformParameters( parameters );

  if (!m_UseSymmetric)
    {
    m_AsymMetric->SetTransformParameters( parameters );
    }
  else
    {
    typedef typename itk::ResampleImageFilter< TMovingImage, HalfwayImageType >    ResampleMovingFilterType;
    typedef typename itk::ResampleImageFilter< TFixedImage, HalfwayImageType >    ResampleFixedFilterType;
    typename ResampleMovingFilterType::Pointer resamplemoving = ResampleMovingFilterType::New();
    typename ResampleFixedFilterType::Pointer resamplefixed = ResampleFixedFilterType::New();

    resamplemoving->SetInput( m_MovingImage );
    resamplefixed->SetInput( m_FixedImage );
  
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



    TransformPointer fixedtran = TransformType::New();
    TransformPointer movingtran = TransformType::New();
    //fixedtran.SetIdentity();
    //movingtran.SetIdentity();
  
    this->GetHalfTransform( Y, movingtran);
    this->GetHalfTransform( Yinv, fixedtran);
  
  
    itkDebugMacro(  "Callback for parameters: " << parameters );
    itkDebugMacro( "Fixed half transform parameters are: " << fixedtran->GetParameters() );
    itkDebugMacro( "Moving half transform parameters are: " << movingtran->GetParameters() );

    resamplefixed->SetTransform( fixedtran );
    resamplemoving->SetTransform( movingtran );
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
         fixedhalfim->GetLargestPossibleRegion() );
  
    // MattesMutualInformationImageToImageMetric requires initialization every time
    if (!strcmp( m_AsymMetric->GetNameOfClass(), "MattesMutualInformationImageToImageMetric"))
      {
      m_AsymMetric->Initialize();
      }
  }
  MeasureType metric = m_AsymMetric->GetValue( m_AsymMetric->GetTransform()->GetParameters() );
  //if (this->GetDebug())
    std::cout << "metric " << metric << " parameters: " << parameters << std::endl; 
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
  if this->GetDebug()
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
    m_AsymMetric->SetTransform (idtran);
  else
    m_AsymMetric->SetTransform( this->GetTransform() );
  m_AsymMetric->SetFixedImage( m_FixedImage );
  m_AsymMetric->SetMovingImage( m_MovingImage );
  m_AsymMetric->SetFixedImageRegion(
       m_FixedImage->GetLargestPossibleRegion() );
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
  hwspaceimg->SetLargestPossibleRegion(this->m_FixedImage->GetLargestPossibleRegion());
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
::Flip_LPS_to_RAS(itk::Matrix<double,FixedImageDimension+1,FixedImageDimension+1> &matrix, itk::Matrix<double,FixedImageDimension,FixedImageDimension> &amat, itk::Vector<double, FixedImageDimension> &aoff) const
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
::Flip_RAS_to_LPS(itk::Matrix<double,FixedImageDimension+1,FixedImageDimension+1> &matrix, itk::Matrix<double,FixedImageDimension,FixedImageDimension> &amat, itk::Vector<double, FixedImageDimension> &aoff) const
{

    amat.GetVnlMatrix().update(
      matrix.GetVnlMatrix().extract(FixedImageDimension, FixedImageDimension));
    aoff.GetVnlVector().update(
      matrix.GetVnlMatrix().get_column(FixedImageDimension).extract(FixedImageDimension));


    // Extrernal matrices are assumed to be RAS to RAS, so we must convert to LPS to LPS
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
