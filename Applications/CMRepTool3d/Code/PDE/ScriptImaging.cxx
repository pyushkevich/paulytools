#include "ScriptInterface.h"
#include "ITKImageWrapper.h"

#include <itkCastImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkUnaryFunctorImageFilter.h>
#include <itkDerivativeImageFilter.h>

namespace medialpde {

FloatImage::FloatImage()
{
  xImage = ITKImageWrapperFactory<float>::NewImageWrapper();
  xGradient[0] = ITKImageWrapperFactory<float>::NewImageWrapper();
  xGradient[1] = ITKImageWrapperFactory<float>::NewImageWrapper();
  xGradient[2] = ITKImageWrapperFactory<float>::NewImageWrapper();
}
  
FloatImage::~FloatImage()
{
  delete xImage;
  delete xGradient[0];
  delete xGradient[1];
  delete xGradient[2];
}

void FloatImage::LoadFromFile(const char *file)
{
  // Load an image from ITK
  xImage->LoadFromFile(file);

  // Clear the gradient images
  for(unsigned int i = 0; i < 3; i++)
    xGradient[i]->SetInternalImage(NULL);
}

void FloatImage::LoadGradientFromFile(unsigned int iComponent, const char *file)
{
  // Load an image from ITK
  xGradient[iComponent]->LoadFromFile(file);
}

void FloatImage::SaveToFile(const char *file)
{
  // Save the ITK image
  xImage->SaveToFile(file);
}

void FloatImage::SaveGradientToFile(unsigned int iComponent, const char *file)
{
  // Save the ITK image
  xGradient[iComponent]->SaveToFile(file);
}

void FloatImage::SetToGradientMagnitude(BinaryImage *imgSource, double xSigma)
{
  // Blur the image with a Gaussian
  typedef itk::Image<unsigned char, 3> BinaryImageType;
  typedef itk::Image<float, 3> FloatImageType;
  typedef itk::CastImageFilter<BinaryImageType, FloatImageType> CasterType;
  typedef itk::DiscreteGaussianImageFilter<FloatImageType, FloatImageType>
    GaussianType;
  typedef itk::GradientMagnitudeImageFilter<FloatImageType, FloatImageType>
    GradMagType;

  // Create a pipeline of filters
  CasterType::Pointer fltCaster = CasterType::New();
  fltCaster->SetInput(imgSource->xImage->GetInternalImage());

  GaussianType::Pointer fltGaussian = GaussianType::New();
  fltGaussian->SetInput(fltCaster->GetOutput());
  fltGaussian->SetVariance(xSigma * xSigma);

  GradMagType::Pointer fltGrad = GradMagType::New();
  fltGrad->SetInput(fltGaussian->GetOutput());
  fltGrad->Update();

  // Set the internal image
  xImage->SetInternalImage(fltGrad->GetOutput());
}

bool FloatImage::IsGradientAvailable()
{
  return xGradient[0]->IsImageLoaded();
}

void FloatImage::SetOutsideValue(float xOutsideValue)
{
  this->xOutsideValue = xOutsideValue;
}

class BinaryFloatFunctor 
{
public:
  float operator() (unsigned char x)
    { return (x == 0) ? -1.0f : 1.0f; } 
};

void FloatImage
::SetToBlurredBinary(BinaryImage *imgSource, double xSigma)
{
  // Blur the image with a Gaussian
  typedef itk::Image<unsigned char, 3> BinaryImageType;
  typedef itk::Image<float, 3> FloatImageType;
  typedef itk::UnaryFunctorImageFilter<
    BinaryImageType, FloatImageType, BinaryFloatFunctor> MapperType;
  typedef itk::DiscreteGaussianImageFilter<FloatImageType, FloatImageType>
    GaussianType;
  typedef itk::DerivativeImageFilter<FloatImageType, FloatImageType>
    DerivativeType;

  // Create a pipeline of filters to blur the image
  MapperType::Pointer fltMapper = MapperType::New();
  fltMapper->SetInput(imgSource->xImage->GetInternalImage());

  GaussianType::Pointer fltGaussian = GaussianType::New();
  fltGaussian->SetInput(fltMapper->GetOutput());
  fltGaussian->SetVariance(xSigma * xSigma);
  fltGaussian->Update();

  // Set the internal image
  xImage->SetInternalImage(fltGaussian->GetOutput());
  
  // Now compute the gradient magnitude by convolution with a Gaussian directional
  // derivative filter in each cardinal direction
  for(unsigned int d = 0; d < 3; d++)
    {
    DerivativeType::Pointer fltDerivative = DerivativeType::New();
    fltDerivative->SetInput(fltGaussian->GetOutput());
    fltDerivative->SetOrder( 1 );
    fltDerivative->SetDirection( d );
    fltDerivative->Update();
    xGradient[d]->SetInternalImage(fltDerivative->GetOutput());
    }
}

BinaryImage::BinaryImage()
{
  xImage = ITKImageWrapperFactory<unsigned char>::NewImageWrapper();
}
  
BinaryImage::~BinaryImage()
{
  delete xImage;
}

void BinaryImage::LoadFromFile(const char *file)
{
  // Load an image from ITK
  xImage->LoadFromFile(file);
}

} // Namespace
