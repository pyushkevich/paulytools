#include "imagesurface.h"

ImageSurface
::ImageSurface(ImageType *image)
{
  m_Image = image;
  
  m_Interpolator = InterpolatorType::New();
  m_Interpolator->SetInputImage(image);
}

double
ImageSurface
::getFunction(double x,double y,double z) 
{
  // Make an 'index' for itk
  InterpolatorType::ContinuousIndexType index;
  index[0] = x; index[1] = y; index[2] = z;
  
  // Interpolate the image linearly to get the value of the function
  return m_Interpolator->EvaluateAtContinuousIndex(index);
}
  
