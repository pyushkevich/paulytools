#ifndef __ImageSurface_h_
#define __ImageSurface_h_

#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "blobmodl.h"

/** 
 * class ImageSurface
 *
 * A class derived from ISurface that is used to compute implicit surfaces
 * in floating point images such as the images used in level set methods
 */
class ImageSurface : public ISurface 
{
public:
  typedef itk::Image<float,3> ImageType;
  typedef itk::LinearInterpolateImageFunction<ImageType,double> InterpolatorType;

  /** Constructor that takes an image as input */
  ImageSurface(ImageType *image);

  /** Compute the value of the function at a point in space */
  double getFunction(double x,double y,double z);
  
private:
  /** The image whose level set we will extract */
  ImageType::Pointer m_Image;

  /** The linear interpolator */
  InterpolatorType::Pointer m_Interpolator;
};

#endif
