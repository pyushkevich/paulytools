#ifndef __ITKImageWrapper_h_
#define __ITKImageWrapper_h_

namespace itk {
  template <typename TPixel, unsigned int VDim> class Image;
};

template<typename TPixel>
class ITKImageWrapper
{
public:
  // The wrapped ITK image
  typedef itk::Image<TPixel, 3> ImageType;
  
  // Load the image from a file
  virtual void LoadFromFile(const char *file) = 0;

  // Save the image to a file
  virtual void SaveToFile(const char *file) = 0;

  // Set the internal image to some ITK image, can be NULL
  virtual void SetInternalImage(ImageType *xImage) = 0;

  // Check whether an image is loaded into the wrapper
  virtual bool IsImageLoaded() = 0;

  // Get the basic image properties
  virtual unsigned int GetImageSize(unsigned int d) = 0;
  virtual double GetImageSpacing(unsigned int d) = 0;
  virtual double GetImageOrigin(unsigned int d) = 0;

  // Interpolate the image at a continuous index, or return background if out of bounds
  virtual float Interpolate(float x, float y, float z, float xBackground) = 0;

  // Interpolate the image at a continuous index, or return background if out of bounds
  virtual float InterpolateNearestNeighbor(float x, float y, float z, float xBackground) = 0;

  // Get the internal ITK image pointer
  virtual ImageType *GetInternalImage() = 0;
};

/**
 * A factory class for creating new instances of image wrappers
 */
template<class TPixel>
class ITKImageWrapperFactory
{
public:
  static ITKImageWrapper<TPixel> *NewImageWrapper();
};

#endif // __ITKImageWrapper_h_
