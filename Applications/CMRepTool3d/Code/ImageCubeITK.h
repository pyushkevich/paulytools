#ifndef __ImageCubeITK_h_
#define __ImageCubeITK_h_

#include "imaging.h"
#include "itkImage.h"

template <typename TPixel> 
class ImageCubeITK : public ImageCube<TPixel>
{
public:
  typedef itk::Image<TPixel, 3> ImageType;

  // Set the background
  ImageCubeITK()
    { xBackgroundValue = 0.0; }

  // Copy the contents of an ITK image into the cube
  void SetImage(ImageType *imgSource, const TPixel &xBackground);

  // Handle the out-of-bounds situation
  virtual float handleOutsideVoxel(float x, float y, float z)
    { return xBackgroundValue; }

private:
  float xBackgroundValue;
};

#endif // __ImageCubeITK_h_
