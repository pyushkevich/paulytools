#ifndef __ImageCubeITK_h_
#define __ImageCubeITK_h_

#include "imaging.h"
#include "itkImage.h"

template <typename TPixel> 
class ImageCubeITK : public ImageCube<TPixel>
{
public:
  typedef itk::Image<TPixel, 3> ImageType;

  // Copy the contents of an ITK image into the cube
  void SetImage(ImageType *imgSource, const TPixel &xBackground);
};

#endif // __ImageCubeITK_h_
