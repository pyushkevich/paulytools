#include "blobmodl.h"
#include "imagesurface.h"
#include <itkImageFileReader.h>

/**
 * PROGRAM bettermesh.cpp
 * 
 * This program computes a better implicit surface from a floating point image
 * than does a method such as marching cubes. The reason it's better is because
 * the resulting mesh follows the ridges of the structure of interest. These ridges
 * are based on smoothing the image with a small Gaussian kernel
 */

int main(int argc, char *argv[])
{
  // The user will pass in the image and the kernel size
  const char *fnImage = argv[1];  
  double alpha = atof(argv[2]);

  // Read the image 
  typedef itk::Image<float,3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(fnImage);
  fltReader->Update();

  // Compute the implicit surface of the image
  ISurface *srfImage = new ImageSurface(fltReader->GetOutput());
  srfImage->recompute(false);
  cout << "Starting surface has " << srfImage->triangles.getSize() << " triangles " << endl;

  return 0;
}


