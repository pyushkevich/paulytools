/**
 * Program: MaskImage
 * Author: Paul Yushkevich
 * Purpose: Convert a binary image to a VTK mesh
 */

#include "itkMaskImageFilter.h"
#include "ReadWriteImage.h"

using namespace std;

int usage()
{
  cout << "maskimage: applies a binary mask to an image" << endl;
  cout << "usage: " << endl;
  cout << "   maskimage [options] input.img mask.img output.img" << endl;
  cout << "options: " << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  // Command line options
  double aaParm = 0.024;
  double xScale = 1.0;

  // Check parameters  
  if(argc < 3) return usage();
  
  // Read filename parameters
  const char *fnSource = argv[argc - 3];
  const char *fnMask   = argv[argc - 2];
  const char *fnOutput = argv[argc - 1];

  // Define images
  typedef itk::Image<unsigned short, 3> ImageType;
  ImageType::Pointer imgSource, imgMask, imgOutput;
  
  // Read images
  ReadImage(imgSource, fnSource);
  ReadImage(imgMask, fnMask);

  // Apply filter
  typedef itk::MaskImageFilter<ImageType, ImageType, ImageType> FilterType;
  FilterType::Pointer fltMask = FilterType::New();
  fltMask->SetInput(0, imgSource);
  fltMask->SetInput(1, imgMask);
  fltMask->Update();

  // Write output
  imgOutput = fltMask->GetOutput();
  WriteImage(imgOutput, fnOutput);
}
  
