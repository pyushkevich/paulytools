/**
 * Program: ThresholdImage
 * Author: Paul Yushkevich
 * Purpose: Threshhold an image
 */

#include "itkBinaryThresholdImageFilter.h"
#include "ReadWriteImage.h"

using namespace std;

int usage()
{
  cout << "threshimage: applies a threshold to an image" << endl;
  cout << "usage: " << endl;
  cout << "   threshimage input.img output.img tLower tUpper vInside vOutside" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  // Check parameters  
  if(argc < 6) return usage();
  
  // Read filename parameters
  const char *fnSource = argv[argc - 6];
  const char *fnOutput = argv[argc - 5];

  // Command line options
  short tLower = atoi(argv[argc - 4]);
  short tUpper = atoi(argv[argc - 3]);
  short vInside = atoi(argv[argc - 2]);
  short vOutside = atoi(argv[argc - 1]);

  // Define images
  typedef itk::Image<short, 3> ImageType;
  ImageType::Pointer imgSource, imgOutput;
  
  // Read images
  ReadImage(imgSource, fnSource);

  // Apply filter
  typedef itk::BinaryThresholdImageFilter<ImageType,ImageType> FilterType;
  FilterType::Pointer fltThresh = FilterType::New();
  fltThresh->SetInput(imgSource);
  fltThresh->SetLowerThreshold(tLower);
  fltThresh->SetUpperThreshold(tUpper);
  fltThresh->SetInsideValue(vInside);
  fltThresh->SetOutsideValue(vOutside);

  fltThresh->Update();

  // Write output
  imgOutput = fltThresh->GetOutput();
  WriteImage(imgOutput, fnOutput);
}
  
