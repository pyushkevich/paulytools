/**
 * Program: ThresholdImage
 * Author: Paul Yushkevich
 * Purpose: Threshhold an image
 */

#include "itkBinaryThresholdImageFilter.h"
#include "itkVoxBoCUBImageIOFactory.h"
#include "ReadWriteImage.h"

using namespace std;

int usage()
{
  cout << "threshimage: applies a threshold to an image" << endl;
  cout << "usage: " << endl;
  cout << "   threshimage [options] input.img output.img tLower tUpper vInside vOutside" << endl;
  cout << "options: " << endl;
  cout << "  -f              Floating point input and output (def: short)" << endl;
  return -1;
}

template<class TPixel>
inline TPixel atox(char *data)
{
  istringstream iss(data);
  TPixel x;
  iss >> x;
  return x;
}

template <class TPixel>
inline int thresh(int argc, char *argv[])
{
  
  // Read filename parameters
  const char *fnSource = argv[argc - 6];
  const char *fnOutput = argv[argc - 5];

  // Command line options
  TPixel tLower = atox<TPixel>(argv[argc - 4]);
  TPixel tUpper = atox<TPixel>(argv[argc - 3]);
  TPixel vInside = atox<TPixel>(argv[argc - 2]);
  TPixel vOutside = atox<TPixel>(argv[argc - 1]);

  // Define images
  typedef itk::Image<TPixel, 3> ImageType;
  typename ImageType::Pointer imgSource, imgOutput;
  
  // Read images
  ReadImage(imgSource, fnSource);

  // Apply filter
  typedef itk::BinaryThresholdImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer fltThresh = FilterType::New();
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

int main(int argc, char *argv[])
{
  // Enable support for VoxBo
  itk::ObjectFactoryBase::RegisterFactory(itk::VoxBoCUBImageIOFactory::New());
  
  // Check parameters  
  if(argc < 6) return usage();

  // Read command line parameters
  string type = "short";
  for(size_t i = 1; i < argc-6;i++)
    if(!strcmp(argv[i], "-f"))
      type = "float";

  // Run the program
  if(type == "short")
    return thresh<short>(argc, argv);
  else if(type == "float")
    return thresh<float>(argc, argv);
}
  
