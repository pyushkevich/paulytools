/**
 * FilterSulci - take a binary image of the grey matter and construct a 
 * mask image that has a high positive value at sulcal curves. Sulci are
 * ridges of curvature: the smaller principal curvature attains a negative 
 * minimum.   
 */ 

#include "ReadWriteImage.h"
#include <iostream>

#include "itkDiscreteGaussianImageFilter.h"
#include "itkLaplacianImageFilter.h"

using namespace std;

int usage()
{
  cout << "usage sulcalfilter [options] input.img output.img\n" << endl;
  cout << "options: " << endl;
  cout << "   -s X.XX         Standard deviation of the smoothing filter" << endl;
  return -1;  
}

int main(int argc, char *argv[])
{
  // Usage
  if(argc < 3) return usage();
  
  // Process the options
  const char *fnInput = argv[-2];
  const char *fnOutput = argv[-1];
  double sigma = 1.0;

  for(unsigned int iArg=1;iArg<argc-2;iArg++)
    {
    if(!strcmp(argv[iArg],"-s"))
      sigma = atof(argv[++iArg]);
    else
      return usage();
    }

  // Define the float based image
  typedef itk::Image<float,3> ImageType;
  ImageType::Pointer imgInput;
  
  // Read the image
  ReadImage(imgInput,fnInput);

  // Compute the laplacian of the image
  typedef itk::DiscreteGaussianImageFilter<ImageType,ImageType> GaussianFilter;
  typedef itk::LaplacianImageFilter<ImageType,ImageType> LaplacianFilter;

  GaussianFilter::Pointer fltGaussian = GaussianFilter::New();
  fltGaussian->SetInput(imgInput);
  fltGaussian->SetVariance(sigma * sigma);

  LaplacianFilter::Pointer fltLaplacian = LaplacianFilter::New();
  fltLaplacian->SetInput(fltGaussian->GetOutput());
  fltLaplacian->Update();
  
  ImageType::Pointer imgOutput = fltLaplacian->GetOutput();

  WriteImage(imgOutput, fnOutput);
}
