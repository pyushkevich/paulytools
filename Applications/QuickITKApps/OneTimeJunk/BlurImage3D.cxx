#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkVoxBoCUBImageIOFactory.h"
#include <itkDiscreteGaussianImageFilter.h>

#include <iostream>

using namespace std;
using namespace itk;

int usage()
{
  cout << "USAGE: blurimg [flags] input_image output_image" << endl;
  cout << "FLAGS: " << endl;
  cout << "  -s X.X              Sigma of the blurring kernel (mm)" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  // Enable support for VoxBo
  itk::ObjectFactoryBase::RegisterFactory(itk::VoxBoCUBImageIOFactory::New());

  // Parameters
  double xStdDev = 1.0;
  
  // Parse parameters
  if(argc < 3) return usage();
  for(int i = 1; i < argc-2; i++)
    {
    if(!strcmp(argv[i], "-s"))
      xStdDev = atof(argv[++i]);
    else return usage();
    }

  // Read the image
  typedef Image<float, 3> ImageType;
  typedef ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[argc-2]);

  // Report read errors
  try { reader->Update(); }
  catch(ExceptionObject &exc)
    { cerr << exc; return usage(); }

  // Create the gaussian filter
  typedef DiscreteGaussianImageFilter<ImageType, ImageType> Gaussian;
  Gaussian::Pointer gaussian = Gaussian::New();
  gaussian->SetInput(reader->GetOutput());
  gaussian->SetVariance(xStdDev * xStdDev);
  gaussian->Update();

  // Write the image
  typedef ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(gaussian->GetOutput());
  writer->SetFileName(argv[argc-1]);
  writer->Update();
}
