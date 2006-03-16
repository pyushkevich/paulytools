#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkVoxBoCUBImageIOFactory.h"
#include <itkCropImageFilter.h>

#include <iostream>

using namespace std;
using namespace itk;

int usage()
{
  cout << "USAGE: cropimg [flags] input_image output_image" << endl;
  cout << "FLAGS: " << endl;
  cout << "  -i N N N              Index of the crop region" << endl;
  cout << "  -s N N N              Size of the crop region" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  // Enable support for VoxBo
  itk::ObjectFactoryBase::RegisterFactory(itk::VoxBoCUBImageIOFactory::New());

  // Parameters
  unsigned long sz[] = {0, 0, 0}, idx[] = {0, 0, 0};
  
  // Parse parameters
  if(argc < 3) return usage();
  for(int i = 1; i < argc-2; i++)
    {
    if(!strcmp(argv[i], "-s"))
      {
      sz[0] = atoi(argv[++i]);
      sz[1] = atoi(argv[++i]);
      sz[2] = atoi(argv[++i]);
      }
    if(!strcmp(argv[i], "-i"))
      {
      idx[0] = atoi(argv[++i]);
      idx[1] = atoi(argv[++i]);
      idx[2] = atoi(argv[++i]);
      }
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
