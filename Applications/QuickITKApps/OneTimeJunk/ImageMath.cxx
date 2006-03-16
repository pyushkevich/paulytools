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
  cout << "USAGE: imagemath [command] output_image" << endl;
  cout << "COMMANDS: " << endl;
  cout << "  -add   img1 img2 " << endl;
  cout << "  -mult  img1 img2 " << endl;
  cout << "  -scale img const " << endl;
  return -1;
}

typedef Image<float, 3> ImageType;

:q
typedef ImageFileReader<ImageType> ReaderType;

int main(int argc, char *argv[])
{
  // Enable support for VoxBo
  itk::ObjectFactoryBase::RegisterFactory(itk::VoxBoCUBImageIOFactory::New());

  // Parse parameters
  if(argc < 3) return usage();

  // Apply commands to the image
  for(int i = 1; i < argc-2; i++)
    {
    if(!strcmp(argv[i], "-add"))
      {
      }
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
