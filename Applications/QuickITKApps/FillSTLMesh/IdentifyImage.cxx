#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkVoxBoCUBImageIOFactory.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
  // Enable support for VoxBo
  itk::ObjectFactoryBase::RegisterFactory(itk::VoxBoCUBImageIOFactory::New());
  
  if(argc < 2)
    {
    cerr << "Usage: identimg image_file" << endl;
    return -1;
    }

  typedef itk::Image<float, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(argv[1]);
  fltReader->Update();

  ImageType::Pointer image = fltReader->GetOutput();

  cout << "Image: " << argv[1] << endl;
  cout << "Dimensions: " << image->GetBufferedRegion() << endl;
  cout << "Spacing: " << image->GetSpacing() << endl;
  cout << "Origin: " << image->GetOrigin() << endl;

  // Compute the image histogram
  double iMin, iMax, iSum = 0;
  bool flagFirst = true;
  itk::ImageRegionIterator<ImageType> it(image, image->GetBufferedRegion());
  for( ; !it.IsAtEnd(); ++it)
    {
    double val = it.Value();
    if(!vnl_math_isfinite(val))
      continue;
    
    if(flagFirst)
      { 
      iMin = iMax = iSum = val;
      flagFirst = false;
      }
    else
      {
      iMax = val > iMax ? val : iMax;
      iMin = val < iMin ? val : iMin;
      iSum += val;
      }
    }

  cout << "Intensity Range: [" << iMin << ", " << iMax << "]" << endl;
  cout << "Intensity Mean: " << 
    iSum / image->GetBufferedRegion().GetNumberOfPixels() << endl;

  return 0;
}
