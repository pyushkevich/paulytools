#include "itkImageFileReader.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
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

  return 0;
}
