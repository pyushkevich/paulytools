#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <iostream>

using namespace std;
using namespace itk;

int main(int argc, char *argv[])
{
  cout << "USAGE: program input_image output_image" << endl;

  // Read the input image
  typedef Image<short, 3> ImageType;
  typedef ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(argv[1]);
  fltReader->Update();

  // Save the output image
  typedef ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer fltWriter = WriterType::New();
  fltWriter->SetInput(fltReader->GetOutput());
  fltWriter->SetFileName(argv[2]);
  fltWriter->Update();
}
