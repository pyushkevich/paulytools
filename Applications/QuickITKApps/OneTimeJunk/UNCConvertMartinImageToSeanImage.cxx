#include <itkImage.h>
#include <itkImageRegionIteratorWithIndex.h>
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

  // Create an output image region
  ImageType::Pointer imgMartin = fltReader->GetOutput();
  ImageType::RegionType rgnInput = imgMartin->GetBufferedRegion();
  ImageType::RegionType rgnOutput = rgnInput;
  rgnOutput.SetSize(1, rgnInput.GetSize(2));
  rgnOutput.SetSize(2, rgnInput.GetSize(1));

  // Allocate the output image
  ImageType::Pointer imgSean = ImageType::New();
  imgSean->SetRegions(rgnOutput);
  imgSean->Allocate();
  imgSean->FillBuffer(0);
  
  ImageType::SpacingType spcInput = imgMartin->GetSpacing();
  ImageType::SpacingType spcOutput = spcInput;
  spcOutput[1] = spcInput[2];
  spcOutput[2] = spcInput[1];
  imgSean->SetSpacing(spcOutput);

  // Copy the image contents
  typedef ImageRegionIteratorWithIndex<ImageType> IteratorType;
  IteratorType it(imgMartin, rgnInput);
  while(!it.IsAtEnd())
    {
    ImageType::IndexType idx1 = it.GetIndex();
    ImageType::IndexType idx2 = idx1;
    idx2.SetElement(1, rgnInput.GetSize(2) - idx1.GetElement(2) - 1);
    idx2.SetElement(2, idx1.GetElement(1));
    if(rgnOutput.IsInside(idx2))
      imgSean->SetPixel(idx2, it.Get());
    ++it;
    }

  // Save the output image
  typedef ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer fltWriter = WriterType::New();
  fltWriter->SetInput(imgSean);
  fltWriter->SetFileName(argv[2]);
  fltWriter->Update();
}
