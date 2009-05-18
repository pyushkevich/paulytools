#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
  cout << "usage: coordinate_images input.nii outx.nii outy.nii outz.nii";

  typedef itk::Image<double,3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  // Read the file
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[1]);
  reader->Update();
  ImageType::Pointer image = reader->GetOutput();

  // Create coordinate images
  ImageType::Pointer xyz[3];
  for(size_t i = 0; i < 3; i++)
    {
    xyz[i] = ImageType::New();
    xyz[i]->SetRegions(image->GetBufferedRegion());
    xyz[i]->SetOrigin(image->GetOrigin());
    xyz[i]->SetSpacing(image->GetSpacing());
    xyz[i]->SetDirection(image->GetDirection());
    xyz[i]->Allocate();

    for(itk::ImageRegionIteratorWithIndex<ImageType> it(
      xyz[i], xyz[i]->GetBufferedRegion()); !it.IsAtEnd(); ++it)
      {
      it.Set(it.GetIndex()[i]);
      }

    WriterType::Pointer writer = WriterType::New();
    writer->SetInput(xyz[i]);
    writer->SetFileName(argv[2+i]);
    writer->Update();
    }
}
