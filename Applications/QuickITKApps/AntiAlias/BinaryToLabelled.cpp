#include <iostream>
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkAnalyzeImageIO.h"

using namespace itk;
using namespace std;

template <class TImageType> 
void ReadImage(SmartPointer<TImageType> &target, const char *file)
{
  typename ImageFileReader<TImageType>::Pointer reader = 
    ImageFileReader<TImageType>::New();
  reader->SetFileName(file);
  reader->Update();
  target = reader->GetOutput();
}

template <class TImageType> 
void WriteImage(SmartPointer<TImageType> image, const char *file)
{
  typename ImageFileWriter<TImageType>::Pointer writer = 
    ImageFileWriter<TImageType>::New();
  writer->SetFileName(file);
  writer->SetInput(image);
  writer->Update();
}

int main(int argc, char *argv[])
{
  cout << "usage: lcombine target.hdr source1.hdr source2.hdr ... sourcen.hdr" << endl;

  typedef itk::Image<unsigned short,3> OutImageType;
  typedef itk::Image<unsigned char,3> InImageType;
  typedef itk::ImageRegionConstIterator<InImageType> InIteratorType;
  typedef itk::ImageRegionIterator<OutImageType> OutIteratorType;

  // Create input array
  unsigned int nImages = argc - 2;

  // Create output image
  OutImageType::Pointer imgOutput = OutImageType::New();
  typedef itk::ImageFileReader<InImageType> ReaderType;

  // Read the input images
  for(unsigned int i=2;i<argc;i++)
    {
    cout << "Reading image " << argv[i] << endl;

    itk::AnalyzeImageIO::Pointer aio = itk::AnalyzeImageIO::New();
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(argv[i]);
    reader->SetImageIO(aio);
    reader->Update();

    InImageType::Pointer img = reader->GetOutput();
    if(imgOutput->GetBufferedRegion().GetNumberOfPixels() == 0)
      {
      imgOutput->SetRegions(img->GetBufferedRegion());
      imgOutput->Allocate();
      imgOutput->FillBuffer(0);
      }

    // Create an output iterator
    OutIteratorType itOut(imgOutput, imgOutput->GetBufferedRegion());
    InIteratorType itIn(img,img->GetBufferedRegion());

    // Place the labels in the image
    while(!itOut.IsAtEnd())
      {
        if(itIn.Get() != 0) itOut.Set(i - 2 + 35);
        ++itIn; ++itOut;
      }
    }
    
  // Write the output image
  WriteImage(imgOutput,argv[1]);
};
