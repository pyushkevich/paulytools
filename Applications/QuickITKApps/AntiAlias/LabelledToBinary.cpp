#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImage.h"

#include <iostream>
#include <sstream>
#include <cstring>
#include <string>
#include <vector>

using namespace itk;
using namespace std;

int main(int argc, char **argv)
{
  // This program takes an image whose intensities are labels 
  // and decomposes it into binary images for each label

  if(argc!=3) 
    {
    cout << "Usage lsplit file.img outrootname" << endl;
    return -1;
    }

  // Load the input image
  typedef itk::Image<short,3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  
  ReaderType::Pointer fltReader = ReaderType::New();
  ImageType::Pointer imgInput = ImageType::New();
  
  fltReader->SetFileName(argv[1]);
  imgInput = fltReader->GetOutput();
  fltReader->Update();

  // Allocate a histogram of the image
  vector<unsigned int> xHistogram(0x10000,0);

  // Create an iterator
  typedef itk::ImageRegionConstIterator<ImageType> IteratorType;
  IteratorType itInput(imgInput,imgInput->GetBufferedRegion());
  
  // Use the iterator to construct a histogram of the image
  while(!itInput.IsAtEnd())
    {
    // Update the count
    short val = itInput.Get();
    xHistogram[val]++;
    ++itInput;
    }

  // For each histogram entry, construct an image
  for(unsigned int i=0;i<xHistogram.size();i++)
    {
    if(xHistogram[i] > 0) 
      {
      cout << "Pixel " << i << " Count " << xHistogram[i] << endl;
      
      typedef itk::Image<unsigned char,3> BinImageType;
      typedef itk::ImageRegionIterator<BinImageType> BinIteratorType;

      BinImageType::Pointer imgBin = BinImageType::New();
      imgBin->SetRegions(imgInput->GetBufferedRegion());
      imgBin->Allocate();
      imgBin->FillBuffer(0);
      
      IteratorType itSource(imgInput,imgInput->GetBufferedRegion());
      BinIteratorType itTarget(imgBin,imgBin->GetBufferedRegion());
      
      while(!itSource.IsAtEnd())
        {
        if(itSource.Get() == i)
          itTarget.Set(127);
        ++itSource;++itTarget;
        }

      ostringstream fn;
      fn << argv[2] << "." << (i / 100) << ((i % 100) / 10) << (i % 10) << ".img";
      cout << "Writing " << fn.str().c_str() << endl;

      typedef itk::ImageFileWriter<BinImageType> WriterType;
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName(fn.str().c_str());
      writer->SetInput(imgBin);
      writer->Update();
      }
    }
}
