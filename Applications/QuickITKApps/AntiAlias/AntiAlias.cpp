#include "itkAntiAliasBinaryImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkBinaryThresholdImageFilter.h"
        
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;
        
class HistoEntry {
public:
  unsigned int count;
  short val;
  int lb[3], ub[3]; 
};


int main(int argc, char **argv)
{
  if(argc != 2)
    {
    cout << "USAGE: aalias inputimage outputimage" << endl;
    cout << "No more help" << endl;
    return -1;
    }

  // Load the input image
  typedef itk::Image<short,3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  
  ReaderType::Pointer fltReader = ReaderType::New();
  ImageType::Pointer imgInput = ImageType::New();
  
  fltReader->SetFileName(argv[1]);
  fltReader->Update();

  // Create the default histogram entry
  HistoEntry dummy;
  dummy.count = 0;
  for(unsigned int i=0;i<3;i++)
    {
    dummy.lb[i] = (int) imgInput->GetBufferedRegion().GetSize(i);
    dummy.ub[i] = 0;
    }

  // Allocate a histogram of the image
  vector<HistoEntry> xHistogram(0x10000,dummy);

  // Create an iterator
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;
  IteratorType itInput(imgInput,imgInput->GetBufferedRegion());
  
  // Use the iterator to construct a histogram of the image
  while(!itInput.IsAtEnd())
    {
    // Update the count
    short val = itInput.Get();
    xHistogram[val].count++;
    xHistogram[val].val = val;

    // Update the bounding box
    for(unsigned int i=0;i<3;i++)
      {
      int idx = itInput.GetIndex().GetIndex()[i];
      if(idx < xHistogram[val].lb[i])
        xHistogram[val].lb[i] = idx;
      if(idx > xHistogram[val].ub[i])
        xHistogram[val].ub[i] = idx;
      }
    }

  // For each non-zero type, apply the antialias application
  for(vector<HistoEntry>::iterator it=xHistogram.begin();it!=xHistogram.end();it++)
    { 
    if(it->count > 0)
      {
      
      // Create a region corresponding to the bounding box
      itk::ImageRegion<3> roi;
      for(unsigned int i=0;i<3;i++)
        {
        roi.SetIndex(i,it->lb[i]);
        roi.SetSize(i,1 + it->ub[i] - it->lb[i]);
        }

      // Pad the region sufficiently and crop
      roi.PadByRadius(10);
      roi.Crop(imgInput->GetBufferedRegion());
      
      // Extract a binary image for the region of interest of the image
      typedef itk::Image<unsigned char, 3> ByteImageType;
      typedef itk::BinaryThresholdImageFilter<ImageType,ByteImageType> BinarizerType;
      BinarizerType::Pointer fltBinarizer = BinarizerType::New();
      fltBinarizer->SetInput(imgInput);
      fltBinarizer->SetLowerThreshold(dummy.val);
      fltBinarizer->SetUpperThreshold(dummy.val);
      fltBinarizer->SetInsideValue(255);
      fltBinarizer->SetOutsideValue(0);
      
      ByteImageType::Pointer imgBinary = fltBinarizer->GetOutput();
      imgBinary->SetRequestedRegion(roi);
      imgBinary->Update();

      // Apply antialiasing to the binary image
      typedef itk::AntiAliasBinaryImageFilter<ByteImageType,ByteImageType>
        AntiFilterType;
      AntiFilterType::Pointer fltAnti = AntiFilterType::New();
      fltAnti->SetInput(imgBinary);
      fltAnti->SetMaximumRMSError(0.07);
      fltAnti->Update();

      // Save the output
      std::ostringstream oss;
      oss << argv[2] << "." << it->val << ".gipl";
      
      typedef itk::ImageFileWriter<ByteImageType> WriterType;
      WriterType::Pointer fltWriter = WriterType::New();
      fltWriter->SetInput(fltAnti->GetOutput());
      fltWriter->SetFileName(oss.str().c_str());
      fltWriter->Update();
      }
    }
}
