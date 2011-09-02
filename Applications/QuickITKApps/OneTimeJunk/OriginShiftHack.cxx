#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkVoxBoCUBImageIOFactory.h"
#include <itkImageRegionConstIteratorWithIndex.h>

#include <iostream>
#include <vector>

using namespace std;
using namespace itk;

int usage()
{
  cout << "USAGE: maskmunge input_images" << endl;
  cout << "DESCRIPTION: " << endl;
  cout << "  This program takes a list of images. It computes the" << endl;
  cout << "  smallest rect. region that contains all non-zero pixels" << endl;
  cout << "  in all images. " << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  // Enable support for VoxBo
  itk::ObjectFactoryBase::RegisterFactory(itk::VoxBoCUBImageIOFactory::New());

  // Parse parameters
  if(argc < 3) return usage();

  // Super-region
  ImageRegion<3> xRegion;
  bool regionEmpty = true;

  // Input image array
  typedef Image<float, 3> ImageType;
  vector<ImageType::Pointer> vInput;
  for(int i = 1; i < argc; i++)
    {
    // Read the image
    typedef ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(argv[i]);

    // Report read errors
    try { reader->Update(); }
    catch(ExceptionObject &exc)
      { cerr << exc; return usage(); }

    // Store the image
    ImageType::Pointer img = reader->GetOutput();

    // Find the extent of non-zero pixels in the image
    typedef ImageRegionConstIteratorWithIndex<ImageType> Iter;
    Iter it(img,img->GetBufferedRegion());

    for( ; !it.IsAtEnd(); ++it)
      {
      if(it.Value() != 0)
        {
        if(regionEmpty) 
          {
          for(int j = 0; j < 3; j++)
            { xRegion.SetSize(j, 1); xRegion.SetIndex(j, it.GetIndex()[j]); }
          regionEmpty = false;
          }
        else
          {
          for(int j = 0; j < 3; j++)
            {
            if(xRegion.GetIndex(j) > it.GetIndex()[j])
              {
              xRegion.SetSize(j, xRegion.GetSize(j) + it.GetIndex()[j] - xRegion.GetIndex(j));
              xRegion.SetIndex(j, it.GetIndex()[j]);
              }
            if(xRegion.GetIndex(j) + xRegion.GetSize(j) <= it.GetIndex()[j])
              xRegion.SetSize(j, 1 + it.GetIndex()[j] - xRegion.GetIndex(j));
            }
          }
        }
      }
    cout << "Image " << i << "; Region is " << xRegion << endl;
    }

  // Print out the region
  cout << "Total Mask Region: " << xRegion << endl;
}
