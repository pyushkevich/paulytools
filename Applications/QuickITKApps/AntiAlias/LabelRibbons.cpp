#include <iostream>
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkAnalyzeImageIO.h"
#include "itkNeighborhoodIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

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

/**
 * This program will relabel the grey scale with ribbons image, assigning to each 
 * ribbon pixel the label that is closest to it
 **/
int main(int argc, char *argv[])
{
  cout << "usage: labelribbons grey_with_ribbon.hdr segmentation.hdr" << endl;

  // Typedefs
  typedef itk::Image<short,3> ImageType;
  typedef ImageType::Pointer ImagePointer;

  // Read the two images
  ImagePointer imgRibbon, imgLabel;
  ReadImage(imgRibbon,argv[1]);
  ReadImage(imgLabel,argv[2]);

  // Create iterators
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  Iterator it(imgLabel, imgLabel->GetBufferedRegion());

  // Iterate through the image
  it.GoToBegin();
  while(!it.IsAtEnd())
    {
    // Are we on a ribbon now?
    if(it.Get() > 250) 
      {
      // Are we on the grey matter
      if(imgRibbon->GetPixel(it.GetIndex()) > 0)
        {
        // Create a histogram for different labels
        vector<unsigned int> histogram(1024, 0);

        // Keep track of the largest member
        unsigned int iLargest = 255;
        itk::Index<3> idx = it.GetIndex();

        //cout << idx << " : " << flush;

        for(int x=-3;x<=3;x++)
          for(int y=-3;y<=3;y++)
            for(int z=-3;z<=3;z++)
              {
              itk::Index<3> i2 = idx;
              i2[0] += x;i2[1] += y;i2[2] += z;

              if(imgRibbon->GetPixel(i2) > 0) 
                {
                short pix = imgLabel->GetPixel(i2);

                histogram[pix]++;
                if(histogram[pix] > histogram[iLargest] && pix < 250)
                  iLargest = pix;

                //cout << pix << " " << flush;
                }
              }

        // Record the located value
        if(iLargest < 250)
          it.Set(iLargest);

        // cout << " : " << iLargest << endl;
        }
      else
        {
        it.Set(0);
        }
      }

    ++it;
    }

  // 

  // Write the image
  WriteImage(imgLabel, "blended.hdr");
}






  
