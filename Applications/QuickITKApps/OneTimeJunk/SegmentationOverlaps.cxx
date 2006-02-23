#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageRegionConstIterator.h>
#include "itkVoxBoCUBImageIOFactory.h"

#include <iostream>
#include <vector>
#include <string>
#include <set>

using namespace std;
using namespace itk;

int usage()
{
  cout << "overlaps - compute pairwise overlaps between a set of segmentations" << endl;
  cout << "usage: " << endl;
  cout << "  overlaps FILES " << endl;
  cout << "options: " << endl;
  cout << "  FILES    A list of images for which overlaps will be computed." << endl;
  cout << "           In each image, pixels with each non-zero label are treated" << endl;
  cout << "           as a single structure. Overlaps are computed label-wise." << endl;
  return -1;
}

// Image type information
typedef Image<short, 3> ImageType;
typedef ImageType::Pointer ImagePointer;
typedef ImageFileReader<ImageType> ReaderType;
typedef ImageRegionConstIterator<ImageType> IteratorType;
typedef set<short> LabelSet;
typedef LabelSet::const_iterator LabelIterator;

/**
 * Structure representing an image
 */
struct ImageEntry
{
  ImagePointer image;
  set<short> labels;
  string file;
};

int main(int argc, char *argv[])
{
  // Enable support for VoxBo
  itk::ObjectFactoryBase::RegisterFactory(itk::VoxBoCUBImageIOFactory::New());
  
  // Must have at least two images in input
  if(argc < 3) return usage();

  // Read the image files into memory
  vector<ImageEntry> data;
  for(size_t i = 1; i < argc; i++)
    {
      // Read the image
      ReaderType::Pointer fltReader = ReaderType::New();
      fltReader->SetFileName(argv[i]);
    
      // Catch read exceptions
      try { fltReader->Update(); }
      catch(...)
        {
        cerr << "Can't read image " << argv[i] << endl;
        return -1;
        }

      // Create an entry
      ImageEntry entry;
      entry.image = fltReader->GetOutput();
      entry.file = argv[i];
      
      // Compute the list of labels in the image
      IteratorType it(entry.image, entry.image->GetBufferedRegion());
      for( ; !it.IsAtEnd(); ++it)
        if(it.Get() != 0)
          entry.labels.insert(it.Get());

      // Store the entry
      data.push_back(entry);

      cout << "Image : " << entry.file << " has labels ";
      for(LabelIterator it2 = entry.labels.begin(); it2 != entry.labels.end(); ++it2)
        cout  << *it2 << " ";
      cout << endl;
    }

  // Iterate over all pairwise combinations of images
  for(size_t j = 0; j < data.size() - 1; j++)
    for(size_t k = j + 1; k < data.size(); k++)
      {
      // Get the images
      ImagePointer i1 = data[j].image;
      ImagePointer i2 = data[k].image;
      
      // Check that the images have the same size
      if(i1->GetBufferedRegion() != i2->GetBufferedRegion())
        {
        cerr << "Images " << data[j].file << " and " << data[k].file 
          << " have different sizes!";
        return -1;
        }

      // Check that the set of labels is common between the images
      if(data[j].labels != data[k].labels)
        {
        cerr << "Images " << data[j].file << " and " << data[k].file 
          << " have different sets of labels!";
        return -1;
        }
          
      // Iterate over the images for each label
      LabelSet &lset = data[j].labels;
      for(LabelIterator ilab = lset.begin(); ilab != lset.end(); ilab++)
        {
        short label = *ilab;
        unsigned long n1 = 0, n2 = 0, nOverlap = 0;
        IteratorType it1(i1, i1->GetBufferedRegion());
        IteratorType it2(i2, i2->GetBufferedRegion());
        while(!it1.IsAtEnd())
          {
          bool m1 = (label == it1.Get()), m2 = (label == it2.Get());
          if(m1) n1++;
          if(m2) n2++;
          if(m1 && m2) nOverlap++;
          ++it1; ++it2;
          }
        
        // Print out results as a comma separated file
        cout << label << ", ";
        cout << data[j].file << ", ";
        cout << data[k].file << ", ";
        cout << n1 << ", " << n2 << ", " << nOverlap << endl;
        }
      }
}
  
