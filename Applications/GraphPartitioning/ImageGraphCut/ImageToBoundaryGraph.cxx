/**
 * Convert image to a graph of boundary pixels. Requires selecting pixel corners 
 * that are not surrounded by eight pixels of the same intensity
 */
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkNeighborhoodIterator.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <iostream>
#include <fstream>
#include <vector>


using namespace std;
using namespace itk;

int usage()
{
  cerr << "Usage: gcut_makepts input.img out.txt" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  // Get command line parameters
  if(argc < 3)
    return usage();

  char *input = argv[argc-2];
  char *output = argv[argc-1];

  cerr << "Processing " << input << " and " << output << endl;
  
  if(!input || !output)
    return usage();

  typedef Image<unsigned char, 3> ImageType;
  typedef Image<short, 3> ShortImageType;
  typedef ImageFileReader<ImageType> ReaderType;

  // Read the image
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(input);
  fltReader->Update();
  ImageType::Pointer imgBinary = fltReader->GetOutput();

  // Count the number of voxels - in order to report the potential graph size
  unsigned int nVertices = 0;
  typedef ImageRegionConstIterator<ImageType> SillyIterator;
  SillyIterator itCount(imgBinary,imgBinary->GetBufferedRegion());
  while(!itCount.IsAtEnd())
    {
    if(itCount.Value() != 0) nVertices++;
    ++itCount;
    }

  // Report the size
  cout << "There are " << nVertices << " non-zero voxels in the image " << endl;

  // Compute the connected components
  typedef ConnectedComponentImageFilter<ImageType,ShortImageType> CCFilter;
  CCFilter::Pointer fltConComp = CCFilter::New();
  fltConComp->SetInput(imgBinary);
  fltConComp->SetFullyConnected(true);
  fltConComp->Update();

  // Relabel the components
  typedef RelabelComponentImageFilter<ShortImageType,ShortImageType> RCFilter;
  RCFilter::Pointer fltRelabel = RCFilter::New();
  fltRelabel->SetInput(fltConComp->GetOutput());
  fltRelabel->Update();

  // Print the statistics about the connected components
  cout << "There are " << fltRelabel->GetNumberOfObjects() << " connected components." << endl;
  cout << "Largest component has " << fltRelabel->GetSizeOfObjectInPixels(1) << " pixels." << endl;

  // Define the region of iteration
  ImageType::RegionType rgnIter = imgBinary->GetBufferedRegion();
  rgnIter.SetSize(0, rgnIter.GetSize(0) - 1);
  rgnIter.SetSize(1, rgnIter.GetSize(1) - 1);
  rgnIter.SetSize(2, rgnIter.GetSize(2) - 1);
  
  // Create an iterator to traverse the image
  typedef itk::NeighborhoodIterator<ShortImageType> IterType;
  itk::Size<3> szRadius = {{1,1,1}};
  IterType it(szRadius, fltRelabel->GetOutput(), rgnIter);

  // Get the offsets that are tested
  unsigned int sx = it.GetStride(0);
  unsigned int sy = it.GetStride(1);
  unsigned int sz = it.GetStride(2);
  unsigned int c = sx + sy + sz;
  unsigned int off[] = 
    {c, c+sx, c+sy, c+sz, c+sx+sy, c+sx+sz, c+sy+sz, c+sx+sy+sz};

  typedef Point<double,3> PointType;
  vector<PointType> vPoints;

  // Iterate over the image
  while(!it.IsAtEnd())
    {
    // Get the index
    ImageType::IndexType idx = it.GetIndex();

    // Compute the sum of the eight adjacent pixels
    short sum = 0;
    for(unsigned int i=0;i<8;i++)
      sum += it.GetPixel(off[i]) == 1 ? 1 : 0;

    // If the sum is not 0 or 8 we have a border pixel
    if(sum != 0 && sum != 8)
      {
      PointType pt;
      imgBinary->TransformIndexToPhysicalPoint(idx,pt);
      
      vPoints.push_back(pt);
      }

    ++it;
    }

  cout << "The contour consists of " << vPoints.size() << " vertices." << endl;

  ofstream fout(output, ios_base::out);
  fout << "3" << endl;
  fout << vPoints.size() << endl;
  for(unsigned int j=0;j<vPoints.size();j++)
    fout << vPoints[j][0] << " " << vPoints[j][1] << " " << vPoints[j][2] << endl;
}
