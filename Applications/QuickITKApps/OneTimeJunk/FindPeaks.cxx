#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkVoxBoCUBImageIOFactory.h"
#include "itkConstNeighborhoodIterator.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;
using namespace itk;

int usage()
{
  cout << "USAGE: findpeaks [flags] input_image" << endl;
  cout << "FLAGS: " << endl;
  cout << "  -r N         Radius of the neighborhood for peak extraction" << endl;
  cout << "  -t X.XX      Set the lower threshold for the peaks (def: 0.0)" << endl;
  cout << "  -d file.txt  Write Talairach coordinates to a file that can be" << endl;
  cout << "               used with the Talairach deamon software" << endl;
  cout << "  -i image     Save an image of peak voxels (others set to zero) " << endl;
  return -1;
}

typedef Image<float, 3> ImageType;
struct Peak
{
  ImageType::IndexType index;
  ImageType::PointType point;
};

ImageType::PointType MNI2Tal(ImageType::PointType &p)
{
  ImageType::PointType q;
  if(p[2] >= 0.0) 
    {
    q[0] = 0.9900 * p[0];
    q[1] = 0.9688 * p[1] + 0.046 * p[2];
    q[2] = -0.0485 * p[1] + 0.9189 * p[2];
    }
  else
    {
    q[0] = 0.9900 * p[0];
    q[1] = 0.9688 * p[1] + 0.042 * p[2];
    q[2] = -0.0485 * p[1] + 0.8390 * p[2];
    }
  return q;
}

bool operator < (const Peak &p1, const Peak &p2)
{
  return true;
}

int main(int argc, char *argv[])
{
  // Parameter input
  size_t xRad = 1;
  float xThresh = 0.0f;
  string fnTal = "", fnImage = "";

  if(argc < 2) return usage();
  for(size_t j = 1; j < argc-1; j++)
    {
    if(!strcmp(argv[j],"-r"))
      xRad = atoi(argv[++j]);
    else if(!strcmp(argv[j],"-t"))
      xThresh = atof(argv[++j]);
    else if(!strcmp(argv[j],"-d"))
      fnTal = argv[++j];
    else if(!strcmp(argv[j],"-i"))
      fnImage = argv[++j];
    else return usage();
    }
  
  // Enable support for VoxBo
  itk::ObjectFactoryBase::RegisterFactory(itk::VoxBoCUBImageIOFactory::New());

  // Read the image
  typedef ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[argc-1]);

  // Report read errors
  try { reader->Update(); }
  catch(ExceptionObject &exc)
    { cerr << exc; return usage(); }

  // Create a radius size object
  ImageType::SizeType sz; sz.Fill(xRad);

  // Create a peak array
  vector< pair<float, Peak> > lPeaks;

  // Create a neighborhood iterator for searching around in the image
  ImageType::Pointer img = reader->GetOutput();
  ConstNeighborhoodIterator<ImageType> it(sz, img, img->GetBufferedRegion());
  for( ; !it.IsAtEnd(); ++it)
    {
    bool ismax = true;

    // Only admin finite pixels
    if(!finite(it.GetCenterPixel()) || it.GetCenterPixel() < xThresh) continue;
    
    // Check if the value is maximum
    for(size_t i = 0; i < it.GetNeighborhood().Size(); ++i)
      {
      // Skip the center pixel
      if(i == it.GetNeighborhood().GetCenterNeighborhoodIndex()) continue;
      
      bool inb = false;
      float z = it.GetPixel(i, inb);
      if(!isfinite(z) || (inb && z >= it.GetCenterPixel()))
        {
        ismax = false;
        break;
        }
      }

    if(ismax)
      {
      Peak peak;
      peak.index = it.GetIndex();
      
      // Adjust coordinate values
      img->TransformIndexToPhysicalPoint(peak.index, peak.point);
      for(size_t l = 0; l < 3; l++)
        peak.point[l] -= img->GetOrigin()[l] * 2;
      
      lPeaks.push_back(make_pair(it.GetCenterPixel(), peak));
      cerr << "." << flush;
      }
    }

  // Sort the peak list
  cerr << endl;
  sort(lPeaks.begin(), lPeaks.end());
  for(size_t k = 0; k < lPeaks.size(); k++)
    {
    cout << "Peak " << setw(5) << k << ": ";
    cout << "Value = " << setw(8) << lPeaks[k].first << "; ";
    cout << "Index = " << lPeaks[k].second.index << "; ";
    cout << "MNI = " << lPeaks[k].second.point << "; ";
    cout << "TAL = " << MNI2Tal(lPeaks[k].second.point) << ";" << endl;
    }

  // Open the talairach file
  if(fnTal.size() > 0)
    {
    ofstream fout(fnTal.c_str());
    for(size_t k = 0; k < lPeaks.size(); k++)
      {
      ImageType::PointType p = MNI2Tal(lPeaks[k].second.point);
      fout << p[0] << " " << p[1] << " " << p[2] << endl;
      }
    fout.close();
    }

  // Save the image if asked to
  if(fnImage.size() > 0) 
    {
    ImageType::Pointer peakimg = ImageType::New();
    peakimg->SetRegions(img->GetBufferedRegion());
    peakimg->Allocate();
    peakimg->SetOrigin(img->GetOrigin());
    peakimg->SetSpacing(img->GetSpacing());
    peakimg->FillBuffer(0.0f);

    for(size_t k = 0; k < lPeaks.size(); k++)
      {
      peakimg->SetPixel(lPeaks[k].second.index, lPeaks[k].first);
      }

    typedef ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput(peakimg);
    writer->SetFileName(fnImage.c_str());
    writer->Update();
    }
    
}

  
