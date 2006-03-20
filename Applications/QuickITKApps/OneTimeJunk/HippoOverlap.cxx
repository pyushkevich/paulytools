#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#include <iostream>

using namespace std;
using namespace itk;

int usage()
{
  cout << "  hippov: Computes overlap between segmentations in shared space " << endl;
  cout << "  usage: " << endl;
  cout << "    hippov [options] statmap.img" << endl;
  cout << "  options; " << endl;
  cout << "    -i img1 .. imgN    Hippo segmentation image list" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  if(argc < 4) return usage();

  // Process the command line options
  string fnStatmap = argv[argc-1];
  vector<string> fnInput;
  string sLastOption = "";
  for(size_t k = 1; k < argc-1; k++)
    {
    if(argv[k][0] == '-')
      if(argv[k][1] == 'i')
        sLastOption = "i";
      else
        return usage();
    else if(sLastOption == "i")
      fnInput.push_back(argv[k]);
    else
      return usage();
    }

  // Print report
  size_t n = fnInput.size();
  cout << "Performing analysis on " << n << " images" 
    << "; statmap is " << fnStatmap << endl;

  // Typedefs
  typedef itk::Image<float, 3> ImageType;
  typedef ImageType::Pointer ImagePointer;
  typedef itk::ImageFileReader<ImageType> ImageReader;
  
  // Read the input images
  vector<ImagePointer> imgInput;
  for(size_t i = 0; i < n; i++)
    {
    // Read the input image
    ImageReader::Pointer fltReader = ImageReader::New();
    fltReader->SetFileName(fnInput[i].c_str());
    fltReader->Update();
    imgInput.push_back(fltReader->GetOutput());
    }

  // Read the statmap image
  ImageReader::Pointer fltReader = ImageReader::New();
  fltReader->SetFileName(fnStatmap.c_str());
  fltReader->Update();
  ImagePointer imgStat = fltReader->GetOutput();
  
  // Create an interpolator
  typedef itk::NearestNeighborInterpolateImageFunction<
    ImageType, double> InterpolatorType;
  InterpolatorType::Pointer fInterp = InterpolatorType::New();
  fInterp->SetInputImage(imgStat);

  // Compute the overlap and union area maximal values
  size_t nOverlap = 0, nUnion = 0;
  float maxOverlap = 0.0, maxUnion = 0.0;

  itk::ImageRegionIteratorWithIndex<ImageType> it(
    imgInput[0], imgInput[0]->GetBufferedRegion());
  
  for(; !it.IsAtEnd(); ++it)
    {
    bool flagAnyPixel = false, flagAllPixels = true;
    ImageType::IndexType idx = it.GetIndex();
    for(size_t i = 0; i < n; i++)
      {
      float x = imgInput[i]->GetPixel(idx);
      if(x > 0) flagAnyPixel = true;
      else flagAllPixels = false;
      }

    ImageType::PointType point;
    itk::ContinuousIndex<double, 3> cidx;
    imgInput[0]->TransformIndexToPhysicalPoint(idx, point);
    imgStat->TransformPhysicalPointToContinuousIndex(point, cidx);
    float val = fInterp->EvaluateAtContinuousIndex(cidx);
    
    if(flagAnyPixel)
      { nUnion++; if(val > maxUnion) maxUnion = val; }
    if(flagAllPixels)
      { nOverlap++; if(val > maxOverlap) maxOverlap = val; }
    }
    
  cout << "Overlap : " << nOverlap << endl;
  cout << "Union : " << nUnion << endl;

  cout << "tMax Overlap : " << maxOverlap << endl;
  cout << "tMax Union : " << maxUnion << endl;
}

