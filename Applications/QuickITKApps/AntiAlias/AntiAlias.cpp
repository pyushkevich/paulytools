#include "itkAntiAliasBinaryImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkCommand.h"

#include "itkIdentityTransform.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"

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

void ProgressCommand(itk::Object *dummy, const itk::EventObject &, void *data)
{
  cout << "." << flush;
}

int usage()
{
  cout << "aalias - Apply ITK Whitaker's Anti-Aliasing method to a label image" << endl;
  cout << "usage:" << endl;
  cout << "   aalias [options] input.img output_base_name " << endl;
  cout << "options:" << endl;
  cout << "   -iso x.xx    Rescale image to be isotropic w/ voxel size x.xx" << endl;
  cout << "   -rms x.xx    Max RMS error parameter (0.07) " << endl;
  cout << "   -itr nn      Max iterations (omit: as many as needed) " << endl;
  cout << endl;
  return -1;
}

int main(int argc, char **argv)
{
  // Check number of parameters
  if(argc < 3) return usage();

  // Read in the parameters
  unsigned int nMaxIterations = 0;
  bool flagDoResample = false;
  double xMaxRMS = 0.07, xVoxelSize = 1.0;
  string fnInput = argv[argc-2], fnOutput = argv[argc-1];

  try 
    {
    for(unsigned int iArg = 1; iArg < argc -2; iArg++)
      {
      string arg = argv[iArg];
      if(arg == "-iso")
        {
        flagDoResample = true;
        xVoxelSize = atof(argv[++iArg]);
        }        
      else if(arg == "-rms")
        {
        xMaxRMS = atof(argv[++iArg]);
        }        
      else if(arg == "itr")
        {
        nMaxIterations = atoi(argv[++iArg]);
        }        
      else
        {
        cerr << "Bad parameter " << arg << endl;
        return usage();
        }
      }
    }
  catch(...)
    {
    cerr << "Error parsing parameters" << endl;
    return usage();
    }

  // Load the input image
  typedef itk::Image<short,3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  
  ReaderType::Pointer fltReader = ReaderType::New();
  ImageType::Pointer imgInput = ImageType::New();
  
  fltReader->SetFileName(fnInput.c_str());
  imgInput = fltReader->GetOutput();
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
      ++itInput;
    }

  // For each non-zero type, apply the antialias application
  for(vector<HistoEntry>::iterator it=xHistogram.begin();it!=xHistogram.end();it++)
    { 
    if(it->count > 0 && it->val > 0)
      {
      cout << "Processing intensity level " << it->val << endl;
      cout << "     Number of voxels        : " << it->count << endl << endl;
      
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
      typedef itk::Image<float, 3> FloatImageType;
      FloatImageType::Pointer imgBinary = FloatImageType::New();
      imgBinary->SetRegions(imgInput->GetBufferedRegion());
      imgBinary->SetSpacing(imgInput->GetSpacing());
      imgBinary->Allocate();
      imgBinary->FillBuffer(1.0f);

      // Pass the image through, binarizing current label
      itk::ImageRegionConstIterator<ImageType> itBinSource(imgInput,imgInput->GetBufferedRegion());
      itk::ImageRegionIterator<FloatImageType> itBinTarget(imgBinary,imgInput->GetBufferedRegion());
      while(!itBinSource.IsAtEnd())
        {
        itBinTarget.Set(it->val == itBinSource.Get() ? -1.0f : 1.0f);
        ++itBinSource;++itBinTarget;
        }

      // Create progress command
      itk::CStyleCommand::Pointer command = itk::CStyleCommand::New();
      command->SetCallback(ProgressCommand);

      // If requested, resample the binary image
      if(flagDoResample)
        {
        // Resample the image from 5mm z to 1mm z distance
        typedef itk::ResampleImageFilter<FloatImageType,FloatImageType> ResampleFilter;
        typedef itk::LinearInterpolateImageFunction<FloatImageType,double> IFType;
          // typedef itk::BSplineInterpolateImageFunction<FloatImageType,double> IFType;

        ResampleFilter::Pointer fltResample = ResampleFilter::New();
        fltResample->SetInput(imgBinary);
        fltResample->SetTransform(itk::IdentityTransform<double,3>::New());
        fltResample->SetInterpolator(IFType::New());

        // Compute the size
        itk::Size<3> szInput = imgBinary->GetBufferedRegion().GetSize();
        itk::Size<3> szOutput; 
        for(unsigned int j=0;j<3;j++)
          {
          szOutput[j] = 
            (unsigned long) (szInput[j] * imgBinary->GetSpacing()[j] / xVoxelSize);
          }
        fltResample->SetSize(szOutput);

        // Set the spacing
        double spacing[] = {xVoxelSize,xVoxelSize,xVoxelSize};
        fltResample->SetOutputSpacing(spacing);
        fltResample->SetDefaultPixelValue(1.0f);

        // Run the resampler
        fltResample->AddObserver(itk::ProgressEvent(),command);
        fltResample->UpdateLargestPossibleRegion();

        // Set the 'new' binary image 
        imgBinary = fltResample->GetOutput();
        }

      // Apply antialiasing to the binary image
      typedef itk::AntiAliasBinaryImageFilter<FloatImageType,FloatImageType>
        AntiFilterType;
      AntiFilterType::Pointer fltAnti = AntiFilterType::New();
      fltAnti->SetInput(imgBinary);
      fltAnti->SetMaximumRMSError(xMaxRMS);      
      if(nMaxIterations > 0)
        fltAnti->SetNumberOfIterations(nMaxIterations);
      fltAnti->SetIsoSurfaceValue(0.0f);
      fltAnti->AddObserver(itk::ProgressEvent(),command);
      fltAnti->Update();

      // Recast to short
      ImageType::Pointer imgOut = ImageType::New();
      imgOut->SetRegions(fltAnti->GetOutput()->GetBufferedRegion());
      imgOut->Allocate();
      imgOut->FillBuffer(0);

      itk::ImageRegionConstIterator<FloatImageType> itOutSource(fltAnti->GetOutput(),imgOut->GetBufferedRegion());
      itk::ImageRegionIterator<ImageType> itOutTarget(imgOut,imgOut->GetBufferedRegion());
      while(!itOutSource.IsAtEnd())
        {
        itOutTarget.Set(itOutSource.Get() < 0.0 ? it->val : 0);
        ++itOutTarget;++itOutSource;
        }

      // Save the output
      std::ostringstream oss;
      oss << argv[2] << "." << it->val << ".hdr";
      
      typedef itk::ImageFileWriter<ImageType> WriterType;
      WriterType::Pointer fltWriter = WriterType::New();
      fltWriter->SetInput(imgOut);
      fltWriter->SetFileName(oss.str().c_str());
      fltWriter->Update();
      }
    }
}
