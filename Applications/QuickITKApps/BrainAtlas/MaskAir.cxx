#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkFFTWRealToComplexConjugateImageFilter.h"
#include "itkFFTWComplexConjugateToRealImageFilter.h"

#include <iostream>

using namespace itk;
using namespace std;

// Global type definitions
typedef Image<unsigned char, 3> ByteImageType;
typedef Image<float, 3> FloatImageType;
typedef ImageFileReader<FloatImageType> ReaderType;
typedef ImageFileWriter<ByteImageType> WriterType;
typedef std::complex<float> cplex;
typedef Image<float, 1> Image1f;
typedef Image<cplex, 1> Image1c; 
typedef Image<cplex, 3> Image3c; 

int usage()
{
  cout << "maskair - mask air/background pixels in an MRI image" << endl;
  cout << "usage: " << endl;
  cout << "  maskair [options] input.img output.img " << endl;
  cout << "options: " << endl;
  cout << "  -t X.XX        Threshold intensity for air. Default: 75" << endl;
  cout << "  -r N           Radius of the box filter. Default: 16" << endl;
  return -1;
}


int main(int argc, char *argv[])
{
  // Must have enough parameters
  if(argc < 3) return usage();
  
  // These are the parameters that users can set
  const char *fnIn, *fnOut;
  double xThresh = 75.0;
  unsigned int iRadius = 10;

  // Process user input
  fnIn = argv[argc - 2];
  fnOut = argv[argc - 1];
  for(unsigned int iArg = 1; iArg < argc-2; iArg++)
    {
    if(strcmp(argv[iArg],"-t") == 0)
      xThresh = atof(argv[++iArg]);
    else if(strcmp(argv[iArg],"-r") == 0)
      iRadius = atoi(argv[++iArg]);
    else
      cerr << "Unknown parameter " << argv[iArg] << endl;
    }

  // Read the input image
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(fnIn);
  FloatImageType::Pointer imgInput = fltReader->GetOutput();
  try 
    {
    fltReader->Update();
    }
  catch(...)
    {
    cerr << "Error while trying to load the input image. Aborting!";
    return -1;
    }

  // Generate the mask image, fill it with ones.
  FloatImageType::Pointer imgToThresh = FloatImageType::New();
  imgToThresh->SetRegions(imgInput->GetBufferedRegion());
  imgToThresh->Allocate();
  imgToThresh->SetOrigin(imgInput->GetOrigin());
  imgToThresh->SetSpacing(imgInput->GetSpacing());
  imgToThresh->FillBuffer(1);

  // Define the FFT filters in 3D
  typedef FFTWRealToComplexConjugateImageFilter<float, 3> FFTFilter3d;
  typedef FFTWComplexConjugateToRealImageFilter<float, 3> FFTInvFilter3d;
  FFTFilter3d::Pointer fltFFT = FFTFilter3d::New();
  FFTInvFilter3d::Pointer fltInvFFT = FFTInvFilter3d::New();

  // Compute the FFT of the input image
  cout << "Computing the FFT of the input image ... " << flush;
  fltFFT->SetInput(imgInput);
  fltFFT->Update();
  Image3c::Pointer imgFFT = fltFFT->GetOutput();
  cout << "done" << endl;

  // Create a kernel image
  FloatImageType::Pointer imgKernel = FloatImageType::New();
  imgKernel->SetRegions(imgInput->GetBufferedRegion());
  imgKernel->Allocate();

  // There are eight kernels that will be computed.
  for(unsigned int z = 0; z < 2; z++)
    {
    for(unsigned int y = 0; y < 2; y++)
      {
      for(unsigned int x = 0; x < 2; x++)
        {
        // Zero out the kernel
        imgKernel->FillBuffer(0.0);

        // Compute the region where the kernel is set to ones
        FloatImageType::RegionType regBox;
        regBox.SetIndex(0, x * (imgInput->GetBufferedRegion().GetSize(0) - iRadius));
        regBox.SetIndex(1, y * (imgInput->GetBufferedRegion().GetSize(1) - iRadius));
        regBox.SetIndex(2, z * (imgInput->GetBufferedRegion().GetSize(2) - iRadius));
        regBox.SetSize(0,iRadius); regBox.SetSize(1,iRadius); regBox.SetSize(2,iRadius);
        
        // Create an iterator for the kernel image
        typedef ImageRegionIterator<FloatImageType> FloatIter;
        FloatIter itKernel(imgKernel, regBox);
        float xKernelValue = 1.0 / regBox.GetNumberOfPixels();
        while(!itKernel.IsAtEnd())
          {
          itKernel.Set(xKernelValue);
          ++itKernel;
          }

        // Compute the FFT of the kernel
        cout << "Computing the FFT of the kernel [" << x << y << z << "] ... " << flush;
        FFTFilter3d::Pointer fltKernelFFT = FFTFilter3d::New();
        fltKernelFFT->SetInput(imgKernel);
        fltKernelFFT->Update();
        Image3c::Pointer imgKernelFFT = fltKernelFFT->GetOutput();
        cout << "done" << endl;
        
        // Iterate through the input image
        typedef ImageRegionIteratorWithIndex<Image3c> Iter3c;
        Iter3c it1(imgFFT, imgFFT->GetBufferedRegion());
        Iter3c it2(imgKernelFFT, imgKernelFFT->GetBufferedRegion());
        while(!it2.IsAtEnd())
          {
          // Multiply the two values
          cplex xInput = it1.Value();
          cplex xKernel = it2.Value();
          it2.Set(xInput * xKernel);
          ++it1; ++it2;
          }

        // Compute the inverse FFT
        cout << "Computing the Inverse FFT for [" << x << y << z << "] ... " << flush;
        fltInvFFT->SetInput(imgKernelFFT);
        fltInvFFT->Update();
        cout << "done" << endl;

        // Integrate the image with the output image
        typedef ImageRegionIterator<FloatImageType> Iter3f;
        Iter3f iz1(fltInvFFT->GetOutput(), imgToThresh->GetBufferedRegion());
        Iter3f iz2(imgToThresh, imgToThresh->GetBufferedRegion());
        
        // If this is the first iteration in the loop, we just copy the image
        if(x == 0 && y == 0 && z == 0)
          {
          while(!iz2.IsAtEnd())
            {
            iz2.Set(iz1.Value());
            ++iz1; ++iz2;
            }
          }
        else
          {
          while(!iz2.IsAtEnd())
            {
            float a = iz1.Value();
            float b = iz2.Value();
            if( a < b )
              iz2.Set( a );
            ++iz1; ++iz2;
            }
          }
        }
      }
    }
  
  // Finally, compute the threshhold
  typedef BinaryThresholdImageFilter<FloatImageType, ByteImageType> ThreshFilter;
  ThreshFilter::Pointer fltThresh = ThreshFilter::New();
  fltThresh->SetInput(imgToThresh);
  fltThresh->SetUpperThreshold(xThresh);
  fltThresh->SetInsideValue(1);
  fltThresh->SetOutsideValue(2);
  fltThresh->Update();

  // Write the result
  WriterType::Pointer fltWriter = WriterType::New();
  fltWriter->SetInput(fltThresh->GetOutput());
  fltWriter->SetFileName(fnOut);
  fltWriter->Update();
}
