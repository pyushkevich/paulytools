/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkResampleImageFilter.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

using namespace itk;

int usage()
{
  cout << "sincfilter2d [options] input.png output.png" << endl;
  cout << "options: " << endl;
  cout << "   -s X.XX                Scale image isotropically by X.XX " << endl;
  cout << "   -as X.XX Y.YY          Anisotropically scale by X.XX Y.YY " << endl;
  cout << "   --blackman             Use Blackman window function" << endl;
  cout << "   --cosine               Use cosine window function" << endl;
  cout << "   --hamming              Use Hamming window function" << endl;
  cout << "   --lancos               Use Lancos window function" << endl;
  cout << "   --welch                Use Welch window function" << endl;
  cout << "   --linear               Use linear interpolation instead" << endl;
  cout << endl;
  return -1;
}

class MapFunction {
public:
  unsigned char operator()(double x)
    {
    if(x < 0.0) return 0;
    else if(x >= 256.0) return 255;
    else return (unsigned char) x;
    }
};

typedef itk::Image<unsigned char, 2> ImageType;
typedef itk::Image<float, 2> FloatImageType;

template <class TInterpolationFunction> 
void 
ResampleImage(
  ImageType *imgInput, const char *fnOutput, 
  double xScale, double yScale, TInterpolationFunction *fInterp)
{
  // Compute the new size and the new spacing of the output image
  ImageType::SizeType szInput = imgInput->GetBufferedRegion().GetSize();
  ImageType::SizeType szOutput;
  szOutput[0] = (unsigned long) (szInput[0] * xScale);
  szOutput[1] = (unsigned long) (szInput[1] * yScale);
  
  ImageType::SpacingType xSpacing = imgInput->GetSpacing();
  xSpacing[0] /= xScale;
  xSpacing[1] /= yScale;

  // Create a filter for resampling the image
  typedef itk::ResampleImageFilter<ImageType,FloatImageType> ResampleFilterType;
  ResampleFilterType::Pointer fltSample = ResampleFilterType::New();

  // Initialize the resampling filter
  fltSample->SetInput(imgInput);
  fltSample->SetTransform(itk::IdentityTransform<double,2>::New());

  // Assign the interpolation function to the filter
  fltSample->SetInterpolator(fInterp);

  // Set the image sizes and spacing
  fltSample->SetSize(szOutput);
  fltSample->SetOutputSpacing(xSpacing);
  fltSample->SetOutputOrigin(imgInput->GetOrigin());

  // Set the progress bar
  //if(progressCommand)
  //  fltSample->AddObserver(itk::AnyEvent(),progressCommand);

  // Perform resampling
  fltSample->Update();
  
  // Convert the resulting floating point image to the output image
  typedef UnaryFunctorImageFilter<FloatImageType, ImageType, MapFunction>
    RemapFilter;
  RemapFilter::Pointer fltRemap = RemapFilter::New();
  fltRemap->SetInput(fltSample->GetOutput());
  fltRemap->Update();

  // Save the image
  typedef ImageFileWriter<ImageType> WriterFilter;
  WriterFilter::Pointer fltWriter = WriterFilter::New();
  fltWriter->SetInput(fltRemap->GetOutput());
  fltWriter->SetFileName(fnOutput);
  fltWriter->Update();
}

int main(int argc, char *argv[])
{
  // Operation modes
  enum Mode {LINEAR, BLACKMAN, COSINE, HAMMING, LANCOS, WELCH};
  Mode mode = HAMMING;
  
  // Check parameters
  if(argc < 3) return usage();
  const char *fnInput = argv[argc - 2];
  const char *fnOutput = argv[argc - 1];
  double xScale = 2.0, yScale = 2.0;
  bool flagLinear = false;

  // Check options
  for(unsigned int iArg=1;iArg<argc-2;iArg++)
    {
    if(!strcmp(argv[iArg], "-s"))
      {
      xScale = yScale = atof(argv[++iArg]);
      }
    else if(!strcmp(argv[iArg], "-as"))
      {
      xScale = atof(argv[++iArg]);
      yScale = atof(argv[++iArg]);
      }
    else if(!strcmp(argv[iArg], "--blackman"))
      mode = BLACKMAN;
    else if(!strcmp(argv[iArg], "--cosine"))
      mode = COSINE;
    else if(!strcmp(argv[iArg], "--hamming"))
      mode = HAMMING;
    else if(!strcmp(argv[iArg], "--lancos"))
      mode = LANCOS;
    else if(!strcmp(argv[iArg], "--welch"))
      mode = WELCH;
    else if(!strcmp(argv[iArg], "--linear"))
      mode = LINEAR;
    else 
      {
      cerr << "Unknown option " << argv[iArg] << endl;
      return usage();
      }
    }
  
  // Load the input image
  typedef ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(fnInput);
  fltReader->Update();
  ImageType::Pointer imgInput = fltReader->GetOutput();
    
  // Create the interpolation function
  const unsigned int xRadius = 6;
  typedef ConstantBoundaryCondition<ImageType> BoundaryCondition;

  if(mode != LINEAR)
    {
    if(mode == BLACKMAN)
      {
      typedef WindowedSincInterpolateImageFunction< ImageType, xRadius, 
         Function::BlackmanWindowFunction<xRadius>, 
         BoundaryCondition, double> FunctionType;
      FunctionType::Pointer fInterp = FunctionType::New();
      ResampleImage(imgInput, fnOutput, xScale, yScale, fInterp.GetPointer());
      }

    else if(mode == COSINE)
      {
      typedef WindowedSincInterpolateImageFunction< ImageType, xRadius, 
         Function::CosineWindowFunction<xRadius>, 
         BoundaryCondition, double> FunctionType;
      FunctionType::Pointer fInterp = FunctionType::New();
      ResampleImage(imgInput, fnOutput, xScale, yScale, fInterp.GetPointer());
      }
    else if(mode == HAMMING)
      {
      typedef WindowedSincInterpolateImageFunction< ImageType, xRadius, 
         Function::HammingWindowFunction<xRadius>, 
         BoundaryCondition, double> FunctionType;
      FunctionType::Pointer fInterp = FunctionType::New();
      ResampleImage(imgInput, fnOutput, xScale, yScale, fInterp.GetPointer());
      }
    else if(mode == LANCOS)
      {
      typedef WindowedSincInterpolateImageFunction< ImageType, xRadius, 
         Function::LancosWindowFunction<xRadius>, 
         BoundaryCondition, double> FunctionType;
      FunctionType::Pointer fInterp = FunctionType::New();
      ResampleImage(imgInput, fnOutput, xScale, yScale, fInterp.GetPointer());
      }
    else if(mode == WELCH)
      {
      typedef WindowedSincInterpolateImageFunction< ImageType, xRadius, 
         Function::WelchWindowFunction<xRadius>, 
         BoundaryCondition, double> FunctionType;
      FunctionType::Pointer fInterp = FunctionType::New();
      ResampleImage(imgInput, fnOutput, xScale, yScale, fInterp.GetPointer());
      }
    }
  else
    {
    typedef LinearInterpolateImageFunction<ImageType, double>
      FunctionType;
    FunctionType::Pointer fInterp = FunctionType::New();

    // Perform interpolation using that function
    ResampleImage(imgInput, fnOutput, xScale, yScale, fInterp.GetPointer());
    }
}

