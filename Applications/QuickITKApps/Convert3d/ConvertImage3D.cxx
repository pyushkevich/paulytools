// We define manual instantiation to speed up compilation
#define ITK_MANUAL_INSTANTIATION 1

// ITK includes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkExtractImageFilter.h"
#include "itkByteSwapper.h"
#include "itkMetaDataObject.h"
#include "itkIOCommon.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h" 
#include "itkAntiAliasBinaryImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkVoxBoCUBImageIOFactory.h"
#include "itkPovRayDF3ImageIOFactory.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include <iostream>
#include <cctype>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

template<class TPixel, unsigned int VDim>
class ImageConverter
{
public:
  ImageConverter();
  int ProcessCommandLine(int argc, char *argv[]);

private:
  // Image typedef
  typedef itk::Image<TPixel, VDim> ImageType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef typename ImageType::SizeType SizeType;
  typedef vnl_vector_fixed<double, VDim> RealVector;
  typedef vnl_vector_fixed<int, VDim> IntegerVector;
  
  // Internal functions
  void ReadImage(const char *file);
  void WriteImage(const char *file);
  void AntiAliasImage(double iso);
  void CopyImage();
  void ScaleShiftImage(double a, double b);
  void PrintImageInfo(bool flagFullInfo);
  void ThresholdImage(double u1, double u2, double vIn, double vOut);
  void TrimImage(const RealVector &margin);
  int ProcessCommand(int argc, char *argv[]);
  int ProcessResampleCommand(int argc, char *argv[]);
  int ProcessSmoothCommand(int argc, char *argv[]);

  // Read vectors, etc from command line
  SizeType ReadSizeVector(char *vec);
  RealVector ReadRealVector(char *vec);

  // Get bounding box of an image
  void GetBoundingBox(ImageType *image, RealVector &bb0, RealVector &bb1);

  // Templated write function
  template<class TOutPixel> void TemplatedWriteImage(const char *file, double xRoundFactor);

  // Stack of images from the command line
  vector<ImagePointer> m_ImageStack;

  // Typeid of the image to be saved
  string m_TypeId;

  // Interpolation type
  string m_Interpolation;

  // Background intensity
  double m_Background;

  // Rounding factor (not used for float and double) is 0.5 unless -noround was specified
  double m_RoundFactor;

  // Whether SPM extensions are used
  bool m_FlagSPM;

  // Number of iterations for various algorithms
  size_t m_Iterations;

  // Root mean square error for anti-aliasing algorithm
  double m_AntiAliasRMS;

  // Verbose output stream
  std::ostringstream devnull;
  std::ostream *verbose;
};

template<class TPixel, unsigned int VDim>
ImageConverter<TPixel,VDim>
::ImageConverter()
  : verbose(&devnull)
{
  // Initialize to defaults
  m_TypeId = "float";
  m_Interpolation = "Linear";
  m_Background = 0.0;
  m_RoundFactor = 0.5;
  m_FlagSPM = false;
  m_AntiAliasRMS = 0.07;
  m_Iterations = 0;
}

template<class TPixel, unsigned int VDim>
int
ImageConverter<TPixel, VDim>
::ProcessCommand(int argc, char *argv[])
{
  // Get the first command
  string cmd = argv[0];

  // Anti-alias a binary image, turning it into a smoother floating point image;
  // the argument is the iso-surface value
  // This command is affected by -iterations and -rms flags
  if(cmd == "-antialias" || cmd == "-alias")
    {
    AntiAliasImage(atof(argv[1]));
    return 1;
    }
  
  else if(cmd == "-background")
    { m_Background = atof(argv[1]); return 1; }

  // f(x) = (x == xBackground) ? 0 : 1
  else if(cmd == "-binarize")
    {
    ThresholdImage(m_Background, m_Background, 0.0, 1.0);
    return 0;
    }

  else if(cmd == "-info")
    { PrintImageInfo(false); return 0; }

  else if(cmd == "-info-full")
    { PrintImageInfo(true); return 0; }

  else if(cmd == "-interpolation")
    { m_Interpolation = argv[1]; return 1; }

  else if(cmd == "-iterations")
    { m_Iterations = static_cast<size_t>(atoi(argv[1])); return 1; }

  else if(cmd == "-noround")
    { m_RoundFactor = 0.0; return 0; }

  // Enable SPM extensions
  else if(cmd == "-nospm")
    { m_FlagSPM = false; return 0; }

  // Resample command
  else if(cmd == "-resample")
    return ProcessResampleCommand(argc-1, argv+1);

  else if(cmd == "-rms")
    { m_AntiAliasRMS = atof(argv[1]); return 1; }
 
  else if(cmd == "-round")
    { m_RoundFactor = 0.5; return 0; }

  else if(cmd == "-scale")
    {
    double factor = atof(argv[1]);
    ScaleShiftImage(factor, 0.0);
    return 1;
    }

  // Select command: push image on the stack
  /*
  else if (cmd == "-select")
    {
    int number = atoi(argv[1]);
    if(number == 0)
      { cerr << "Wrong parameter for -select!" << endl; throw -1; }
    if(abs(number) > m_ImageStack.size())
      { cerr << "Parameter for -select greater than stack size!" << endl; throw -1; }
    if(number > 0)
      m_ImageStack.push_back(m_ImageStack[number-1]);
    else 
      m_ImageStack.push_back(m_ImageStack[m_ImageStack.size()+number]);
    return 1;
    }
  */

  else if(cmd == "-shift")
    {
    double x = atof(argv[1]);
    ScaleShiftImage(1.0, x);
    return 1;
    }

  // Gaussian smoothing command
  else if(cmd == "-smooth")
    return ProcessSmoothCommand(argc-1, argv+1);

  // Enable SPM extensions
  else if(cmd == "-spm")
    { m_FlagSPM = true; return 0; }

  // Stretch the intensity range
  else if(cmd == "-stretch")
    {
    double u1 = atof(argv[1]);
    double u2 = atof(argv[2]);
    double v1 = atof(argv[3]);
    double v2 = atof(argv[4]);
    double a = (v2 - v1) / (u2 - u1);
    double b = v1 - u1 * a;
    ScaleShiftImage(a, b);
    return 4;
    }

  // Thresholding
  else if(cmd == "-threshold" || cmd == "-thresh")
    {
    double u1 = strcmp(argv[1],"-inf") == 0 ? -vnl_huge_val(0.0) : atof(argv[1]);
    double u2 = strcmp(argv[2],"inf") == 0 ? vnl_huge_val(0.0) : atof(argv[2]);
    double v1 = atof(argv[3]);
    double v2 = atof(argv[4]);
    ThresholdImage(u1, u2, v1, v2);
    return 4;
    }

  // Trim the image (trim background values from the margin)
  else if(cmd == "-trim")
    {
    // Read the size of the wrap region
    RealVector margin = ReadRealVector(argv[1]);

    // Trim the image accordingly
    TrimImage(margin);

    // Return the number of arguments consumed
    return 1;
    }

  // Output type specification
  else if(cmd == "-type")
    { 
    m_TypeId = argv[1]; 
    int (*pf)(int) = tolower;
    transform(m_TypeId.begin(), m_TypeId.end(), m_TypeId.begin(), pf);
    return 1; 
    }

  // Verbose mode
  else if(cmd == "-verbose")
    { verbose = &std::cout; return 0; }

  // Unknown command
  else
    { cerr << "Unknown command " << cmd << endl; throw -1; }

  cerr << "Fell through!" << endl;
  throw -1;
}


template<class TPixel, unsigned int VDim>
int
ImageConverter<TPixel, VDim>
::ProcessCommandLine(int argc, char *argv[])
{
  // The last parameter in the command line is the output file
  string fnOutput = argv[argc-1];
  
  // Command line processing
  for(size_t i = 1; i < argc; ++i)
    {
    string cmd = argv[i];
    if(cmd[0] == '-')
      {
      // A command has been supplied
      i += ProcessCommand(argc-(i+1), argv+i);
      }
    else
      {
      // An image file name has been provided. If this image is followed by commands
      // read it and push in the pipeline.
      if(i != argc-1)
        ReadImage(argv[i]);
      else
        WriteImage(argv[i]);
      }
    }
}

template<class TPixel, unsigned int VDim> 
template<class TOutPixel>
void 
ImageConverter<TPixel, VDim>
::TemplatedWriteImage(const char *file, double xRoundFactor)
{
  // Get the input image
  if(m_ImageStack.size() == 0)
    { cerr << "No data has been generated! Can't write to " << file << endl; throw -1; }
  ImagePointer input = m_ImageStack.back();
  
  // Create the output image 
  typedef itk::Image<TOutPixel, VDim> OutputImageType;
  typename OutputImageType::Pointer output = OutputImageType::New();
  output->SetRegions(input->GetBufferedRegion());
  output->SetSpacing(input->GetSpacing());
  output->SetOrigin(input->GetOrigin());
  output->SetMetaDataDictionary(input->GetMetaDataDictionary());
  output->Allocate();

  // Describe what we are doing
  *verbose << "Writing #" << m_ImageStack.size() << " to file " << file << endl;
  *verbose << "  Output voxel type: " << m_TypeId << "[" << typeid(TOutPixel).name() << "]" << endl;
  *verbose << "  Rounding off: " << (xRoundFactor == 0.0 ? "Disabled" : "Enabled") << endl;
  
  // Set the SPM originator header
  if(m_FlagSPM)
    {
    size_t i;
    string originator;
    originator.resize(VDim * 2);
    
    // Compute the SPM-style origin of the image
    *verbose << "  Setting SPM origin field to:";
    for(i = 0; i < VDim; i++)
      {
      short ospm = (short)(0.5 - input->GetOrigin()[i] / input->GetSpacing()[i]);
      originator[2*i] = (char)(ospm & 0x00ff);
      originator[2*i+1] = (char)(ospm >> 8);
      *verbose << ospm << " ";
      }
    originator[2*i] = 0;
    *verbose << endl;

    itk::EncapsulateMetaData<string>(
      output->GetMetaDataDictionary(),itk::ITK_FileOriginator,originator);
    }

  // Copy everything, rounding if the pixel type is integer
  size_t n = input->GetBufferedRegion().GetNumberOfPixels();
  for(size_t i = 0; i < n; i++)
    output->GetBufferPointer()[i] = (TOutPixel) (input->GetBufferPointer()[i] + xRoundFactor);

  // Write the image out
  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput(output);
  writer->SetFileName(file);
  try { writer->Update(); }
  catch (itk::ExceptionObject &exc) {
    cerr << "Error writing image to " << file << endl;
    cerr << "ITK Exception: " << exc << endl;
    throw -1;
  }
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::WriteImage(const char *file)
{
  if(m_TypeId == "char" || m_TypeId == "byte")
    TemplatedWriteImage<char>(file, m_RoundFactor);
  if(m_TypeId == "uchar" || m_TypeId == "ubyte")
    TemplatedWriteImage<unsigned char>(file, m_RoundFactor);
  
  if(m_TypeId == "short") 
    TemplatedWriteImage<short>(file, m_RoundFactor);
  if(m_TypeId == "ushort")
    TemplatedWriteImage<unsigned short>(file, m_RoundFactor);

  if(m_TypeId == "int") 
    TemplatedWriteImage<int>(file, m_RoundFactor);
  if(m_TypeId == "uint")
    TemplatedWriteImage<unsigned int>(file, m_RoundFactor);

  if(m_TypeId == "float") 
    TemplatedWriteImage<float>(file, 0.0);
  if(m_TypeId == "double")
    TemplatedWriteImage<double>(file, 0.0);
}
  
template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::ReadImage(const char *file)
{
  // Set up the reader
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(file);

  // Report
  *verbose << "Reading #" << (1 + m_ImageStack.size()) << " from " << file << endl;
  
  // Try reading this file
  try { reader->Update(); }
  catch(itk::ExceptionObject &exc)
    {
    cerr << "Error reading image: " << file << endl;
    cerr << "ITK Exception: " << exc << endl;
    throw -1;
    }
  ImagePointer image = reader->GetOutput();

  // See if the originator field is present
  if(m_FlagSPM)
    {
    string temp;
    if(itk::ExposeMetaData<std::string>(
        image->GetMetaDataDictionary(), itk::ITK_FileOriginator, temp))
      {
      typename ImageType::PointType oitk = image->GetOrigin();
      typename ImageType::SpacingType sitk = image->GetSpacing();
      
      // Read the SPM-style origin
      *verbose << "  Applying SPM origin :";
      for(size_t i=0; i < VDim; i++)
        {
        short xspm = (temp[2*i+1] << 8) + temp[2*i];
        oitk[i] = -sitk[i] * xspm;
        *verbose << xspm << " ";
        }
      *verbose << endl;

      image->SetOrigin(oitk);
      }

    }

  // Push the image into the stack
  m_ImageStack.push_back(image);
}
bool str_ends_with(const std::string &s, const std::string &pattern)
{
  size_t ipos = s.rfind(pattern);
  return(ipos == s.size() - pattern.size());
}

template<class TPixel, unsigned int VDim>
typename ImageConverter<TPixel, VDim>::RealVector
ImageConverter<TPixel, VDim>
::ReadRealVector(char *vec)
{
  size_t i;
  string vecbk = vec;

  // Output vector
  RealVector x, scale;

  // Check the type of the vector
  if(str_ends_with(vec,"mm")) 
    scale.fill(1.0);
  else if(str_ends_with(vec,"vox"))
    for(i = 0; i < VDim; i++)
      scale[i] = m_ImageStack.back()->GetSpacing()[i];
  else
    { 
    cerr << "Invalid vector specification " << vec 
      << " (does not terminate with 'mm' or 'vox')" << endl;
    throw -1;
    }
  
  // Find all the 'x' in the string
  char *tok = strtok(vec, "xmvo");
  for(i = 0; i < VDim && tok != NULL; i++)
    {
    x[i] = atof(tok);
    tok = strtok(NULL, "xmvo");
    } 

  // Check if there is only one number
  if(i == 1)
    x.fill(x[0]);
  else if(i != VDim)
    {
    cerr << "Invalid vector specification " << vecbk
      << " (must have 1 or " << VDim << " components)" << endl;
    throw -1;
    }

  // Scale the vector
  for(i = 0; i < VDim; i++)
    x[i] *= scale[i];

  return x;
}

template<class TPixel, unsigned int VDim>
typename ImageConverter<TPixel, VDim>::ImageType::SizeType
ImageConverter<TPixel, VDim>
::ReadSizeVector(char *vec)
{
  string vecbk = vec;
  size_t i;

  typename ImageType::SizeType sz;

  // Check if the string ends with %
  if(str_ends_with(vec, "%"))
    {
    // Read floating point size
    RealVector factor;
    char *tok = strtok(vec, "x%");
    for(i = 0; i < VDim && tok != NULL; i++)
      {
      factor[i] = atof(tok);
      if(factor[i] <= 0)
        { cerr << "Non-positive percent size specification: " << vecbk << endl; throw -1; }      
      tok = strtok(NULL, "x%");
      } 

    if(i == 1)
      factor.fill(factor[0]);

    // Get the size of the image in voxels
    for(size_t i = 0; i < VDim; i++)
      sz[i] = (unsigned long)(m_ImageStack.back()->GetBufferedRegion().GetSize(i) * 0.01 * factor[i] + 0.5);
    }
  else
    {
    // Find all the 'x' in the string
    char *tok = strtok(vec, "x");
    for(size_t i = 0; i < VDim; i++)
      {
      if(tok == NULL)
        { cerr << "Invalid size specification: " << vecbk << endl; throw -1; }      
      int x = atoi(tok);
      if(x <= 0)
        { cerr << "Non-positive size specification: " << vecbk << endl; throw -1; }      
      sz[i] = (unsigned long)(x);
      tok = strtok(NULL, "x");
      } 
    }

  
  

  return sz;
}

template<class TPixel, unsigned int VDim>
int
ImageConverter<TPixel, VDim>
::ProcessSmoothCommand(int argc, char *argv[])
{
  // Get the input image
  ImagePointer input = m_ImageStack.back();
  
  // Read the parameter vector
  RealVector stdev = ReadRealVector(argv[0]);

  // Describe what we are doing
  *verbose << "Smoothing #" << m_ImageStack.size() << " with std.dev. " << stdev << endl;

  // Create a smoothing kernel and use it
  typedef itk::DiscreteGaussianImageFilter<ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  typename FilterType::ArrayType variance;

  for(size_t i = 0; i < VDim; i++)
    variance[i] = stdev[i] * stdev[i];
  
  filter->SetInput(input);
  filter->SetVariance(variance);
  filter->SetUseImageSpacingOn();
  filter->Update();

  // Save the output
  m_ImageStack.pop_back();
  m_ImageStack.push_back(filter->GetOutput());

  // Return the number of parameters used
  return 1;
}

template<class TPixel, unsigned int VDim>
int
ImageConverter<TPixel, VDim>
::ProcessResampleCommand(int argc, char *argv[])
{
  // Read the parameter vector
  typename ImageType::SizeType sz = ReadSizeVector(argv[0]);

  // Get the input image
  ImagePointer input = m_ImageStack.back();
  
  // Build the resampling filter
  typedef itk::ResampleImageFilter<ImageType,ImageType> ResampleFilterType;
  typename ResampleFilterType::Pointer fltSample = ResampleFilterType::New();

  // Initialize the resampling filter with an identity transform
  fltSample->SetInput(input);
  fltSample->SetTransform(itk::IdentityTransform<double,VDim>::New());

  // Typedefs for interpolators
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double> NNInterpolatorType;
  typedef itk::LinearInterpolateImageFunction<ImageType,double> LinearInterpolatorType;
  typedef itk::BSplineInterpolateImageFunction<ImageType,double> CubicInterpolatorType;

  // Create the interpolator depending on parameter
  if(m_Interpolation == "NearestNeighbor")
    fltSample->SetInterpolator(NNInterpolatorType::New());
  else if(m_Interpolation == "Linear")
    fltSample->SetInterpolator(LinearInterpolatorType::New());
  else if(m_Interpolation == "Cubic")
    fltSample->SetInterpolator(CubicInterpolatorType::New());
  else
    { cerr << "Unknown interpolation type: " << m_Interpolation << endl; throw -1; }

  // Compute the spacing of the new image
  typename ImageType::SpacingType spc = input->GetSpacing();
  for(size_t i = 0; i < VDim; i++)
    spc[i] *= input->GetBufferedRegion().GetSize()[i] * 1.0 / sz[i];

  // Get the bounding box of the input image
  RealVector bb0, bb1;
  GetBoundingBox(input, bb0, bb1);
  
  // Set the image sizes and spacing.
  fltSample->SetSize(sz);
  fltSample->SetOutputSpacing(spc);
  fltSample->SetOutputOrigin(bb0.data_block());

  // Set the unknown intensity to positive value
  fltSample->SetDefaultPixelValue(m_Background);

  // Describe what we are doing
  *verbose << "Resampling #" << m_ImageStack.size() << " to have" << sz << " voxels." << endl;
  *verbose << "  Interpolation method: " << m_Interpolation << endl;
  *verbose << "  Background intensity: " << m_Background << endl;

  // Perform resampling
  fltSample->UpdateLargestPossibleRegion();
    
  // Change the source to the output 
  m_ImageStack.pop_back();
  m_ImageStack.push_back(fltSample->GetOutput());

  // Return the number of parameters read
  return 1;
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::CopyImage()
{
  // Get the input image
  ImagePointer input = m_ImageStack.back();
  
  // Simply make a copy of the input image on the stack
  ImagePointer output = ImageType::New();
  output->SetRegions(input->GetBufferedRegion());
  output->SetSpacing(input->GetSpacing());
  output->SetOrigin(input->GetOrigin());
  output->SetMetaDataDictionary(input->GetMetaDataDictionary());
  output->Allocate();

  // Copy everything
  size_t n = input->GetBufferedRegion().GetNumberOfPixels();
  for(size_t i = 0; i < n; i++)
    output->GetBufferPointer()[i] = input->GetBufferPointer()[i];

  // Put on the end of the stack
  m_ImageStack.pop_back();
  m_ImageStack.push_back(output);
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::ThresholdImage(double u1, double u2, double vIn, double vOut)
{
  // Get the input image
  ImagePointer input = m_ImageStack.back();
  
  // Say what we are doing
  *verbose << "Thresholding #" << m_ImageStack.size() << endl;
  *verbose << "  Mapping range [" << u1 << ", " << u2 << "] to " << vIn << endl;
  *verbose << "  Values outside are mapped to " << vOut << endl;

  // Do the thresholding
  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(input);
  filter->SetLowerThreshold(u1);
  filter->SetUpperThreshold(u2);
  filter->SetInsideValue(vIn);
  filter->SetOutsideValue(vOut);
  filter->Update();

  m_ImageStack.pop_back();
  m_ImageStack.push_back(filter->GetOutput());
}

template<unsigned int VDim>
void ExpandRegion(itk::ImageRegion<VDim> &region, const itk::Index<VDim> &idx)
{
  if(region.GetNumberOfPixels() == 0)
    {
    region.SetIndex(idx);
    for(size_t i = 0; i < VDim; i++)
      region.SetSize(i, 1);
    }
  else {
    for(size_t i = 0; i < VDim; i++)
      {
      if(region.GetIndex(i) > idx[i])
        {
        region.SetSize(i, region.GetSize(i) + (region.GetIndex(i) - idx[i]));
        region.SetIndex(i, idx[i]);
        }
      else if(region.GetIndex(i) + region.GetSize(i) <= idx[i])
        {
        region.SetSize(i, 1 + idx[i] - region.GetIndex(i));
        }
      }
  }
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::TrimImage(const RealVector &margin)
{
  // Get the input image
  ImagePointer input = m_ImageStack.back();

  // Debugging info
  *verbose << "Trimming #" << m_ImageStack.size() << endl;
  *verbose << "  Wrapping non-background voxels with margin of " << margin << " voxels." << endl;

  // Initialize the bounding box
  typename ImageType::RegionType bbox;

  // Find the extent of the non-background region of the image
  itk::ImageRegionConstIteratorWithIndex<ImageType> it(input, input->GetBufferedRegion());
  for( ; !it.IsAtEnd(); ++it)
    if(it.Value() != m_Background)
      ExpandRegion(bbox, it.GetIndex());

  // Pad the region by radius specified by user
  typename ImageType::SizeType radius;
  for(size_t i = 0; i < VDim; i++)
    radius[i] = (int) ceil(margin[i] / input->GetSpacing()[i]);
  bbox.PadByRadius(radius);

  // Make sure the bounding box is within the contents of the image
  bbox.Crop(input->GetBufferedRegion());

  // Report the bounding box size
  *verbose << "  Extracting bounding box " << bbox.GetIndex() << " " << bbox.GetSize() << endl;

  // Chop off the region
  typedef itk::ExtractImageFilter<ImageType, ImageType> TrimFilter;
  typename TrimFilter::Pointer fltTrim = TrimFilter::New();
  fltTrim->SetInput(input);
  fltTrim->SetExtractionRegion(bbox);
  fltTrim->Update();

  // Unfortunately we have to adjust the index and origin of the trimmed region
  RealVector bb0, bb1;
  GetBoundingBox(fltTrim->GetOutput(), bb0, bb1);
  typename ImageType::RegionType outreg;
  outreg.SetSize(fltTrim->GetOutput()->GetBufferedRegion().GetSize());
  fltTrim->GetOutput()->SetRegions(outreg);
  fltTrim->GetOutput()->SetOrigin(bb0.data_block());

  // Update the image stack
  m_ImageStack.pop_back();
  m_ImageStack.push_back(fltTrim->GetOutput());
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::ScaleShiftImage(double a, double b)
{
  // Get the input image
  ImagePointer input = m_ImageStack.back();
  
  // Say what we are doing
  *verbose << "Scaling #" << m_ImageStack.size() << " by " << a << " and adding " << b << endl;

  // If a=0, this means setting the image to a constant
  if(a == 0.0)
    {
    CopyImage();
    m_ImageStack.back()->FillBuffer(b);
    return;
    }

  // Create and run filter
  typedef itk::ShiftScaleImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(input);
  filter->SetScale(a);
  filter->SetShift(b / a);
  filter->Update();

  m_ImageStack.pop_back();
  m_ImageStack.push_back(filter->GetOutput());
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::AntiAliasImage(double xIsoSurface)
{
  // Get the input image
  ImagePointer input = m_ImageStack.back();

  // Report what the filter is doing
  *verbose << "Anti-aliasing #" << m_ImageStack.size() << endl;
  *verbose << "  Root Mean Square error: " << m_AntiAliasRMS << endl;
  *verbose << "  Iterations: "; 
  if(m_Iterations == 0) 
    *verbose << "Unlimited" << endl; 
  else 
    *verbose << m_Iterations << endl;

  // Apply antialiasing to the image
  typedef itk::AntiAliasBinaryImageFilter<ImageType,ImageType> AntiFilterType;
  typename AntiFilterType::Pointer fltAnti = AntiFilterType::New();
  fltAnti->SetInput(input);
  fltAnti->SetMaximumRMSError(m_AntiAliasRMS);
  if(m_Iterations > 0)
    fltAnti->SetNumberOfIterations(m_Iterations);
  fltAnti->SetIsoSurfaceValue(xIsoSurface);
  // fltAnti->AddObserver(itk::ProgressEvent(),command);
  fltAnti->Update();

  m_ImageStack.pop_back();
  m_ImageStack.push_back(fltAnti->GetOutput());
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::GetBoundingBox(ImageType *image, RealVector &bb0, RealVector &bb1)
{
  for(size_t i = 0; i < VDim; i++)
    {
    bb0[i] = image->GetOrigin()[i];
    bb1[i] = bb0[i] + image->GetSpacing()[i] * image->GetBufferedRegion().GetSize()[i];
    }
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::PrintImageInfo(bool full)
{
  // Get the input image
  ImagePointer image = m_ImageStack.back();

  // Print the image number
  cout << "Image #" << m_ImageStack.size() << ":";

  // Compute the bounding box
  RealVector bb0, bb1, ospm;
  for(size_t i = 0; i < VDim; i++)
    {
    bb0[i] = image->GetOrigin()[i];
    bb1[i] = bb0[i] + image->GetSpacing()[i] * image->GetBufferedRegion().GetSize()[i];
    ospm[i] = -image->GetOrigin()[i] / image->GetSpacing()[i];
    }

  // Compute the intensity range of the image
  size_t n = image->GetBufferedRegion().GetNumberOfPixels();
  double *vox = image->GetBufferPointer();
  double iMax = vox[0], iMin = vox[0], iMean = vox[0];
  for(size_t i = 1; i < n; i++)
    {
    iMax = (iMax > vox[i]) ? iMax : vox[i];
    iMin = (iMin < vox[i]) ? iMin : vox[i];
    iMean += vox[i];
    }
  iMean /= n;

  // Short or long?
  if(!full) 
    {
    cout << " dim = " << image->GetBufferedRegion().GetSize() << "; ";
    cout << " bb = {[" << bb0 << "], [" << bb1 << "]}; ";
    cout << " vox = " << image->GetSpacing() << "; ";
    cout << " range = [" << iMin << ", " << iMax << "]; ";
    cout << endl;
    }
  else
    {
    cout << endl;
    cout << "  Image Dimensions: " << image->GetBufferedRegion().GetSize() << endl;
    cout << "  Bounding Box: " << "{[" << bb0 << "], [" << bb1 << "]}" << endl;
    cout << "  BufferedRegion: " << image->GetBufferedRegion() << endl;
    cout << "  LargestPossibleRegion: " << image->GetLargestPossibleRegion() << endl;
    cout << "  Voxel Size: " << image->GetSpacing() << endl;
    cout << "  Direction Matrix: " << image->GetDirection() << endl;
    cout << "  Position of Origin in Voxel Units (SPM Origin): " << ospm << endl;
    cout << "  Intensity Range: [" << iMin << ", " << iMax << "]" << endl;
    cout << "  Mean Intensity: " << iMean << endl;
    cout << "  Image Metadata: " << endl;

    // Print metadata
    itk::MetaDataDictionary &mdd = image->GetMetaDataDictionary();
    itk::MetaDataDictionary::ConstIterator itMeta;
    for(itMeta = mdd.Begin(); itMeta != mdd.End(); ++itMeta)
      {
      // Get the metadata as a generic object
      string key = itMeta->first, value;
      if(itk::ExposeMetaData<string>(mdd, key, value))
        {
        // For some weird reason, some of the strings returned by this method
        // contain '\0' characters. We will replace them by spaces
        std::ostringstream sout("");
        for(unsigned int i=0;i<value.length();i++)
          if(value[i] >= ' ') sout << value[i];
        value = sout.str();

        // Make sure the value has more than blanks
        if(value.find_first_not_of(" ") != value.npos)
          cout << "    " << key << " = " << value << endl;
        }
      }
    }
}


  
int main(int argc, char *argv[])
{
  if(argc == 1)
    {
    cerr << "PICSL convert3d tool" << endl;
    cerr << "usage: http://milesdavis.uphs.upenn.edu/mediawiki-1.4.10/index.php?title=Convert3D_Command_Line" << endl;
    return -1;
    }

  // Load the ITK factories
  itk::ObjectFactoryBase::RegisterFactory(itk::VoxBoCUBImageIOFactory::New());
  itk::ObjectFactoryBase::RegisterFactory(itk::PovRayDF3ImageIOFactory::New());


  // Create the converter and run it
  try {
    ImageConverter<double, 3> convert;
    convert.ProcessCommandLine(argc, argv);
    return 0;
  }
  catch (...) { return -1; }
}

