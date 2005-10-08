/** Stackem - a program to stack together (and optionally coregister) a
 * bunch of 2D image files */

#include <iostream>
#include <sstream>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkPNGImageIO.h"
#include "itkJPEGImageIO.h"
#include "itkTIFFImageIO.h"
#include "itkPNGImageIO.h"
#include "itksys/Directory.hxx"
#include "itksys/SystemTools.hxx"

#include "itkImageRegistrationMethod.h"
#include "itkTranslationTransform.h"
#include "itkAffineTransform.h"
#include "itkAffineTransform.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkGradientDescentOptimizer.h"

using namespace std;
using namespace itksys;

int usage(const char *message = NULL)
{
  if(message)
    cerr << message << endl << endl;
  
  cout << "Stack'em: stack a list of 2D images into an Analyze 3D image." << endl;
  cout << "USAGE:" << endl;
  cout << "   stackem [options] input_dir output_file.hdr" << endl;
  cout << "REQUIRED PARAMETERS:" << endl;
  cout << "   input_dir    A directory containing 2D slices as PNG, JPG or TIFF files" << endl;
  cout << "   output_file  Name of the output file" << endl;
  cout << "OPTIONS:" << endl;
  // cout << "   -q           Quiet mode: suppress info messages" << endl;
  // cout << "   -r           Attempt to co-register slices" << endl;
  cout << "   -x x.xx      Pixel spacing (mm) in x direction [default: 1.0]" << endl;
  cout << "   -y x.xx      Pixel spacing (mm) in y direction [default: 1.0]" << endl;
  cout << "   -z x.xx      Spacing (mm) between consecutuve slices [default: 1.0]" << endl;
  cout << "EXAMPLE:" << endl;
  cout << "   stackem -z 1.5 ../slices out.hdr" << endl;
  return -1;
}

typedef itk::Image<short,2> SliceType;

void registerSlices(SliceType *img1, SliceType *img2)
  {
  typedef SliceType FixedImageType;
  typedef SliceType MovingImageType;
  
  // Transform Type
  typedef itk::TranslationTransform< double, 2 > TransformType;
  typedef TransformType::ParametersType             ParametersType;

  // Optimizer Type
  typedef itk::GradientDescentOptimizer                  OptimizerType;

  // Metric Type
  typedef itk::MeanSquaresImageToImageMetric< 
                                    FixedImageType, 
                                    MovingImageType >    MetricType;

  // Interpolation technique
  typedef itk:: LinearInterpolateImageFunction< 
                                    MovingImageType,
                                    double >             InterpolatorType;

  // Registration Method
  typedef itk::ImageRegistrationMethod< 
                                    FixedImageType, 
                                    MovingImageType >    RegistrationType;

  //typedef itk::CommandIterationUpdate<  
  //                                OptimizerType >    CommandIterationType;


  MetricType::Pointer         metric        = MetricType::New();
  TransformType::Pointer      transform     = TransformType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  TransformType::Pointer      trasform      = TransformType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  FixedImageType::ConstPointer     fixedImage    = img1;
  MovingImageType::ConstPointer    movingImage   = img2;

  //
  // Connect all the components required for Registratio
  //
  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );
  registration->SetTransform(     transform     );
  registration->SetFixedImage(    fixedImage    );
  registration->SetMovingImage(   movingImage   );
  registration->SetInterpolator(  interpolator  );

  // Select the Region of Interest over which the Metric will be computed
  // Registration time will be proportional to the number of pixels in this region.
  metric->SetFixedImageRegion( fixedImage->GetBufferedRegion() );

  // Instantiate an Observer to report the progress of the Optimization
  //CommandIterationType::Pointer iterationCommand = CommandIterationType::New();
  //iterationCommand->SetOptimizer(  optimizer.GetPointer() );

  // Scale the translation components of the Transform in the Optimizer
  OptimizerType::ScalesType scales( transform->GetNumberOfParameters() );
  scales.Fill( 1.0 );
  
  unsigned long   numberOfIterations =  100;
  double          translationScale   = 1e-6;
  double          learningRate       = 1e-8;

  for( unsigned int i=0; i<2; i++)
    {
    scales[ i + 2 * 2 ] = translationScale;
    }

  optimizer->SetScales( scales );
  optimizer->SetLearningRate( learningRate );
  optimizer->SetNumberOfIterations( numberOfIterations );
  optimizer->MinimizeOn();

  // Start from an Identity transform (in a normal case, the user 
  // can probably provide a better guess than the identity...
  transform->SetIdentity();
  registration->SetInitialTransformParameters( transform->GetParameters() );

  // Initialize the internal connections of the registration method. 
  // This can potentially throw an exception
  try
    {
    registration->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << e << std::endl;
    }

  ParametersType finalParameters  = registration->GetLastTransformParameters();

  //
  //  Get the transform as the Output of the Registration filter
  //
  RegistrationType::TransformOutputConstPointer transformDecorator = 
                                                        registration->GetOutput();

  TransformType::ConstPointer finalTransform = 
    static_cast< const TransformType * >( transformDecorator->Get() ); 

  cout << "Trasform computed:" << finalParameters << endl;
  }

int main(int argc, char *argv[])
{
  // Type definitions
  typedef itk::Image<short, 3> ImageType;
  typedef ImageType::SpacingType SpacingType;
  typedef list<string> FileList;
  typedef FileList::iterator FileIterator; 
  
  // Parameters to the program
  const char *fnOutput = NULL;
  const char *fnDirectory = ".";
  SpacingType xSpacing(1.0);
  bool flagCoregister = false;
  bool flagVerbose = true;
  
  // Check the number of options
  if(argc < 3) return usage("Incorrect number of parameters!");

  // Get the required options
  fnDirectory = argv[argc - 2];
  fnOutput = argv[argc - 1];

  // Read in the required and optional parameters
  for(unsigned int i = 1; i < argc - 2; i++)
    {
    if(!strcmp(argv[i], "-x")) 
      xSpacing[0] = atof(argv[++i]);
    else if(!strcmp(argv[i], "-y")) 
      xSpacing[1] = atof(argv[++i]);
    else if(!strcmp(argv[i], "-z")) 
      xSpacing[2] = atof(argv[++i]);
    else if(!strcmp(argv[i],"-q"))
      flagVerbose = false;
    else if(!strcmp(argv[i],"-r"))
      flagCoregister = true;
    else
      return usage("Unknown option encountered!");
    }

  // Open the input directory 
  Directory dirInput;
  if(!dirInput.Load(fnDirectory))
    return usage("Could not open the input directory!");

  // Create a PNG IO object to check readability
  itk::PNGImageIO::Pointer ioPNG  = itk::PNGImageIO::New();
  itk::TIFFImageIO::Pointer ioTIFF = itk::TIFFImageIO::New();
  itk::JPEGImageIO::Pointer ioJPEG = itk::JPEGImageIO::New();
  
  // Parse the directory contents, checking files
  FileList fnInput;
  for(unsigned long iFile = 0; iFile < dirInput.GetNumberOfFiles(); iFile++)
    {
    string fnCandidate = string(fnDirectory) + "/" + dirInput.GetFile(iFile);
    if(ioPNG->CanReadFile(fnCandidate.c_str())
      || ioTIFF->CanReadFile(fnCandidate.c_str())
      || ioJPEG->CanReadFile(fnCandidate.c_str())) fnInput.push_back(fnCandidate);
    }

  // Check that there are input images
  if(fnInput.size() < 2) return usage("The directory must contain at least 2 input files!");

  // Check that there is an output image
  if(!fnOutput) return usage("Please specify the output image!");

  // Set up the verbose output
  ostream &vout = flagVerbose ? cout : clog;

  // Print current options
  vout << "Stack'em: read command line parameters" << endl;
  vout << "   input images :   " << fnInput.size() << endl;
  vout << "      first :       " << fnInput.front() << endl;
  vout << "      last :        " << fnInput.back() << endl;
  vout << "   output image :   " << fnOutput << endl;
  vout << "   output spacing : " << xSpacing << endl ;
  vout << "   coregistration : " << (flagCoregister ? "yes" : "no") << endl;
  vout << endl;

  // A region object for keeping track of the size of the images
  typedef itk::Image<short,2> SliceType;
  typedef itk::ImageFileReader<SliceType> SliceReader;
  typedef SliceType::RegionType SliceRegion;
  SliceRegion rgnCommon;

  // A list of loaded images
  vector<SliceType::Pointer> imgSlices;
  
  // Load each of the images
  FileIterator itFile = fnInput.begin();
  while(itFile != fnInput.end())
    {
    try 
      {
      // Read the image
      SliceReader::Pointer fltReader = SliceReader::New();
      fltReader->SetFileName(itFile->c_str());
      fltReader->Update();

      // Get the image extent
      SliceType::Pointer imgSlice = fltReader->GetOutput();
      SliceRegion rgn = imgSlice->GetBufferedRegion();
      for(unsigned int d = 0; d < 2; d++)
        if(rgnCommon.GetSize(d) < rgn.GetSize(d))
          rgnCommon.SetSize(d, rgn.GetSize(d));

      // Add the image to the list of stored images
      imgSlices.push_back(imgSlice);
      }
    catch(...)
      {
      cerr << "Error reading file " << *itFile << "; continuing... " << endl;
      }   
    ++itFile;
    }

  // Report on the common region
  vout << "Stack'em: loaded " << imgSlices.size() << " slice images" << endl;
  vout << "   common slice size: " << rgnCommon.GetSize() << endl; 

  if(flagCoregister)
    {
    for(unsigned int iSlice = 1; iSlice < imgSlices.size(); iSlice++)
      {
      registerSlices(imgSlices[iSlice-1],imgSlices[iSlice]);
      }
    }

  // If no registration is required, simply combine the slices and save the image
  if(!flagCoregister)
    {
    // Create the 3D region
    ImageType::RegionType rgnVolume;
    rgnVolume.SetSize(0, rgnCommon.GetSize(0));
    rgnVolume.SetSize(1, rgnCommon.GetSize(1));
    rgnVolume.SetSize(2, (unsigned long) imgSlices.size() );

    // Create an output volume to store stacked images
    ImageType::Pointer imgOutput = ImageType::New();
    imgOutput->SetRegions(rgnVolume);
    imgOutput->SetSpacing(xSpacing);
    imgOutput->Allocate();
    imgOutput->FillBuffer(0);
    
    // Create an indexable iterator into the output image
    typedef itk::ImageRegionIteratorWithIndex<ImageType> OutputIterator;
    OutputIterator itVolume(imgOutput, rgnVolume);
    
    // An index for the iterator
    OutputIterator::IndexType idxVolume;
    
    // Progress reporting
    vout << "Stack'em: stacking slices into the volume" << endl;
    
    // Read each slice (slow but so what)
    for(unsigned int iSlice = 0; iSlice < imgSlices.size(); iSlice++)
      {
      // Create an iterator for the slice
      typedef itk::ImageRegionConstIteratorWithIndex<SliceType> SliceIterator;
      SliceIterator itSlice(
        imgSlices[iSlice], imgSlices[iSlice]->GetBufferedRegion());

      // Set the z-coordinate in volume iterator
      idxVolume[2] = iSlice;
      
      // Copy the pixels
      while(!itSlice.IsAtEnd())
        {
        // Update the index
        idxVolume[0] = itSlice.GetIndex()[0];
        idxVolume[1] = itSlice.GetIndex()[1];

        // Set the pixel
        itVolume.SetIndex(idxVolume);
        itVolume.Set(itSlice.Value());
        
        // Go to next pixel
        ++itSlice;
        }

      // Report progress
      vout << "   added slice " << iSlice << endl;
      }

    // Save the 3D image
    typedef itk::ImageFileWriter<ImageType> ImageWriter;
    ImageWriter::Pointer fltWriter = ImageWriter::New();
    fltWriter->SetFileName(fnOutput);
    fltWriter->SetInput(imgOutput);

    try 
      {
      fltWriter->Update();
      vout << "Stack'em: saved output image " << fnOutput << endl; 
      }
    catch(...)
      {
      cerr << "Error writing output image " << fnOutput << endl;
      }
    }
}


