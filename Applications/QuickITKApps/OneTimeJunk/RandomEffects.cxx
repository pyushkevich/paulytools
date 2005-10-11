#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <iostream>

using namespace std;
using namespace itk;

int usage()
{
  cout << "  randfx: Computes random effects model given N images " << endl;
  cout << "  usage: " << endl;
  cout << "    randfx [options] output.img" << endl;
  cout << "  options; " << endl;
  cout << "    -i img1 .. imgN    Input image list" << endl;
  cout << "    -m img1 .. imgN    Mask image list " << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  if(argc < 4) return usage();

  // Process the command line options
  string fnOutput = argv[argc-1];
  vector<string> fnInput, fnMasks;
  string sLastOption = "";
  for(size_t k = 1; k < argc-1; k++)
    {
    if(argv[k][0] == '-')
      if(argv[k][1] == 'i')
        sLastOption = "i";
      else if(argv[k][1] == 'm')
        sLastOption = "m";
      else
        return usage();
    else if(sLastOption == "i")
      fnInput.push_back(argv[k]);
    else if(sLastOption == "m")
      fnMasks.push_back(argv[k]);
    else
      return usage();
    }

  // Make sure the number of masks matches the images
  size_t n = fnInput.size();
  bool flagMasks = fnMasks.size() > 0;
  if(flagMasks && fnMasks.size() != n)
    {
    cout << "Number of images and masks do not match!" << endl;
    return usage();
    }

  // Print report
  cout << "Performing analysis on " << fnInput.size() << " images, " 
    << fnMasks.size() << " masks; writing to " << fnOutput << endl;

  // Typedefs
  typedef itk::Image<float, 3> ImageType;
  typedef ImageType::Pointer ImagePointer;
  typedef itk::ImageFileReader<ImageType> ImageReader;
  typedef itk::Image<short, 3> MaskImageType;
  typedef MaskImageType::Pointer MaskImagePointer;
  typedef itk::ImageFileReader<MaskImageType> MaskImageReader;
  
  // Read the input images
  vector<ImagePointer> imgInput;
  vector<MaskImagePointer> imgMask;
  for(size_t i = 0; i < n; i++)
    {
    // Read the input image
    ImageReader::Pointer fltReader = ImageReader::New();
    fltReader->SetFileName(fnInput[i].c_str());
    fltReader->Update();
    imgInput.push_back(fltReader->GetOutput());

    // Read the mask image
    if(flagMasks)
      {
      MaskImageReader::Pointer fltMaskReader = MaskImageReader::New();
      fltMaskReader->SetFileName(fnMasks[i].c_str());
      fltMaskReader->Update();
      imgMask.push_back(fltMaskReader->GetOutput());
      }
    }

  // Allocate the output image
  ImagePointer imgOutput = ImageType::New();
  imgOutput->SetRegions(imgInput[0]->GetBufferedRegion());
  imgOutput->Allocate();
  imgOutput->SetSpacing(imgInput[0]->GetSpacing());
  imgOutput->SetOrigin(imgInput[0]->GetOrigin());
  imgOutput->FillBuffer(0.0f);

  // Just compute the t-test on the pixel arrays
  size_t m = imgOutput->GetBufferedRegion().GetNumberOfPixels();
  double sqrt_n = sqrt(1.0 * n);
  for(size_t j = 0; j < m; j++)
    {
    // Compute the sum and sum of squares of the values
    int N = 0;
    double sx = 0.0, sx2 = 0.0;
    for(size_t i = 0; i < n; i++)
      {
      // Is the pixel masked?
      if(!flagMasks || imgMask[i]->GetBufferPointer()[j] > 0)
        {
        float x = imgInput[i]->GetBufferPointer()[j];
        sx += x; sx2 += x * x; N++;
        }
      }

    // Only proceed if there is full overlap
    if(N < n) continue;

    // Compute the estimate of the mean
    double xbar = sx / n;
    
    // Compute the estimate of variance
    double ev = ( sx2 - n * xbar * xbar ) / (n - 1);

    // Compute the t statistic
    double t = xbar / (ev / sqrt_n);
    
    // Set the output intensity to be that
    imgOutput->GetBufferPointer()[j] = (float) t;
    }
  
  // Save the t-test image
  typedef ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer fltWriter = WriterType::New();
  fltWriter->SetFileName(fnOutput.c_str());
  fltWriter->SetInput(imgOutput);
  fltWriter->Update();
}
  
  

