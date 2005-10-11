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
  cout << "    randfx [options] in1.img in2.img in3.img output.img " << endl;
  cout << "  options; " << endl;
  cout << "    none " << endl;
  
  return -1;
}

int main(int argc, char *argv[])
{
  if(argc < 4) return usage();

  // Get the number of input images
  size_t n = argc - 2;

  // Typedefs
  typedef itk::Image<float, 3> ImageType;
  typedef ImageType::Pointer ImagePointer;
  typedef itk::ImageFileReader<ImageType> ImageReader;
  
  // Read the input images
  vector<ImagePointer> imgInput;
  for(size_t i = 0; i < n; i++)
    {
    ImageReader::Pointer fltReader = ImageReader::New();
    fltReader->SetFileName(argv[i+1]);
    fltReader->Update();
    imgInput.push_back(fltReader->GetOutput());
    }

  // Allocate the output image
  ImagePointer imgOutput = ImageType::New();
  imgOutput->SetRegions(imgInput[0]->GetBufferedRegion());
  imgOutput->Allocate();
  imgOutput->SetSpacing(imgInput[0]->GetSpacing());
  imgOutput->SetOrigin(imgInput[0]->GetOrigin());

  // Just compute the t-test on the pixel arrays
  size_t m = imgOutput->GetBufferedRegion().GetNumberOfPixels();
  double sqrt_n = sqrt(1.0 * n);
  for(size_t j = 0; j < m; j++)
    {
    // Compute the sum and sum of squares of the values
    double sx = 0.0, sx2 = 0.0;
    for(size_t i = 0; i < n; i++)
      {
      float x = imgInput[i]->GetBufferPointer()[j];
      sx += x; sx2 += x * x;
      }

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
  fltWriter->SetFileName(argv[argc-1]);
  fltWriter->SetInput(imgOutput);
  fltWriter->Update();
}
  
  

