#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkVoxBoCUBImageIOFactory.h"

#include <vnl/vnl_random.h>

#include <vector>
#include <algorithm>
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
  cout << "    -p N               Permutation test (N trials)" << endl;
  cout << "    -hb                Compute a mask of hotelling T2 along" << endl;
  cout << "                       the cm-rep boundary" << endl;
  cout << "    -mb                Compute a mask of max inward t-stat along" << endl;
  cout << "                       the cm-rep boundary" << endl;
  return -1;
}

// Typedefs
typedef itk::Image<float, 3> ImageType;
typedef ImageType::Pointer ImagePointer;
typedef itk::ImageFileReader<ImageType> ImageReader;
typedef itk::Image<short, 3> MaskImageType;
typedef MaskImageType::Pointer MaskImagePointer;
typedef itk::ImageFileReader<MaskImageType> MaskImageReader;

void ComputePointwiseT(
  vector<ImagePointer> &imgInput, 
  vector<MaskImagePointer> &imgMask, 
  bool flagMasks,
  ImageType *imgOutput)
{
  // Just compute the t-test on the pixel arrays
  size_t m = imgOutput->GetBufferedRegion().GetNumberOfPixels();
  size_t n = imgInput.size();
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
    double t = xbar / sqrt(ev / n);
    
    // Set the output intensity to be that
    imgOutput->GetBufferPointer()[j] = (float) t;
    }
}

double TTest(const vnl_vector<double > &x)
{
  double mu = 0, s = 0;
  for(size_t i = 0; i < x.size(); i++)
    {
    mu += x[i];
    s += x[i] * x[i];
    }
  mu /= x.size();
  s = (s - x.size() * mu * mu) / (x.size() - 1);
  return mu / sqrt(s / x.size());
}

double MaxT(vnl_matrix<double> &A)
{
  double tMax = 0.0;
  for(size_t i = 0; i < A.columns(); i++)
    {
    double t = TTest(A.get_column(i));
    tMax = (tMax < t) ? t : tMax;
    }
  return tMax;
}

double HotellingT2(vnl_matrix<double> &A)
{
  size_t n = A.rows(), m = A.columns(), i, j;
  
  // Compute the mean
  vnl_vector<double> mu(m, 0.0);
  for(i = 0; i < n; i++)
    mu += A.get_row(i);
  mu /= 1.0 * n;

  // Subtract the mean from A
  vnl_matrix<double> B(n, m, 0.0);
  for(i = 0; i < n; i++)
    B.set_row(i, A.get_row(i) - mu);

  // Compute the sample covariance matrix
  vnl_matrix<double> S = (B.transpose() * B) / (n - 1);

  // Compute the statistic
  vnl_matrix_inverse<double> iS(S);
  double T2 = dot_product(mu, iS * mu) * n;

  return T2;
}

void ComputeHotellingBoundaryStat(
  vector<ImagePointer> &imgInput, 
  ImageType *imgOutput)
{
  size_t nx = imgOutput->GetBufferedRegion().GetSize()[0];
  size_t ny = imgOutput->GetBufferedRegion().GetSize()[1];
  size_t nz = imgOutput->GetBufferedRegion().GetSize()[2];
  size_t N = imgInput.size();
  size_t k = 1 + nz / 2;

  // All computations repeated over x, y
  size_t x, y, z, i, j;
  for(x = 0; x < nx; x++) 
    {
    for(y = 0; y < ny; y++)
      {
      vnl_matrix<double> A1(N, k), A2(N, k);
      itk::Index<3> idx; idx[0] = x; idx[1] = y;

      // Compute the data matrix
      for(i = 0; i < N; i++) for(j = 0; j < k; j++)
        {
        idx[2] = j; A1(i,j) = imgInput[i]->GetPixel(idx);
        idx[2] = nz - (1 + j); A2(i,j) = imgInput[i]->GetPixel(idx);
        }

      // Compute the sample covariance matrix
      double t1 = HotellingT2(A1);
      double t2 = HotellingT2(A2);

      // Fill the output image 
      for(j = 0; j < k; j++) 
        {
        idx[2] = j; imgOutput->GetPixel(idx) = t1;
        idx[2] = nz - (1 + j); imgOutput->GetPixel(idx) = t2;
        }
      }
    }
}

void ComputeMaximumTStat(
  vector<ImagePointer> &imgInput, 
  ImageType *imgOutput)
{
  size_t nx = imgOutput->GetBufferedRegion().GetSize()[0];
  size_t ny = imgOutput->GetBufferedRegion().GetSize()[1];
  size_t nz = imgOutput->GetBufferedRegion().GetSize()[2];
  size_t N = imgInput.size();
  size_t k = 1 + nz / 2;

  // All computations repeated over x, y
  size_t x, y, z, i, j;
  for(x = 0; x < nx; x++) 
    {
    for(y = 0; y < ny; y++)
      {
      vnl_matrix<double> A1(N, k), A2(N, k);
      itk::Index<3> idx; idx[0] = x; idx[1] = y;

      // Compute the data matrix
      for(i = 0; i < N; i++) for(j = 0; j < k; j++)
        {
        idx[2] = j; A1(i,j) = imgInput[i]->GetPixel(idx);
        idx[2] = nz - (1 + j); A2(i,j) = imgInput[i]->GetPixel(idx);
        }

      // Compute the sample covariance matrix
      double t1 = MaxT(A1);
      double t2 = MaxT(A2);

      // Fill the output image 
      for(j = 0; j < k; j++) 
        {
        idx[2] = j; imgOutput->GetPixel(idx) = t1;
        idx[2] = nz - (1 + j); imgOutput->GetPixel(idx) = t2;
        }
      }
    }
}


int main(int argc, char *argv[])
{
  // Enable support for VoxBo
  itk::ObjectFactoryBase::RegisterFactory(itk::VoxBoCUBImageIOFactory::New());
  
  if(argc < 4) return usage();

  // Process the command line options
  string fnOutput = argv[argc-1];
  vector<string> fnInput, fnMasks;
  string sLastOption = "";
  bool doHotellingBoundary = false;
  bool doMaxTOnBoundary = false;
  size_t nPermutations = 0;

  for(size_t k = 1; k < argc-1; k++)
    {
    if(argv[k][0] == '-')
      if(strstr(argv[k],"-i") == 0)
        sLastOption = "i";
      else if(strstr(argv[k],"-m") == 0)
        sLastOption = "m";
      else if(strcmp(argv[k],"-hb") == 0)
        {
        doHotellingBoundary = true;
        sLastOption = "x";
        }
      else if(strcmp(argv[k],"-mb") == 0)
        {
        doMaxTOnBoundary = true;
        sLastOption = "x";
        }
      else if(strcmp(argv[k],"-p") == 0)
        {
        nPermutations = atoi(argv[++k]);
        }
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

  // Read the input images
  vector<ImagePointer> imgInput;
  vector<MaskImagePointer> imgMask;
  for(size_t i = 0; i < n; i++)
    {
    // Read the input image
    try 
      {
      ImageReader::Pointer fltReader = ImageReader::New();
      fltReader->SetFileName(fnInput[i].c_str());
      fltReader->Update();
      imgInput.push_back(fltReader->GetOutput());
      }
    catch (itk::ExceptionObject &exc)
      {
      cerr << "Exception loading image " << fnInput[i] << endl;
      cerr << exc << endl;
      return -1;
      }

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

  // Compute the hotelling statistics if needed
  if(doHotellingBoundary)
    ComputeHotellingBoundaryStat(imgInput, imgOutput);
  else if(doMaxTOnBoundary)
    ComputeMaximumTStat(imgInput, imgOutput);
  else
    ComputePointwiseT(imgInput, imgMask, flagMasks, imgOutput);

  // Has the user requested permutation testing?
  if(nPermutations > 0)
    {
    // Create a random number generator
    vnl_random randy;
    
    // Allocate the permutation statistic image
    ImagePointer imgTrial = ImageType::New();
    imgTrial->SetRegions(imgInput[0]->GetBufferedRegion());
    imgTrial->Allocate();
    imgTrial->SetSpacing(imgInput[0]->GetSpacing());
    imgTrial->SetOrigin(imgInput[0]->GetOrigin());
    imgTrial->FillBuffer(0.0f);

    // Create a histogram
    vector<float> xHistogram;

    size_t n = imgTrial->GetBufferedRegion().GetNumberOfPixels();
    for(size_t i = 0; i < nPermutations; i++)
      {
      // For each permutation, we need to flip signs of all input data
      for(size_t j = 0; j < imgInput.size(); j++)
        {
        short q = 0; unsigned long r32 = randy.lrand32();
        float *pix = imgInput[j]->GetBufferPointer();
        for(size_t k = 0; k < n; k++)
          {
          // Flip the sign at random
          if((r32 & 0x01) == 1)
            pix[k] = -pix[k];

          // Update the random number source
          q = (q + 1) % 32;
          r32 = (q == 0) ? randy.lrand32() : r32 >> 1;
          }
        }

      // Now, compute the statistics image
      if(doHotellingBoundary)
        ComputeHotellingBoundaryStat(imgInput, imgTrial);
      else if(doMaxTOnBoundary)
        ComputeMaximumTStat(imgInput, imgOutput);
      else
        ComputePointwiseT(imgInput, imgMask, flagMasks, imgTrial);

      // Compute the maximum of the statistical field
      float xMax = 0.0; float *pix = imgTrial->GetBufferPointer();
      for(size_t k = 0; k < n; k++)
        if(xMax < pix[k]) xMax = pix[k];

      // Stick this maximum into a histogram
      xHistogram.push_back(xMax);

      // Indicate progress
      cout << "." << " " << xMax << " " << flush;
      }

    // Sort the histogram
    sort(xHistogram.begin(), xHistogram.end());
    size_t iPos = (size_t) (0.95 * xHistogram.size());

    // Report the 95th percentile from the histogram
    cout << endl;
    cout << "95th Percentile is " << xHistogram[iPos] << endl;
    }
  
  // Save the t-test image
  typedef ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer fltWriter = WriterType::New();
  fltWriter->SetFileName(fnOutput.c_str());
  fltWriter->SetInput(imgOutput);
  fltWriter->Update();
}
  
  

