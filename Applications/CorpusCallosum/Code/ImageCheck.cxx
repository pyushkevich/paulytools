#include <iostream>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkLinearInterpolateImageFunction.h"
#include "optima.h"

using namespace std;

typedef float PixelType;
typedef itk::Image< PixelType, 2 >	ImageType;
typedef itk::ImageFileReader< ImageType >	ImageReader;
typedef itk::LinearInterpolateImageFunction< ImageType, double > ImageInterpolator;

int main(int argc, char *argv[]) {
	if (argc < 2) {
		cerr << "Not enough arguments" << endl;
		cerr << "Usage: " << argv[0] << " inputImageFile" << endl;
		exit (1);
	}

  const char *inputImageFile = argv[1];
	
  ImageReader::Pointer reader = ImageReader::New();
  reader->SetFileName(inputImageFile);
  reader->Update();
  ImageInterpolator::Pointer im = ImageInterpolator::New();
  im->SetInputImage(reader->GetOutput());
  
  double sumPixelValue = 0.0;
  double sumNoneZeroPixel = 0.0;
  double lowerBound = 10000.0;
  double higherBound = -10000.0;
  for (int i = 30; i < 220; ++i) {
    for (int j = 30; j < 220; ++j) {
      double index[2] = {i,j};
      double tmp = im->Evaluate(index);
      sumPixelValue += im->Evaluate(index);
      if(tmp > 0.5) {
	sumNoneZeroPixel += 1.0;
      }
      if (lowerBound > tmp) {lowerBound = tmp;}
      if (higherBound < tmp) {higherBound = tmp;}   
    }
  }

  int counterlowerBound = 0;
  int counterhigherBound = 0;
  for (int i = 0; i < 256; ++i) {
    for (int j = 0; j < 256; ++j) {
      double index[2] = {i,j};
      double tmp = im->Evaluate(index);
      if (tmp == lowerBound) { counterlowerBound++;}
      if (tmp == higherBound) { counterhigherBound++;}
    }
  }


  double index1[2] = {128, 128};
  double index2[2] = {10,10};
  cout << "Sum of Pixel Value = " << sumPixelValue << endl;
  cout << " Sum of none zero pixel numbers = " << sumNoneZeroPixel << endl;
  cout << "lower Bound of image pixel value = " << lowerBound << endl;
  cout << "higher Bound of image pixel value = " << higherBound << endl;
  cout << "pixel value at (0 0) = " << im->Evaluate(index1) << endl;
  cout << "pixel value at (-128 -128) = " << im->Evaluate(index2) << endl;
  cout << "lower bound pixel number = " << counterlowerBound << endl;
  cout << "higher bound pixel number = " << counterhigherBound << endl;

  Vector ve;
  Vector ye;
  ve.setSize(0);
  ye = ve;
  int size = ye.size();
  cout << "size of vector = " << size << endl;
  return 0;

}
