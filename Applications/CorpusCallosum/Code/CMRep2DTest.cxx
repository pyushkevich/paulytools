#include <iostream>
#include <cstdlib>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkDiscreteGaussianImageFilter.h"

#include "CMRep2D.h"
#include "TemplatedWaveletRep.h"

using namespace std;

int main (int argc, char *argv[]) {
	typedef float PixelType;
	const int		Dimension = 2;
	typedef itk::Image< PixelType, Dimension >	ImageType;
	typedef itk::ImageFileReader< ImageType >	ImageReader;
	typedef itk::ImageFileWriter< ImageType >	ImageWriter;
	typedef itk::LinearInterpolateImageFunction< ImageType, double > ImageInterpolator;
	typedef itk::DiscreteGaussianImageFilter< ImageType, ImageType > ImageSmoother;
	
	if (argc < 3) {
		cerr << "Not enough arguments" << endl;
		cerr << "Usage: " << argv[0] << " inputImageFile outputImageFile" << endl;
//		cerr << "Usage: " << argv[0] << " inputImageFile dim" << endl;
		exit (1);
	}
	
	const char *inputImageFile = argv[1];
	const char *outputImageFile = argv[2];
	
	cout << "input Image File is " << inputImageFile << endl;
	cout << "output Image File is " << outputImageFile << endl;
	
	ImageReader::Pointer reader = ImageReader::New();
	reader->SetFileName( inputImageFile );
	reader->Update();
	
	// Gaussian smooth the image
	ImageSmoother::Pointer smooth = ImageSmoother::New();
	smooth->SetVariance(4.0);
	smooth->SetInput(reader->GetOutput());
	smooth->Update();
	
	ImageWriter::Pointer writer = ImageWriter::New();
	writer->SetFileName( outputImageFile );
	writer->SetInput( smooth->GetOutput() );
	writer->Update();
	
/*	
	ImageInterpolator::Pointer interpolator = ImageInterpolator::New();
	interpolator->SetInputImage( reader->GetOutput() );

	double areaOfImage = 0.0;
	for (int i = 6; i < 250; ++i) {
		for (int j = 6; j < 250; ++j) {
			double index[2] = {i,j};
			areaOfImage += interpolator->Evaluate(index);
		}
	}
	
	const double tXCoeff[] = {
	  0.0,
	  0.0, 0.0,
	  0.0, 0.0,
	  0.0, 0.0, 0.0, 0.0,
	  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	};
	
	const double tYCoeff[] = {
	  0.0,
	  0.0, 0.0,
	  0.0, 0.0,
	  0.0, 0.0, 0.0, 0.0,
	  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	};
	
	const double tRhoCoeff[] = {
	  0.0,
	  0.0, 0.0,
	  0.0, 0.0,
	  0.0, 0.0, 0.0, 0.0,
	  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	};
	
	double dXCoeff[] = {
		50.0,
		-100.0, 0.0,
		0.0, 0.0,
		0.0, 0.0, 0.0, 0.0
	};
	
	double dYCoeff[] = {
		-10.0,
		0.0, 100.0,
		0.0, 0.0,
		0.0, 0.0, 0.0, 0.0
	};
	
	double dRhoCoeff[] = {
		-0.025,
		0.0, 0.0,
		0.0, 0.0,
		0.0, 0.0, 0.0, 0.0
	};
	
	// set up fx, fy and frho
	const int tJMax = 3;
	const int dJMax = 2;
	TemplatedWaveletRep fx(tXCoeff, tJMax, dXCoeff, dJMax);
	TemplatedWaveletRep fy(tYCoeff, tJMax, dYCoeff, dJMax);
	TemplatedWaveletRep frho(tRhoCoeff, tJMax, dRhoCoeff, dJMax);
	cout << "fx(0.4) = " << fx.get(0.4) << endl;
	cout << "fy(0.4) = " << fy.get(0.4) << endl;
	cout << "frho(0.4) = " << frho.get(0.4) << endl;
	cout << "fx_1stD(0.4) = " << fx.get1stDeriv(0.4) << endl;
	cout << "fy_1stD(0.4) = " << fy.get1stDeriv(0.4) << endl;
	cout << "frho_1stD(0.4) = " << frho.get1stDeriv(0.4) << endl;
	cout << "fx_2ndD(0.4) = " << fx.get2ndDeriv(0.4) << endl;
	cout << "fy_2ndD(0.4) = " << fy.get2ndDeriv(0.4) << endl;
	cout << "frho_2ndD(0.4) = " << frho.get2ndDeriv(0.4) << endl;

	const int dim = atoi(argv[2]);
	CMRep2D cm(dim);
	
	double* phi = new double [dim + 3];
	for (int i = 0; i < dim + 3; ++i) {
		phi[i] = 1.0;
	}
	
	// solve for the cmrep
	double minPhi = cm.buildCMRep2D(fx, fy, frho, phi);
	cout << "minPhi = " << minPhi << endl;
	// compute the area overlap to white background
	double areaOverlap = cm.computeAreaOverlap(interpolator, 0.8, 10, 5);
	cout << "area overlap = " << areaOverlap << endl;
	double areaCMRep = cm.getAreaOfCMRep();
	cout << "area of CMRep = " << areaCMRep << endl; 
	double check = cm.checkBoundaryFold();
	cout << "check = " << check << endl;
*/	
	return 0;
}

