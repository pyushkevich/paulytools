#include <iostream>
#include <cstdlib>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkDiscreteGaussianImageFilter.h"

#include "CMRep2D.h"
#include "TemplatedWaveletRep.h"

#include "/mnt/data1/huiz/research/eclipse/workspace/Mala_3D/src/numerics/numerics.h"

using namespace std;
using namespace numerics;

typedef float PixelType;
typedef itk::Image< PixelType, Dimension >	ImageType;
typedef itk::ImageFileReader< ImageType >	ImageReader;
typedef itk::LinearInterpolateImageFunction< ImageType, double > ImageInterpolator;
typedef itk::DiscreteGaussianImageFilter< ImageType, ImageType > ImageSmoother;

double compute (double p[]);
void grad (double p[], double xi[]);
double both (double p[], double xi[]);

// globally defined variables
const int tJMax = 3;
const int tWDim = 17;
const int dJMax = 2;
const int dWDim = 9;
ImageInterpolator::Pointer interpolator = ImageInterpolator::New();
double areaOfImage = 0.0;
CMRep2D *cm = 0;
double *phi = 0;

const double scale[dWDim*3] = {
	0.1,
	0.1, 0.1,
	0.1, 0.1,
	1.0, 1.0, 1.0, 1.0,
	
	0.1,
	0.1, 0.1,
	0.1, 0.1,
	1.0, 1.0, 1.0, 1.0,
	
	0.01,
	0.01, 0.01,
	0.01, 0.01,
	0.1, 0.1, 0.1, 0.1
};

double pv[dWDim*3];

const double p0[tWDim*3] = {
	-28.2364,
	60.1971, -0.377662,
	-4.78251, 4.24725,
	4.21383, 1.00549, -0.33139, -1.87177,
	0.568533, -0.603949, 0.154689, -0.1541, 0.072736, -0.123597, 0.104211, -0.666733,
	
	-26.1098,
	29.7462, -0.00113348,
	2.65794, 2.91456,
	1.34698, 0.12066, 0.22392, 2.53658,
	-0.783693, 0.131142, 0.28515, 0.0791865, -0.0151148, 0.130991, -0.0664204, 0.20153,
	
	-0.2780468875676537,
	0.8430449632231002, -0.4481581004039171,
    -0.3168084491489501, 0.4390913581310029,
	-0.5411271156693166, -0.09268973864155224, -0.09077243805937582, 0.1103090808314945,
	-0.6752267991252191, -0.02694546649788497, -0.05732010382287213, -0.03807655204148059, -0.1325252333762741, -0.21076001234790967, -0.31639484874831914, -0.39152479130559137
};

int main (int argc, char *argv[]) {
	
	if (argc < 3) {
		cerr << "Not enough arguments" << endl;
		cerr << "Usage: " << argv[0] << " inputImageFile dim" << endl;
		exit (1);
	}
	
	const char *inputImageFile = argv[1];
	
	cout << "input Image File is " << inputImageFile << endl;
	
	ImageReader::Pointer reader = ImageReader::New();
	reader->SetFileName(inputImageFile);
	reader->Update();
	// Gaussian smooth the image
/*	ImageSmoother::Pointer smooth = ImageSmoother::New();
	smooth->SetVariance(4.0);
	smooth->SetInput(reader->GetOutput());
	smooth->Update();
*/	
	// interpolator defined globally
//	interpolator->SetInputImage(smooth->GetOutput());
	interpolator->SetInputImage(reader->GetOutput());
	
	// areaOfImage defined globally
	areaOfImage = 0.0;
	for (int i = 0; i < 256; ++i) {
		for (int j = 0; j < 256; ++j) {
			double index[2] = {i,j};
			areaOfImage += interpolator->Evaluate(index);
		}
	}
	
	const int dim = atoi(argv[2]);
	// phi defined globally
	phi = new double[dim + 3];
	for (int i = 0; i < dim + 3; ++i) {
		phi[i] = 1.0;
	}
	// cm defined globally
	cm = new CMRep2D(dim);
	
	double ftol = 0.0001;
	double pDelta[dWDim*3] = {
		0.0,
		0.0, 0.0,
		0.0, 0.0,
		0.0, 0.0, 0.0, 0.0,
		
		0.0,
		0.0, 0.0,
		0.0, 0.0,
		0.0, 0.0, 0.0, 0.0,
		
		0.0,
		0.0, 0.0,
		0.0, 0.0,
		0.0, 0.0, 0.0, 0.0
	};
	
	Function func(dWDim*3, compute);
	DirectionSetMinimizer dsm(true);
	double min = dsm.run(pDelta, ftol, func);
	
//	Gradient func(dWDim*3, compute, grad, both);
//	ConjugateGradientMinimizer cgm(true);
//	double min = cgm.run(pDelta, ftol, func);
	
	// output
	cout << "xCoeff = {";
	for (int i = 0; i < dWDim - 1; ++i) {
		cout << scale[i] * pDelta[i] << ",";
	}
	cout << scale[dWDim - 1] * pDelta[dWDim - 1] << "}" << endl;
	
	cout << "yCoeff = {";
	for (int i = 0; i < dWDim - 1; ++i) {
		cout << scale[i + dWDim] * pDelta[i + dWDim] << ",";
	}
	cout << scale[2*dWDim - 1] * pDelta[2*dWDim - 1] << "}" << endl;
	
	cout << "rhoCoeff = {";
	for (int i = 0; i < dWDim - 1; ++i) {
		cout << scale[i + 2*dWDim] * pDelta[i + 2*dWDim] << ",";
	}
	cout << scale[3*dWDim - 1] * pDelta[3*dWDim - 1] << "}" << endl;
	
	// clean up
	delete[] phi;
	delete cm;
	
	return 0;
}

int iter = 0;

double compute (double p[]) {
	iter++;
	for (int i = 0; i < dWDim*3; ++i) {
		pv[i] = p[i];
		pv[i] *= scale[i];
	}
	TemplatedWaveletRep fx(p0, tJMax, pv, dJMax);
	TemplatedWaveletRep fy(p0 + tWDim, tJMax, pv + dWDim, dJMax);
	TemplatedWaveletRep frho(p0 + 2*tWDim, tJMax, pv + 2*dWDim, dJMax);
	bool nonnegative = cm->buildCMRep2D(fx, fy, frho, phi);
	if (!nonnegative) {
		return 1.0e20;
	}
	double check = cm->checkBoundaryFold();
	check *= 1000.0;
	double overlap = cm->computeAreaOverlapRatio(interpolator, areaOfImage, 0.8, 10, 5);
	double sum = overlap + check;
	return overlap + check;
} 

double both (double p[], double xi[]) {
	grad(p,xi);
	return compute(p);
}

int iterGrad = 0;

void grad (double p[], double xi[]) {
	++iterGrad;
	double pvv[dWDim * 3];
	for (int i = 0; i < dWDim*3; ++i) {
		pvv[i] = p[i];
	}
	double delta = 0.1;
	for (int i = 0; i < dWDim*3; ++i) {
		const double center = pvv[i];
		// backward
		pvv[i] = center - delta;
		double tmp1 = compute(pvv);
		// forward
		pvv[i] = center + delta;
		double tmp2 = compute(pvv);
		// compute the numerical derivative
		xi[i] = tmp2 - tmp1;
		xi[i] /= 2.0 * delta;
		// reset the value
		pvv[i] = center;
	}
	cout << "iterGrad = " << iterGrad << endl;
}

