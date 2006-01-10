#include "WaveletRep.h"
#include "Wavelet.h"
#include <iostream>
using namespace std;

WaveletRep::WaveletRep (const double _coeff[], const int _dim) : dim(_dim) {
	// size of coeff should be 2^(jMax + 1) + 1
	int tmp = dim-1;
	jMax = 0;
	for (; tmp > 1; jMax++) {
	  tmp /= 2;
	}
	jMax -= 1;

	coeff = new double [dim];
	setCoeff(_coeff);
}

void WaveletRep::setCoeff (const double _coeff[]) {
	for (int i = 0; i < dim; ++i) {
		coeff[i] = _coeff[i];
	}
}

WaveletRep::~WaveletRep () {
	delete[] coeff;
}

double WaveletRep::get (const double x) const {
	// the constant term
	double sum = coeff[0];
	
	// the mother wavelet basis with translational freedom
	double tmp = Wavelet::psi(x - coeff[2]);
	tmp*= coeff[1];
	sum += tmp;
	
	// the sum over the daughter basis
	int index = 3;
	for (int j = 1; j <= jMax; ++j) {
		int kMax = 1;
		kMax <<= j;
		for (int k = 0; k < kMax; ++k) {
			double tmp = Wavelet::wavelet(x, j, k);
			tmp *= coeff[index];
			sum += tmp;
			++index;
		}
	}
	
	return sum;
}

double WaveletRep::get1stDeriv (const double x) const {
	double sum = 0.0;
	
	// derivatives of the mother wavelet basis with translational freedom
	double tmp = Wavelet::psi1stDeriv(x - coeff[2]);
	tmp*= coeff[1];
	sum += tmp;
	
	// derivatives of the sum over the daughter basis
	int index = 3;
	for (int j = 1; j <= jMax; ++j) {
		int kMax = 1;
		kMax <<= j;
		for (int k = 0; k < kMax; ++k) {
			double tmp = Wavelet::wavelet1stDeriv(x, j, k);
			tmp *= coeff[index];
			sum += tmp;
			++index;
		}
	}
	
	return sum;
}

double WaveletRep::get2ndDeriv (const double x) const {
	double sum = 0.0;
	
	// derivatives of the mother wavelet basis with translational freedom
	double tmp = Wavelet::psi2ndDeriv(x - coeff[2]);
	tmp*= coeff[1];
	sum += tmp;
	
	// derivatives of the sum over the daughter basis
	int index = 3;
	for (int j = 1; j <= jMax; ++j) {
		int kMax = 1;
		kMax <<= j;
		for (int k = 0; k < kMax; ++k) {
			double tmp = Wavelet::wavelet2ndDeriv(x, j, k);
			tmp *= coeff[index];
			sum += tmp;
			++index;
		}
	}
	
	return sum;
}

