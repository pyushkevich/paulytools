#include <cmath>
#include "Wavelet.h"

const double Wavelet::MUL = sqrt(sqrt(M_PI))/sqrt(3.0);

double Wavelet::psi (const double x) {
	double xx = 2*x - 1;
	xx *= xx;
	
	double result = 1 - xx;
	result *= MUL;
	result *= 2.0;
	result *= exp(-0.5*xx);

	return result;
}

double Wavelet::psi1stDeriv (const double x) {
	const double x1 = 2*x - 1;
	const double x2 = x1*x1;
	double result = x2 - 3;
	result *= MUL;
	result *= exp(-0.5*x2);
	result *= 4.0;
	result *= x1;
	
	return result;
}

double Wavelet::psi2ndDeriv (const double x) {
	const double x1 = 2*x - 1;
	const double x2 = x1*x1;
	double result = -3;
	result += 6 * x2;
	result -= x2 * x2;
	result *= MUL;
	result *= exp(-0.5*x2);
	result *= 8.0;
	
	return result;
}

double Wavelet::wavelet (const double x, const int j, const int k) {
	int j2 = 2;
	j2 <<= j - 1;
	return psi(j2*x - k);
}

double Wavelet::wavelet1stDeriv (const double x, const int j, const int k) {
	int j2 = 2;
	j2 <<= j - 1;
	double result = psi1stDeriv(j2*x - k);
	result *= j2;
	
	return result;
}

double Wavelet::wavelet2ndDeriv (const double x, const int j, const int k) {
	int j2 = 2;
	j2 <<= j - 1;
	int j4 = 2;
	j4 <<= 2*j - 1;
	double result = psi2ndDeriv(j2*x - k);
	result *= j4;

	return result;
}

