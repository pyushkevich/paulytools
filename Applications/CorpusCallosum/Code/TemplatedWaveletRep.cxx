#include "TemplatedWaveletRep.h"

TemplatedWaveletRep::TemplatedWaveletRep (const double tCoeff[], const int tJMax, const double dCoeff[], const int dJMax) : templateShape(tCoeff, tJMax), WaveletRep(dCoeff, dJMax) {
}

double TemplatedWaveletRep::get (const double x) const {
	double value = templateShape.get(x);
	value += WaveletRep::get(x);
	return value;
}

double TemplatedWaveletRep::get1stDeriv (const double x) const {
	double value = templateShape.get1stDeriv(x);
	value += WaveletRep::get1stDeriv(x);
	return value;
}

double TemplatedWaveletRep::get2ndDeriv (const double x) const {
	double value = templateShape.get2ndDeriv(x);
	value += WaveletRep::get2ndDeriv(x);
	return value;
}

