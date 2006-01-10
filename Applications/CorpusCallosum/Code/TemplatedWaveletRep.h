
#ifndef _TemplatedWaveletRep_H
#define _TemplatedWaveletRep_H

#include "WaveletRep.h"

class TemplatedWaveletRep : public WaveletRep {
	private:
	WaveletRep templateShape;
	
	public:
	TemplatedWaveletRep (const double tCoeff[], const int tdim, const double dCoeff[], const int ddim);
	~TemplatedWaveletRep () {}
	
	double get (const double x) const;
	double get1stDeriv (const double x) const;
	double get2ndDeriv (const double x) const;
};

#endif

