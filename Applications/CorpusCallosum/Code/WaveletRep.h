
#ifndef _WaveletRep_H
#define _WaveletRep_H

#include "FunctionRep.h"

class WaveletRep : public FunctionRep {
	private:
	int jMax;
	const int dim;
	double *coeff;
	
	public:
	WaveletRep (const double _coeff[], const int _dim);
	~WaveletRep ();
	
	void setCoeff (const double _coeff[]);
	virtual double get (const double x) const;
	virtual double get1stDeriv (const double x) const;
	virtual double get2ndDeriv (const double x) const;
};

#endif

