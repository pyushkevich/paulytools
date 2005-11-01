
#ifndef _WaveletRep_H
#define _WaveletRep_H

#include "FunctionRep.h"

class WaveletRep : public FunctionRep {
	private:
	const int jMax;
	int dim;
	double *coeff;
	
	void setJMax (const int _jMax);
	
	public:
	WaveletRep (const double _coeff[], const int _jMax);
	~WaveletRep ();
	
	void setCoeff (const double _coeff[]);
	virtual double get (const double x) const;
	virtual double get1stDeriv (const double x) const;
	virtual double get2ndDeriv (const double x) const;
};

#endif

