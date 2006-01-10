
#ifndef _SumBSplineCubUni_H
#define _SumBSplineCubUni_H

#include "BSplineCubUni.h"

class SumBSplineCubUni : public BSplineCubUni {
	private:
	BSplineCubUni b1;
	
	public:
	SumBSplineCubUni (const double coeff1[], const int np1,const double coeff2[], const int np2);
	~SumBSplineCubUni () {}
	
	double get (const double x) const;
	double get1stDeriv (const double x) const;
	double get2ndDeriv (const double x) const;
};

#endif

