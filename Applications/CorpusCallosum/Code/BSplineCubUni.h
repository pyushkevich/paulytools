// Cubic uniform B-spline with precalculated blending function.

#ifndef _BSplineCubUni_H
#define _BSplineCubUni_H

#include "FunctionRep.h"

class BSplineCubUni : public FunctionRep {
	private:
	const double* M[4][4];
	int np;
	double *ctrlPt;
	
	public:
	BSplineCubUni (const double _ctrlPt[], const int _np);
	~BSplineCubUni ();
	
	void setCtrlPt (const double _ctrlPt[]);
	double get (const double x) const;
	double get1stDeriv (const double x) const;
	double get2ndDeriv (const double x) const;
};

#endif
