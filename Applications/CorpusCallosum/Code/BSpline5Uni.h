// Cubic uniform B-spline with precalculated blending function.

#ifndef _BSpline5Uni_H
#define _BSpline5Uni_H

#include "FunctionRep.h"

class BSpline5Uni : public FunctionRep {
	private:
	static const double M[6][6];
	int np;
	double *ctrlPt;
	
	public:
	BSpline5Uni () {np=1; ctrlPt[0] = 0;}
	BSpline5Uni (const double _ctrlPt[], const int _np);
	~BSpline5Uni ();

	void setCtrlPt (const double _ctrlPt[]);
	virtual double get (const double x) const;
	virtual double get1stDeriv (const double x) const;
	virtual double get2ndDeriv (const double x) const;
};

#endif
