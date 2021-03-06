// Cubic uniform B-spline with precalculated blending function.

#ifndef _BSplineCubUni_H
#define _BSplineCubUni_H

#include "FunctionRep.h"

class BSplineCubUni : public FunctionRep {
	private:
	static const double M[4][4];
	int np;
	double *ctrlPt;
	
	public:
	BSplineCubUni () {np=1; ctrlPt[0] = 0;}
	BSplineCubUni (const double _ctrlPt[], const int _np);
	~BSplineCubUni ();

	void setCtrlPt (const double _ctrlPt[]);
	virtual double get (const double x) const;
	virtual double get1stDeriv (const double x) const;
	virtual double get2ndDeriv (const double x) const;
};

#endif
