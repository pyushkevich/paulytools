#include <cmath>
#include "BSplineCubUni.h"

const double* BSplineCubUni::M ={
	{1.0/6.0, -0.5, 0.5, -1.0/6.0},
	{2.0/3.0, 0.0, -1.0, 0.5},
	{1.0/6.0, 0.5, 0.5, -0.5},
	{0.0, 0.0, 0.0, 1.0/6.0}
};

BSplineCubUni::BSplineCubUni (const double _ctrlPt[], const int _np) : np(_np) {
	ctrlPt = new double [np];
	setCtrlPt(_ctrlPt);
}
BSplineCubUni::~BSplineCubUni () {
	delete[] ctrlPt;
}
	
void BSplineCubUni::setCtrlPt (const double _ctrlPt[]) {
	for (int i = 0; i < np; ++i) {
		ctrlPt[i] = _ctrlPt[i];
	}
}

double BSplineCubUni::get (const double x) const {
	int i = (int) x*(np-3);
	double tmp = 0.0;
	for (int j = 0; j < 4; j++) {
		double tmp1 = M[j][0];
		tmp1 += x*M[j][1];
		tmp1 += x*x*M[j][2];
		tmp1 += x*x*x*M[j][3];
		tmp1 *= ctrlPt[i+j];
		tmp += tmp1;
	}
	return tmp;
}
double BSplineCubUni::get1stDeriv (const double x) const {
	int i = (int) x*(np-3);
	double tmp = 0.0;
	for (int j = 0; j < 4; j++) {
		double tmp1 = M[j][1];
		tmp1 += 2*x*M[j][2];
		tmp1 += 3*x*x*M[j][3];
		tmp1 *= ctrlPt[i+j];
		tmp += tmp1;
	}
	return tmp;
}
double BSplineCubUni::get2ndDeriv (const double x) const {
	int i = (int) x*(np-3);
	double tmp = 0.0;
	for (int j = 0; j < 4; j++) {
		double tmp1 = 2*M[j][2];
		tmp1 += 6*x*M[j][3];
		tmp1 *= ctrlPt[i+j];
		tmp += tmp1;
	}
	return tmp;
}