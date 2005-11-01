#include <iostream>
#include <cstdlib>
#include "BSplineCubUni.h"

using namespace std;

int main (int argc, char *argv[]) {
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " x" << endl;
		exit (1);
	}
	
	double x = atof(argv[1]);
	
	const int jMax = 6;
	const double coeff[] = {
			0.0,
			0.0, 0.0,
			0.0, 0.0,
			1.0, 0.0, 1.0, 0.0};
	BSplineCubUni fx(ctrlPt, np);
	cout << "bspline = " << fx.get(x) << endl;
	cout << "bspline1stDeriv = " << fx.get1stDeriv(x) << endl;
	cout << "bspline2ndDeriv = " << fx.get2ndDeriv(x) << endl;
	
	return 0;
}

