#include <iostream>
#include <cstdlib>
#include "SumBSplineCubUni.h"

using namespace std;

int main (int argc, char *argv[]) {
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " x" << endl;
		exit (1);
	}
	
	double x = atof(argv[1]);
	
	const double ctrlPt[] ={
	  56.0068, 58.6106, 68.3392, 62.8632, 49.797, 37.2775, 23.3805, 9.1962,
	  -5.46414, -19.7429, -33.4751, -46.5852, -59.3993, -69.4148, -69.686,-68.479
	};
	const int np = 16;
	double dx = 1.0e-5;
	SumBSplineCubUni fx(ctrlPt, np, ctrlPt, np);
	double dfx = (fx.get(x+dx)-fx.get(x-dx))/2/dx;
	double ddfx = (fx.get(x+dx)+fx.get(x-dx)-2*fx.get(x))/dx/dx;
	cout << "bspline = " << fx.get(x) << endl;
	cout << "bspline1stDeriv = " << fx.get1stDeriv(x) << endl;
	cout << "bspline2ndDeriv = " << fx.get2ndDeriv(x) << endl;
	cout << "dfx = " << dfx << endl;
	cout << "ddfx = " << ddfx << endl;
	return 0;
}

