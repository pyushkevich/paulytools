#include <iostream>
#include <cstdlib>
#include "WaveletRep.h"

using namespace std;

int main (int argc, char *argv[]) {
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " x" << endl;
		exit (1);
	}
	
	double x = atof(argv[1]);
	
	const int jMax = 2;
	const double coeff[] = {
			0.0,
			0.0, 0.0,
			0.0, 0.0,
			1.0, 0.0, 1.0, 0.0};
	WaveletRep fx(coeff, jMax);
	cout << "waveletRep = " << fx.get(x) << endl;
	cout << "waveletRep1stDeriv = " << fx.get1stDeriv(x) << endl;
	cout << "waveletRep2ndDeriv = " << fx.get2ndDeriv(x) << endl;
	
	return 0;
}

