#include <iostream>
#include <cstdlib>
#include "TemplatedWaveletRep.h"

using namespace std;

void run (const WaveletRep& frep, const double x) {
	cout << frep.get(x) << endl;
	cout << frep.get1stDeriv(x) << endl;
	cout << frep.get2ndDeriv(x) << endl;
}

int main (int argc, char *argv[]) {
	
	if (argc < 2) {
		cerr << "Not enough arguments" << endl;
		cerr << "Usage: " << argv[0] << " x" << endl;
		exit (1);
	}
	
	const int tJMax = 3;
	const int tDim = 17;
	const int dJMax = 2;
	const int dDim = 9;
	
	const double tCoeff[tDim] = {
		-28.2364,
		60.1971, -0.377662,
		-4.78251, 4.24725,
		4.21383, 1.00549, -0.33139, -1.87177,
		0.568533, -0.603949, 0.154689, -0.1541, 0.072736, -0.123597, 0.104211, -0.666733,
	};
	
	const double dCoeff[dDim] = {
		0.0,
		0.0, 0.0,
		0.0, 0.0,
		0.0, 0.0, 0.0, 0.0,
	};
	
	TemplatedWaveletRep frep(tCoeff, tJMax, dCoeff, dJMax);
	double x = atof(argv[1]);
	run(frep, x);
	
	return 0;
}

