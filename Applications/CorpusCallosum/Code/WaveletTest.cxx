#include <iostream>
#include <cstdlib>
#include "Wavelet.h"

using namespace std;

int main (int argc, char *argv[]) {
	if (argc < 4) {
		cerr << "Usage: " << argv[0] << " j k x" << endl;
		exit (1);
	}
	
	int j = atoi(argv[1]);
	int k = atoi(argv[2]);
	double x = atof(argv[3]);
	
	cout << "psi = " << Wavelet::psi(x) << endl;
	cout << "wavelet = " << Wavelet::wavelet(x, j, k) << endl;
	cout << "psi1stDeriv = " << Wavelet::psi1stDeriv(x) << endl;
	cout << "wavelet1stDeriv = " << Wavelet::wavelet1stDeriv(x, j, k) << endl;
	cout << "psi2ndDeriv = " << Wavelet::psi2ndDeriv(x) << endl;
	cout << "wavelet2ndDeriv = " << Wavelet::wavelet2ndDeriv(x, j, k) << endl;
	
	return 0;
}

