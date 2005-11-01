#include "MedialAtom.h"
#include <iostream>
#include <cstdlib>

using namespace std;

int main (int argc, char *argv[]) {
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " dim" << endl;
		exit (1);
	}
	
	const int dim = atoi(argv[1]);
	
	MedialAtom m[dim + 1];
	double fx[dim + 1];
	double fy[dim + 1];
	double fxt[dim + 1];
	double fyt[dim + 1];
	double phi[dim + 3];
	
	phi[0] = 1;
	phi[dim + 3] = 1;
	
	for (int i = 0; i < dim + 1; ++i) {
		fx[i] = i;
		fy[i] = i;
		fxt[i] = 1;
		fyt[i] = 1;
		phi[i + 1] = 1;
	}
	
	MedialAtom::phiToMedialAtoms(fx, fy, fxt, fyt, phi, dim, m);
	
	for (int i = 0 ; i < dim + 1; ++i) {
		cout << "MedialAtom " << i << endl;
		cout << m[i] << endl;
	}
	
	return 0;
}

