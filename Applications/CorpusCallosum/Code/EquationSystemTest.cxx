#include "EquationSystem.h"
#include "WaveletRep.h"
#include <iostream>
#include <cstdlib>

using namespace std;

int main (int argc, char *argv[]) {
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " dim" << endl;
		exit(1);
	}
	const int dim = atoi(argv[1]);
	EquationSystem c(dim);
	
	const int jMax = 2;
	double xCoeff[] = {
			0.0,
			0.0, 0.0,
			0.0, 0.0,
			1.0, 0.0, 1.0, 0.0};
	WaveletRep fx(xCoeff, jMax);
	
	double yCoeff[] = {
			0.0,
			0.0, 0.0,
			0.0, 0.0,
			1.0, 1.0, 1.0, 1.0};
	WaveletRep fy(yCoeff, jMax);
	
	double rhoCoeff[] = {
			0.0,
			0.0, 0.0,
			0.0, 0.0,
			1.0, 0.0, 1.0, 0.0};
	WaveletRep frho(rhoCoeff, jMax);
	
	c.buildEquationSystem(fx, fy, frho);
	cout << "EquationSystem = " << endl;
	cout << c << endl;
	
	return 0;
}

