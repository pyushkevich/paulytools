#include "NewtonMethod.h"
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
			-28.42,
			61.15, 0.38,
			-5.30, 4.67,
			3.62, 0.67, -0.56, 1.89};
	WaveletRep fx(xCoeff, jMax);
	
	double yCoeff[] = {
			-26.31,
			29.86, 0.01,
			3.37, 2.56,
			1.11, -0.39, 0.07, 2.15};
	WaveletRep fy(yCoeff, jMax);
	
	double rhoCoeff[] = {
			-0.01,
			0.0, 0.0,
			0.0, 0.0,
			0.0, 0.0, 0.0, 0.0};
	WaveletRep frho(rhoCoeff, jMax);
	
	c.buildEquationSystem(fx, fy, frho);
	
	NewtonMethod solver(c);
	double *phi = new double[dim + 3];
	for (int i = 0; i < dim + 3; ++i) {
		phi[i] = 1.0;
	}
	solver.initialize(phi);
	double error = solver.solve(c);
	solver.getSolution(phi);
	cout << "solution = " << endl;
	for (int i = 0; i < dim + 3; ++i) {
		cout << phi[i] << endl;
	}
	cout << "error = " << error << endl;
	
	delete[] phi;
	
	return 0;
}

