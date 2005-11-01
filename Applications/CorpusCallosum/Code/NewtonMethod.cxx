#include "NewtonMethod.h"
#include <cmath>

using namespace std;

NewtonMethod::NewtonMethod (const EquationSystem& c) : dim(c.dim), A(c.dim) {
	iteration = 20;
	tolerance = pow(10.0, -14.0);
	b = new double [dim + 3];
	x = new double [dim + 3];
	
	// fill up sparse matrix A with constant terms
	for (int i = 0; i < dim + 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			A.updateElement(i, j, c.coeff1stOrder[i][j]);
		}
	}
}

NewtonMethod::~NewtonMethod () {
	delete[] b;
	delete[] x;
}

void NewtonMethod::initialize (const double *phi) {
	for (int i = 0; i < dim + 3; ++i) {
		x[i] = phi[i];
	}
}

void NewtonMethod::updateb (const EquationSystem& c) {
	// first element
	b[0] = -c.coeffConst[0];
	b[0] -= c.coeff1stOrder[0][1] * x[1];
	double tmp = x[2] - x[0];
	tmp *= tmp;
	tmp *= c.coeff2ndOrder[0];
	b[0] -= tmp;
	
	// last element
	b[dim + 2] = -c.coeffConst[dim + 2];
	b[dim + 2] -= c.coeff1stOrder[dim + 2][1] * x[dim + 1];
	tmp = x[dim + 2] - x[dim];
	tmp *= tmp;
	tmp *= c.coeff2ndOrder[1];
	b[dim + 2] -= tmp;
	
	// element in between
	for (int i = 1; i < dim + 2; ++i) {
		b[i] = -c.coeffConst[i];
		for (int j = 0; j < 3; ++j) {
			b[i] -= c.coeff1stOrder[i][j] * x[i + j - 1];
		}
	}
}

void NewtonMethod::updateA (const EquationSystem& c) {
	double tmp = 2.0;
	tmp *= c.coeff2ndOrder[0];
	tmp *= x[0] - x[2];
	A.updateElement(0, 0, tmp);
	A.updateElement(0, 2, -tmp);
	
	tmp = 2.0;
	tmp *= c.coeff2ndOrder[1];
	tmp *= x[dim] - x[dim + 2];
	A.updateElement(dim + 2, 0, tmp);
	A.updateElement(dim + 2, 2, -tmp);
}

double NewtonMethod::solve (const EquationSystem& c) {
	double error;
	UnsymmetricRealPARDISO solver;
	double *dx = new double [dim + 3];
	for (int counter = 0; counter < iteration; ++counter) {
		// compute the x dependent elements in A and b
		updateb(c);
		updateA(c);
		
		// solve for the dx vector
		solver.SymbolicFactorization(dim + 3, A.rowIndex, A.columnIndex, A.dataArray);
		solver.NumericFactorization(A.dataArray);
		solver.Solve(b, dx);
		
		// update the x vector using dx
		for(int i = 0; i < dim + 3; ++i){
			x[i] += dx[i];
		}
		
		// compute the error
		error = 0.0;
		for (int i = 0; i < dim + 3; ++i) {
			if(fabs(b[i]) > error) {
				error = fabs(b[i]);
			}
		}
		if (error < tolerance) break;
		
	}
	
	delete[] dx;
	return error;
}

void NewtonMethod::getSolution (double* solution) const {
	for (int i = 0; i < dim + 3; ++i) {
		solution[i] = x[i];
	}
}

