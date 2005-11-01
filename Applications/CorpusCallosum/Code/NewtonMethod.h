
#ifndef _NewtonMethod_H
#define _NewtonMethod_H

#include "EquationSystem.h"
#include "SparseMatrix.h"
#include "PardisoInterface.h"

class NewtonMethod {
	private:
	const int dim;
	int iteration;
	double tolerance;
	// solution vector
	double *x;
	// error vector
	double *b;
	// linear system in sparse matrix form
	SparseMatrix A;
	
	void updateA (const EquationSystem& c);
	void updateb (const EquationSystem& c);
	
	public:
	NewtonMethod(const EquationSystem& c);
	~NewtonMethod();
	
	void initialize (const double *phi);
	double solve (const EquationSystem& c);
	void getSolution (double *solution) const;
};

#endif

