#ifndef _POWELL_H_
#define _POWELL_H_

#include <optima.h>
#include <BrentLinearMethod.h>

// Begin namespace
NAMESPACE_PAULY_START

class PowellMethod : public NumericalMethod {
private:
	// Temporary variables
	int i,ibig,j,iter,n;
	double del,fp,fptt,t,ftol;
	Vector pt,ptt,xit;
	Matrix xi;

	// Position flags
	bool done,inSubLoop;

	// Linear minimizer
	BrentLinearMethod *brent;

	// Function pointer
	Function *problem;

	// Current solution
	NumericalSolution X;

	// Constants
	static const int ITMAX;
	static const double TINY;

public:
	// Constructor
	PowellMethod(Function *problem,const Vector Xstart);

	// Destructor
	virtual ~PowellMethod();

	// Set a different direction matrix
	void setDirectionMatrix(const Matrix &M);

	// Set a different tolerance parameter
	void setFTolerance(double ftol);

	// Perform a single iteration
	void performIteration();

	// Check if the method is finished
	bool isFinished();

	// Get the function
	Function &getFunction();

	// Get best ever position
	const Vector &getBestEverX();

	// Get best ever value
	double getBestEverValue();
};

// End namespace
NAMESPACE_PAULY_END

#endif
