/******************************************************************
 * OPTIMA Library                                                 *
 ******************************************************************
 * Author:						Paul Yushkevich
 *
 * Date:							Apr 1, 1999
 *
 * Description					Multidimensional optimization algorithms
 *									See http://www.cs.unc.edu/~pauly/optima
 *									
 *	Sources:						"Numerical Recepies in C", 
 *									Michaelewitz, "Genetic Algorithms + Data
 *									Structures = Evolutionary Programs"
 *
 * Dependencies:				PY Matrix library, CLAPACK
 ******************************************************************
 * ConjugateGradientMethod.h
 *	-------------------------
 * This method from NRC-10.6 allows us to optimize nicely using first
 * partial derivatives in n-space.
 ******************************************************************/
#ifndef _OPTIMA_CONJGRAD_
#define _OPTIMA_CONJGRAD_

// Include files
#include <optima.h>
#include <BrentLinearMethod.h>

// Begin namespace
NAMESPACE_PAULY_START

/*******************************************************************
 Optimize f(x) using conjugate gradient descent 

 Usage:
	-	Create a ConjugateGradientMethod (with a function and a sample x)
	-	Set options with setXXX methods
	-  Call performIteration() until isFinished==true
	-	Check results with getBestEverXXX()
 ******************************************************************/
class ConjugateGradientMethod : public NumericalMethod {
protected:
	NumericalSolution bestEver,current;
	DifferentiableFunction &p;
	BrentLinearMethod *brent;

	// Input parameters
	double fTolerance;

	// Variables used between iterations
	Vector g,h,xi;

	// Initialize
	void initialize();

	// Are we finished
	bool finished;

public:
	// Initialize the method
	// ConjugateGradientMethod(DifferentiableFunction &problem);

	// Initialize the method
	ConjugateGradientMethod(DifferentiableFunction &problem,const Vector &v);
	virtual ~ConjugateGradientMethod();

	// Start the method, supply a new problem
	void setProblem(DifferentiableFunction &problem,const Vector &v) {
		p = problem;
		setStartPoint(v);
	}

	// Start the method - for first constructor
	void setStartPoint(const Vector &start) {
		current.x = start;
		initialize();
	}

	// Set the value for epsilon, the terminating condition
	void setTolerance(double fTol) {
		fTolerance = fTol;
	}

   // Perform an iteration of the method
   void performIteration();

	// Check termination condition (termination determined by tolerance)
	bool isFinished() { return finished; };
   
	// Solution report methods
	double getBestEverValue() {return bestEver.getValue();}
   const Vector &getBestEverX() {return bestEver.x;}

	// Get the problem on which we are operating
	virtual Function &getFunction() {
		return p;
	}
};

// End namespace
NAMESPACE_PAULY_END

#endif //_OPTIMA_CONJGRAD_

