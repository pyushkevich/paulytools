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
#include <stdlib.h>
#include <memory.h>
#include <iostream>

#include <support.h>
#include <ConjugateGradientMethod.h>
#include <BrentLinearMethod.h>

using namespace std;

// Begin namespace
NAMESPACE_PAULY_START

/*****************************************************************
 ConjugateGradientMethod
 ****************************************************************/
ConjugateGradientMethod::ConjugateGradientMethod(DifferentiableFunction &problem,const Vector &v) : 
p(problem),current(v),xi(v.size()),fTolerance(1.0e-5),bestEver(v)
{
	brent = NULL;

  // Set the step size to a default
  xStepSize = 1.0;
  xBrentTolerance = 2.0e-4;
  flagAdaptStepSize = false;
	
	// Since no optional parameters affect initialization, we perform
	// initialization right here and now.
	initialize();
}

ConjugateGradientMethod::~ConjugateGradientMethod() {
	if(brent) {
		delete brent;
		brent = NULL;
	}
}

void ConjugateGradientMethod::initialize() 
{
	// Compute gradient and function
	current.setValue(p.computeOneJet(current.x,xi));
   
	// Initialize these funny vectors
	g = -xi;h = g;xi = g; 

	// Set best ever
	bestEver = current;

  // Set the working step size to the initial value
  xWorkingStepSize = xStepSize;

	finished = false;

	if(brent) {
		delete brent;
		brent = NULL;
	}
}

void ConjugateGradientMethod::performIteration() {
	// Are we in the brent loop?
	if(brent) {
    
    // Check the current state of the linear minimizer
    NumericalSolution sol(brent->getMinimum());
    sol.setValue(brent->getFAtMinimum());

    // Update the best ever solution
    if(bestEver.getValue() > sol.getValue())
      bestEver = sol;

    // If Brent is done, delete it and check for convergence
		if(brent->isFinished()) {

      // Determine the 'guess' step size for the next iteration
      xWorkingStepSize = xi.dotProduct( brent->getMinDirection() ) / xi.dotProduct(xi);
			xi = brent->getMinDirection();

			delete brent;
			brent = NULL;

			// Check for normal return
			if(2.0 * fabs(sol.getValue() - current.getValue()) <=
				fTolerance*(fabs(sol.getValue())+fabs(current.getValue())+1.0e-10))
			{
				// The solution is equally good at current point and next point...
				finished = true;
				return;
			}

			// Update the current solution
			current = sol;
			current.setValue(p.computeOneJet(current.x,xi));

			// Compute gg, dgg
			double gg = g.dotProduct(g);
			double dgg = xi.dotProduct(xi+g);		// Use xi.dotProduct(xi) for Fletcher-Reeves

			// If gradient is zero, bail
			if(gg==0.0) {
				finished = true;
				return;
			}

			double gam = dgg/gg;

			// Update the vectors
			g = -xi;
			h = g + gam*h;
			xi = h;
		}
		else {
			brent->performIteration();
		}
	}
	else {
    // Choose which step size to use
    double xStep = (flagAdaptStepSize) ? xWorkingStepSize : xStepSize;
    
    // Scale the gradient by the step size before entering Brent's method.
		brent = new BrentLinearMethod(&p,current.x, xStep * xi);
    brent->setTolerance(xBrentTolerance);
	}
}

// End namespace
NAMESPACE_PAULY_END
