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

#include <support.h>
#include <ConjugateGradientMethod.h>
#include <BrentLinearMethod.h>

// Begin namespace
NAMESPACE_PAULY_START

/*****************************************************************
 ConjugateGradientMethod
 ****************************************************************/
ConjugateGradientMethod::ConjugateGradientMethod(DifferentiableFunction &problem,const Vector &v) : 
p(problem),current(v),xi(v.size()),fTolerance(1.0e-5),bestEver(v)
{
	brent = NULL;
	
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

void ConjugateGradientMethod::initialize() {
	// Compute gradient and function
	current.setValue(p.computeOneJet(current.x,xi));
   
	// Initialize these funny vectors
	g = -xi;h = g;xi = g; 

	// Set best ever
	bestEver = current;

	finished = false;

	if(brent) {
		delete brent;
		brent = NULL;
	}
}

void ConjugateGradientMethod::performIteration() {
	// Are we in the brent loop?
	if(brent) {
		if(brent->isFinished()) {
			NumericalSolution next(brent->getMinimum());
			next.setValue(brent->getFAtMinimum());
			xi = brent->getMinDirection();

			delete brent;
			brent = NULL;

			// Update the best ever solution
			if(bestEver.getValue() > next.getValue())
				bestEver = next;

			// Check for normal return
			if(2.0 * fabs(next.getValue() - current.getValue()) <=
				fTolerance*(fabs(next.getValue())+fabs(current.getValue())+1.0e-10))
			{
				// The solution is equally good at current point and next point...
				finished = true;
				return;
			}

			// Update the current solution
			current = next;
			current.setValue(p.computeOneJet(current.x,xi));
			//xi.normalize();

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
		brent = new BrentLinearMethod(&p,current.x,0.001 * xi);
	}
}

// End namespace
NAMESPACE_PAULY_END
