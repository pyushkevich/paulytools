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
 * Powell.h
 *	-------------------------
 * This method from NRC-10.5 allows us to optimize nicely w/o first
 * partial derivatives in n-space.
 ******************************************************************/
#include <stdlib.h>
#include <memory.h>

#include <support.h>
#include <Powell.h>

// Begin namespace
NAMESPACE_PAULY_START

inline double SQR(double a) {
	return a*a;
}

const int PowellMethod::ITMAX = 200;
const double PowellMethod::TINY = 1.0e-25;

PowellMethod::PowellMethod(Function *problem,const Vector Xstart)
: X(Xstart)
{
	// Initialize the parameters
	this->problem = problem;

	// Initialize the matrix XI to identity
	n = Xstart.rows();
	xi.setSize(n,n);
	xi = 1.0;

	// Initialize the other vectors
	pt = X.x;
	ptt.setSize(n);
	xit.setSize(n);
	X.setValue(problem->evaluate(X.x));

	// Tolerance initial value
	ftol = 1.0e-5;

	// Initialize the flags and pointers
	brent = NULL;
	done = false;
	iter = 0;
}

PowellMethod::~PowellMethod() {
	if(brent)
		delete brent;
}

void PowellMethod::setDirectionMatrix(const Matrix &M) {
	xi = M;
}

void PowellMethod::setFTolerance(double ftol) {
	this->ftol = ftol;
}

void PowellMethod::performIteration() {
	if(done)
		return;
	
	if(brent) {
		if(brent->isFinished()) {
			// Store the results of the linmin
			// X.value = brent->getFAtMinimum();
			xit = brent->getMinDirection();
			// X.x += xit;

			// Clear the linmin
			delete brent;
			brent = NULL;

			// Where in the algorithm are we?
			if(inSubLoop) {
				if(fptt-X.value > del) {
					del = fptt-X.value;
					ibig = i;
				}
				if(++i < n) {
					// Stay in the loop
					xit = xi.getColumn(i);
					fptt = X.value;
					brent = new BrentLinearMethod(problem,X.x,xit);
				}
				else {
					// Do the after-the loop stuff
					printf("G ");
					inSubLoop = false;

					// Check for termination
					if(2.0*(fp-X.value) <= ftol * (fabs(fp)+fabs(X.value))+TINY) {
						done = true;
						return;
					}

					// Construct the exterpolated point and the average direction moved
					ptt = 2.0 * X.x - pt;
					xit = X.x - pt;
					pt = X.x;
					fptt = problem->evaluate(ptt);

					// The if loop
					if(fptt < fp) {
						t = 2.0 * (fp-2.0*X.value+fptt)*SQR(fp-X.value-del)-del*SQR(fp-fptt);
						if(t < 0.0) {
							brent = new BrentLinearMethod(problem,X.x,xit);
						}
					}
				}
			} 
			else {
				// Other case of brent
				xi.insertMatrix(0,ibig,xi.getColumn(n-1));
				xi.insertMatrix(0,n-1,xit);
			}
		}
		else {
			// Brent not finished
			brent->performIteration();
			X.x = brent->getMinimum();
			X.value = brent->getFAtMinimum();
		}
	} 
	else {
		// Increase the iteration counter
		if(++iter > ITMAX) {
			done = true;
		}
		else {
			// There is no brent: this is hit in the beginning of a loop
			fp = X.value;
			ibig = -1;
			del = 0.0;
			i = 0;
			inSubLoop = true;
			xit = xi.getColumn(i);
			fptt = X.value;
			brent = new BrentLinearMethod(problem,X.x,xit);
		}
	}
}

bool PowellMethod::isFinished() {
	return done;
}

Function &PowellMethod::getFunction() {
	return *problem;
}

const Vector &PowellMethod::getBestEverX() {
	return X.x;
}

double PowellMethod::getBestEverValue() {
	return X.value;
}

// End namespace
NAMESPACE_PAULY_END
