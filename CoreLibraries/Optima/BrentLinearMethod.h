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
 * BrentLinearMethod.h
 *	-------------------
 * This method from NRC-10.2 allows us to optimize along a vector in
 * n-space.
 ******************************************************************/
#ifndef _OPTIMA_LINMIN_
#define _OPTIMA_LINMIN_

// Include files
#include <optima.h>

// Begin namespace
NAMESPACE_PAULY_START


/*******************************************************************
 Solve a problem of finiding s that optimizes f(P+sn) 

 Use this by creating a BrentLinearMethod and calling run.  Check
 the values of x,n and 'value' after the algorithm has finished.
 ******************************************************************
class BrentLinearMethod {
protected:
	Function &p;

	double func(double t) {		
		xt = n;
		xt *= t;
		xt += x;
		return p.evaluate(xt);
	}

	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc);
	double brent(double ax, double bx, double cx, double tol,double *xmin);

public:
	// Position of the solution
	Vector x,xt;

	// Vector from old solution to new solution
	Vector n;

	// Optimal valur
	double value;

	// Constructor (does not require a solution space because there is a single
	// defined starting point.
	BrentLinearMethod(Function &p,Vector x,Vector n);

	// Runs the method until optimum is found.
	void run();
};*/

class MinBrakRoutine {
private:
	double ulim,u,r,q,fu,dum;
	
	Vector vec;
	Function *problem;
	
	bool done;
	double func(double value);
public:
	// These are the variables we care about
	double ax, bx, cx, fa, fb, fc;

	MinBrakRoutine(Function *problem,double ax, double bx);
	void performIteration();
	bool isFinished();

	double getMinimum();
	double getFAtMinimum();
};

inline double MinBrakRoutine::getFAtMinimum() {
	return fb < fc ? fb : fc;
}

inline double MinBrakRoutine::getMinimum() {
	return fb < fc ? bx : cx;
}

/**
 * An abstract line minimization routine
 */
class LineMinRoutine {
protected:
	Function *problem;
	Vector vec;
	double func(double value);
public:
	LineMinRoutine(Function *problem);

	virtual void performIteration() = 0;
	virtual bool isFinished() = 0;

	// Get the result
	virtual double getMinimum() = 0;
	virtual double getFAtMinimum() = 0;
};


// Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
// between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
// the minimum to a fractional precision of about tol using Brent's method. The abscissa of
// the minimum is returned as xmin, and the minimum function value is returned as brent, the
// returned function value.
class BrentRoutine : public LineMinRoutine {
private:
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double ax, bx, cx, tol, xmin, fxmin;
	double e;								

	bool done;
public:
	BrentRoutine(Function *problem,double ax, double bx, double cx,double tol);
	void performIteration();
	bool isFinished();

	// Get the result
	double getMinimum();
	double getFAtMinimum();

	friend class BrentLinearMethod;
};


class GoldenRoutine : public LineMinRoutine {
private:
	double f1,f2,x0,x1,x2,x3;
	double ax, bx, cx, tol;

	static const double R,C;
	bool done;
public:
	GoldenRoutine(Function *problem,double ax, double bx, double cx,double tol);
	void performIteration();
	bool isFinished();

	// Get the result
	double getMinimum();
	double getFAtMinimum();

	friend class BrentLinearMethod;
};



/**
 * A one dimensional problem wrapped around an n-dimensional problem in a vector direction
 */
class Directional1DFunction : public Function {
private:
	Vector X,N,X1;
	Function *fnd;
public:
	Directional1DFunction(Function *fnd,const Vector &X,const Vector &N);
	
	double evaluate(const Vector &V);
	friend class BrentLinearMethod;
};

class BrentLinearMethod {
protected:
	Directional1DFunction f1d;
	
	MinBrakRoutine *mnbrak;
	
	// BrentRoutine *brent;
	LineMinRoutine *linmin;

	bool done;

public:
	// Constructor (does not require a solution space because there is a single
	// defined starting point.
	BrentLinearMethod(Function *problem,const Vector &x,const Vector &n);
	~BrentLinearMethod();

	// Perform an iteration
	void performIteration();
	
	// Runs the method until optimum is found.
	bool isFinished();

	// Get the minimum if the method is finished
	Vector getMinimum();
	Vector getMinDirection();
	double getFAtMinimum();
};

/*******************************************************************
 Inline methods
*******************************************************************/
inline bool BrentLinearMethod::isFinished() {
	return done;
}



/*******************************************************************
 Solve a problem of finiding s that optimizes f(P+sn) This uses derivatives
 to speed up the process.

 Use this by creating a BrentDerivativeMethod and calling run.  Check
 the values of x,n and 'value' after the algorithm has finished.
 ******************************************************************/
class BrentDerivativeMethod {
protected:
	DifferentiableFunction &p;

	double func(double t) {		
		xt = n;
		xt *= t;
		xt += x;
		return p.evaluate(xt);
	}

	double dfunc(double t,double &fn) {		
		xt = n;
		xt *= t;
		xt += x;
		return p.computeDirectionalJet(xt,n,fn);
	}

	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc);
	double dbrent(double ax, double bx, double cx, double tol,double *xmin);

public:
	// Position of the solution
	Vector x;

	// Vector from old solution to new solution
	Vector n;
	Vector xt;

	// Optimal valur
	double value;

	// Constructor (does not require a solution space because there is a single
	// defined starting point.
	BrentDerivativeMethod(DifferentiableFunction &p,Vector x,Vector n);

	// Runs the method until optimum is found.
	void run();
};


/*******************************************************************
 A one dimensional root finder ala Brent and NRC
 ******************************************************************/
class BrentRootFinder {
protected:
	Function1D &f;

	double zbrent(double x1,double x2,double tol);

public:
	// Constructor (does not require a solution space because there is a single
	// defined starting point.
	BrentRootFinder(Function1D &f);

	// Execute the method
	double findRoot(double x1,double x2,double tol) {
		return zbrent(x1,x2,tol);
	}

	// Possible return value
	const static double ROOT_ERROR;
};



// End namespace
NAMESPACE_PAULY_END

#endif //_OPTIMA_LINMIN_

