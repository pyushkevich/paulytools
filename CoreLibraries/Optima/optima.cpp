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
 * optima.cpp
 *	----------
 * Header declares base classes in the library
 ******************************************************************/
#include <optima.h>

// Begin namespace
NAMESPACE_PAULY_START

/******************************************************************
 Gaussian Solution Space
 ******************************************************************/
// Create a Gaussian solution space
GaussianSS::GaussianSS(const Vector &mean,const Vector &stdev) {
	this->mean = mean;
	this->stdev = stdev;
}

// Get a feasible solution
Vector GaussianSS::getFeasibleSolution() {
	Vector x(mean.size());

   for(int i=0;i<mean.size();i++) {
      x(i) = getGaussianRnd(mean(i),stdev(i));
   }

	return x;
}

// Check the feasibility of the solution (feasibility 0 means that the solution may not 
// exist at a given position
double GaussianSS::getFeasibility(const Vector &x) {
   double feasibility = 1.0;
	for(int i=0;i<mean.size();i++) {
		feasibility *= penaltyFunction(x(i),mean(i),stdev(i));
	}
   return feasibility;
}

/******************************************************************
 Uniform Solution Space
 ******************************************************************/
// Create a uniform solution space
UniformSS::UniformSS(const Vector &center,const Vector &radius) {
	this->center = center;
	this->radius = radius;
}

// Get a feasible solution
Vector UniformSS::getFeasibleSolution() {
	Vector x(center.size());

   for(int i=0;i<center.size();i++) {
      x(i) = center(i) - radius(i) + rand(2*radius(i));
   }

	return x;
}

// Check the feasibility of the solution (feasibility 0 means that the solution may not 
// exist at a given position
double UniformSS::getFeasibility(const Vector &x) {
	for(int i=0;i<center.size();i++) {
		if(x(i) < center(i) - radius(i) || x(i) > center(i) + radius(i))
			return 0;
	}

	return 1.0;
}


/******************************************************************
 Jet approximators
 ******************************************************************/
double approxOneJet(Function &f,const Vector &x,Vector &dxi,const Vector &epsili,double gScale) {
	// The function value
	double fx = f.evaluate(x);

	// The gradient vector
	dxi.setSize(x.size());	
	
	// Compute each partial derivative
   Vector x1 = x;
	for(int i=0;i<x.size();i++) {
		x1(i) += epsili(i);
		double fX = f.evaluate(x1);
		x1(i) -= epsili(i);

		dxi(i) = gScale * (fX-fx) / epsili(i);
	}

	return fx;
}

double approxOneJet(Function &f,const Vector &x,Vector &dxi,double epsilon,double gScale) {
	// The function value
	double fx = f.evaluate(x);

	// The gradient vector
	dxi.setSize(x.size());	
	
	// Compute each partial derivative
   Vector x1 = x;
	for(int i=0;i<x.size();i++) {
		x1(i) += epsilon;
		double fX = f.evaluate(x1);
		x1(i) -= epsilon;

		dxi(i) = gScale * (fX-fx) / epsilon;
	}

	return fx;
}

NumericalFunction::NumericalFunction(Function &f,const Vector &epsili) : function(f) {
	this->epsili = epsili;
	epsilon = 0.0;
	gScale = 1.0;
}

NumericalFunction::NumericalFunction(Function &f,double epsilon) : function(f) {
	assert(epsilon > 0.0);
	this->epsilon = epsilon;
	gScale = 1.0;
}


double NumericalFunction::computeOneJet(const Vector &v,Vector &dxi) {
	if(epsilon > 0.0) {
		return approxOneJet(function,v,dxi,epsilon,gScale);
	}
	else {
		return approxOneJet(function,v,dxi,epsili,gScale);
	}
}

double NumericalFunction::computeDirectionalJet(const Vector &x,const Vector &u,double &fu) {
	double eps = epsilon;
	if(eps==0) {
		// Must compute a new epsilon, lets use epsili dot u
		eps = epsili.dotProduct(u);
	}
	
	// The function value
	double fx = evaluate(x);
	double fX = evaluate(x + eps*u);
	fu = (fX-fx) / eps;

	return fx;
}


/******************************************************************
 Solution stuff
 ******************************************************************/
Solution::Solution() {
	value = 0;
   hasValue = false;
}

Solution::Solution(const Solution &copy) {
	setValue(copy.value);
}

NumericalSolution::NumericalSolution(int dims) : Solution(),x(dims) {
	assert(dims > 0);
}

NumericalSolution::NumericalSolution(const NumericalSolution &ns) : Solution(ns),x(ns.x) {
}

NumericalSolution::NumericalSolution(const Vector &inX) : Solution(),x(inX) {
}

// End namespace
NAMESPACE_PAULY_END
