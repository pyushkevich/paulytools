#ifndef _LP_WRAPPER_H_
#define _LP_WRAPPER_H_

// This file includes the matrix and class definitions
#include "matrix.h"

// Shorthand for matrix and vector classes that doesn't clash with other code
typedef pauly::Matrix Mat;
typedef pauly::Vector Vec;

/**
 * An abstract class that represents generic functionality 
 * provided by an LP solver.  The problems are specified in the
 * following format:
 * 	find \f$\bf x$\f that minimizes \f${\bf c}^T {\bf x}$\f
 * 	subject to: \f$M {\bf x} \ge {\bf b}$\f,
 * 	and \f${\bf l} \le {\bf x} \le {\bf u}$\f
 */
class LPWrapper {
    public:
	// Constructor just initializes the lp solver
	LPWrapper();

	// Destructor cleans up
	~LPWrapper();

	/**
	 * Set up a problem by passing in the objective, and the 
	 * inequalities.
	 * @param c	Objective vector
	 * @param M	Matrix of inequalities
	 * @param b	left-hand-side vector
	 * @param u	upper bound
	 * @param l	lower bound
	 */
	virtual void setProblem(const Vec &c,const Mat &M,
		const Vec &b,const Vec &l, const Vec &u) = 0;

	/**
	 * Set up the problem with default upper and lower bounds 
	 */
	virtual void setProblem(const Vec &c,const Mat &M,const Vec &b);

	/**
	 * Update the objective <b>c</b> without (hopefully) 
	 * reinitializing the problem
	 */
	virtual void updateObjective(const Vec &c) = 0;

	/**
	 * Update the bounds without (hopefully)
	 * reinitializing the problem
	 */
	virtual void updateBounds(const Vec &u,const Vec &l) = 0;

	/**
	 * Update one of the rows in the inequalities without
	 * (hopefully) reinitializing the problem
	 */
	virtual void updateRow(int i,const Vec &row,double bVal) = 0;

	/**
	 * Solve the problem
	 * @param res 	Vector in which to place the solution
	 * @return	Whether or not the solver succeeded
	 */
	virtual bool solve(Vec & res) = 0;	
};

#endif // _LP_WRAPPER_H_
