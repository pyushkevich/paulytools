#ifndef _OPTIMIZATION_H_
#define _OPTIMIZATION_H_

#include <optima.h>
#include <ConjugateGradientMethod.h>
#include <EvolutionaryStrategy.h>
#include <Powell.h>
#include <registry.h>
#include "mspline.h"
#include "imaging.h"
#include "imatch.h"

/**
 * A generic prior or constraint	
 */
class SplinePrior {
protected:
	MSpline *spline;

public:
	// Complete prior probability
	virtual float priorProbabilityFull(SplineDataCache *sdc);

	// Prior probability for a region of control points
	virtual float priorProbability(SplineDataCache *sdc,int cuFirst,int cvFirst,int cuLast,int cvLast);

	// Constructor
	SplinePrior(MSpline *spline) {
		this->spline = spline;
	}
};

/**
 * A patch by patch prior
 */
class SplinePatchPrior : public SplinePrior {
protected:
	// Patch arrays
	Array2D<long> paTimeStamp;
	Array2D<float> paMeasure;
	Array2D<float> paWeight;

	// Compute a prior along a patch
	virtual float priorOverPatch(SplineDataCache *sdc,int iPatch,int jPatch,float &weight) = 0;

public:

	// Prior probability for a region of control points
	virtual float priorProbabilityFull(SplineDataCache *sdc);

	// Constructor
	SplinePatchPrior(MSpline *spline);
};

/**
 * A measure integrated over area of the patch
 */

/**
 * This is a constraint that makes sure that all control points are in range and that
 */
class ControlPointConstraint : public SplinePrior {
public:
	ControlPointConstraint(MSpline *spline) : SplinePrior(spline) {}

	// Prior probability for a region of control points
	float priorProbability(SplineDataCache *sdc,int cuFirst,int cvFirst,int cuLast,int cvLast);
};

/**
 * Compute a markov prior on the deformation of the model relative to the starting point 
 */
class MarkovControlPrior : public SplinePrior {
public:
	// Initialization
	MarkovControlPrior(MSpline *spline);

	// Prior probability for a region of control points
	float priorProbabilityFull(SplineDataCache *sdc);

private:
	// Penalty for a pair of points
	float markovPairPenalty(int i1,int j1,int i2,int j2);

	// Control point values
	Array2D<SMLVec4f> ctlBase;
};

/**
 * Compute a markov prior on the deformation of the model relative to the starting point 
 */
class MarkovCurvaturePrior : public SplinePrior {
public:
	// Initialization
	MarkovCurvaturePrior(MSpline *spline) : SplinePrior(spline) {}

	// Prior probability for a region of control points
	float priorProbabilityFull(SplineDataCache *sdc);
};



/**
 * Compute a curvature integral prior on the spline
 */
class CurvaturePrior : public SplinePatchPrior, public MedialMeasure {
public:
	// Prior probability for a region
	virtual float priorOverPatch(SplineDataCache *sdc,int iPatch,int jPatch,float &weight);

	// Medial measure stuff
	float computeMedialMeasure(const MedialPoint &mp);

	CurvaturePrior(MSpline *spline) : SplinePatchPrior(spline) {}
};

/**
 * Medial Representation Constraint
 */
class MRepConstraint : public SplinePatchPrior, public BoundaryMeasure {
private:
	// The epsili associated with each penalty type.  Penalty is applied beginning
	// at epsilon distance from the boundary
	float epsRadius, epsBlum;

	// A registry object
	Registry *settings;
public:
	MRepConstraint(Registry *settings,MSpline *spline);
	MRepConstraint(MSpline *spline);

	// Prior probability for a region
	virtual float priorOverPatch(SplineDataCache *sdc,int iPatch,int jPatch,float &weight);

	// Medial measure stuff
	float computeBoundaryMeasure(const MedialPoint &mp,int side);
	float computeCrestBoundaryMeasure(const MedialPoint &mp);

};


/**
 * Medial Representation Constraint
 */
class MRepFastConstraint : public SplinePatchPrior {
public:
	MRepFastConstraint(MSpline *spline) : SplinePatchPrior(spline) {};

	// Prior probability for a region
	virtual float priorOverPatch(SplineDataCache *sdc,int iPatch,int jPatch,float &weight);
};


/**
 * Medial Representation Constraint
 */
class MRepRadiusConstraint : public SplinePatchPrior {
public:
	MRepRadiusConstraint(MSpline *spline) : SplinePatchPrior(spline) {};

	// Prior probability for a region
	virtual float priorOverPatch(SplineDataCache *sdc,int iPatch,int jPatch,float &weight);
};

/**
 * Optimization code for m-reps
 */ 
class SplineOptProblem : public Function {
public:
	// Constructor.  Takes a range of control points to optimize over.
	SplineOptProblem(MSpline *spline, SplineDataCache *cache, SplineImageMatcher *match, int cpFirstU, int cpFirstV, int cpRangeU, int cpRangeV,Registry *stg);

	// Evaluate a solution
	double evaluate(const Vector &v);

	// Apply a vector to the spline
	void applyVector(const Vector &v);

	// Get the size of the vector used to drive the problem
	int getVectorSize();

	// Get the vector for current spline state.  Vector should be properly sized
	void makeVector(Vector &v);

	// Get the standard deviation vector
	void makeSDVector(Vector &v);

private:
	SplineDataCache *splineData;
	MSpline *spline;
	SplineImageMatcher *match;

	int cpFirstU, cpFirstV, cpRangeU, cpRangeV, cpLastU, cpLastV;
	int vecSize;
	int uPatchFirst,uPatchLast,vPatchFirst,vPatchLast;

	// Backup vector
	Vector backup;

	// Priors and weights
	// CurvaturePrior priorCrv;
	// MRepConstraint constrMRep;
	// ControlPointConstraint constrControl;
	MRepFastConstraint priorFast;
	MRepRadiusConstraint priorRad;
	MarkovControlPrior priorCtl;
	MarkovCurvaturePrior priorCrv;

	float wgtFast,wgtRadius,wgtPriorCtl,wgtConstrMRep,wgtPriorCrv;
};

class SplineRigidProblem : public Function {
public:
	// Constructor.  Takes a range of control points to optimize over.
	SplineRigidProblem(MSpline *spline, SplineDataCache *cache, SplineImageMatcher *match);

	// Evaluate a solution
	double evaluate(const Vector &v);

	// Apply a vector to the spline
	void applyVector(const Vector &v);

	// Get the size of the vector used to drive the problem
	int getVectorSize();

	// Get the vector for current spline state.  Vector should be properly sized
	void makeVector(Vector &v);

	// Get the standard deviation vector
	void makeSDVector(Vector &v);

private:
	SplineDataCache *splineData;
	MSpline *spline;
	SplineImageMatcher *match;

	// Backup for spline data
	Array2D<SMLVec4f> B;

	// The center of the spline
	SMLVec3f C;

	// Revert the spline
	void revertSpline();
};

/**
 * Optimization driver for spline optimization
 */
class SplineOptDriver {
private:
	NumericalMethod *method;
	SplineOptProblem *problem;

	NumericalFunction *numFunc;
	SolutionSpace *ss;

	SplineDataCache *splineData;
	PatchDataCache *pdc;
	MSpline *spline;
	SplineImageMatcher *match;

	// Configure appropriate method
	void initConjGrad(Registry *settings);
	void initPowell(Registry *settings);
	void initEvolution(Registry *settings);

public:
	SplineOptDriver(MSpline *spline, SplineDataCache *cache, SplineImageMatcher *match, 
				    int cpFirstU, int cpFirstV, int cpRangeU, int cpRangeV, Registry *rSettings);
	~SplineOptDriver();

	// Perform optimization for ms milliseconds or until finished
	bool optimize(int ms);

	// Are we finished?
	bool isFinished();

	// Get the best value
	double getBestValue();

	// Apply the best solution
	void applyBestSolution();

	// Get the cost so far
	double getEvalCost();
};





#endif

