#include <SimplexMethod.h>
#include <EvolutionaryStrategy.h>
#include <ConjugateGradientMethod.h>

#ifndef REG_H_
#define REG_H_


class SplineNodeGroup {
public:
	SplineObject *spline;
	virtual int getVectorSize() = 0;
	virtual Vector getVector() = 0;
	virtual void applyVector(const Vector &v) = 0;
	virtual Vector getStepSize() = 0;

	virtual void updateSample(RegularSplineSample *sample) = 0;

	SplineNodeGroup(SplineObject *spline) {
		this->spline = spline;
	}
};

class SNGSegment : public SplineNodeGroup {
protected:
	int iCurve,iStart,iEnd;
	int vSize;
	vector<int> vmap;
	SplineCurve *curve;

	Vector vStep;

public:
	SNGSegment(SplineObject *spline,int iCurve,int iStart,int iEnd);
	
	virtual int getVectorSize() {
		return vSize;
	}

	virtual Vector getVector();
	virtual Vector getStepSize();
	virtual void applyVector(const Vector &v);
	virtual void updateSample(RegularSplineSample *sample);
};

class SNGBranch : public SplineNodeGroup {
protected:
	// int iCurve,iStart,iEnd;
	// int vSize;
	// vector<int> vmap;
	SplineCurve *curve;
	SplineBranch *branch;

public:
	SNGBranch(SplineObject *spline,SplineBranch *branch);
	
	int getVectorSize() {
		return 9;
	}

	Vector getVector();
	Vector getStepSize();
	void applyVector(const Vector &v);
	virtual void updateSample(RegularSplineSample *sample);
};

/***********************************************************
 Code to perform image match on a medial spline.  This problem 
 applies to k control points at a time
 **********************************************************/
class SplineMatchProblem : public Function {
private:
	// Image space to which computations apply
	ImageSpace *space;

	// The starting point and the number of control points to which
	// the problem applies
	SplineNodeGroup *sng;

	// A spline sample;
	RegularSplineSample *sample;

	// An image match object
	SplineSampleImageMatch *imatch;

	// Weight factors for different penaly terms
	double wgtPriorCurvature, wgtPriorElastic, wgtPenaltyJacobian, wgtPenaltyRadius, wgtPenaltyBranch;

public:
	SplineMatchProblem(SplineNodeGroup *sng,ImageSpace *space);
	~SplineMatchProblem();

	// Get a solution space
	GaussianSS getSolutionSpace();

	// Evaluate a vector
	double evaluate(const Vector &v);

	// Return the size of the vector used to control this problem
	int getSize() {
		return sng->getVectorSize();
	}

	// This method is called to update the spline with a vector value
	void updateSpline(const Vector &v) {
		sng->applyVector(v);
		sng->updateSample(sample);
	}

	// Return the vector for the current state of the spline
	Vector getVector() {
		return sng->getVector();
	}

	// Compute constraints on the spline
	double computeConstraints();

	// Compute the regularization prior
	double computeRegularizationPrior();
};





// A nice method that creates a right optimizer based on global settings
NumericalMethod *createOptAlg(NumericalFunction &problem,SolutionSpace &ss);

// Tools for rigid registration
void runRegistration();
void runNodeOptimization();
void runMultiNodeOpt();

// Tool for figure interpolation
void interpolateFigure(MFigure *figure,int n);
void resampleFigureUniform(MFigure *figure);

// Optimization of splines
void runSplineMatch(SplineObject *mrep, ImageSpace *space);


#endif

