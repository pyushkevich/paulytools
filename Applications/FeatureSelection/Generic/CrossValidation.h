#ifndef __CrossValidation_h_
#define __CrossValidation_h_

#include "LinearSeparation.h"

class FeatureSetCrossValidation {
private:
	LPWrapper *lpWrapper;

public:
	FeatureSetCrossValidation();
	virtual ~FeatureSetCrossValidation();

	// Set the LP provider
	void setLPProvider(LPWrapper *provider);

	// Perform x-validation on a pair of classes
	double leaveSomeOut(const Mat &A, const Mat &B, const Vec &w, double pLeftOut, int nTries);
};

#endif


