#ifndef _FEATURE_SELECTION_H_
#define _FEATURE_SELECTION_H_

#include "WindowSelection.h"

class FeatureSelection : public WindowSelection {
    public:
	FeatureSelection();
	~FeatureSelection();

	// A simpler initialize method
	void initialize(const Mat &A, const Mat &B);
	
	// Initialize is overridden to always take the identity matrix
	void initialize(const Mat &A, const Mat &B, const Mat &Omega);

	void setLambda(double lambda) {
	    setLinearModulation(lambda,0);
	}
};

#endif // _FEATURE_SELECTION_H_
