#ifndef _WINDOW_SELECTION_H_
#define _WINDOW_SELECTION_H_

#include "AbstractLPClassification.h"

class WindowSelection : public AbstractLPClassification {
    public:
	WindowSelection();
	~WindowSelection();
	
	/**
	 * Initialize by passing in class matrices A, B and the window
	 * matrix Omega 
	 */
	virtual void initialize(const Mat &A, const Mat &B, const Mat &Omega);

	// Set the modulating parameter on the number of features
	void setLinearModulation(double lambda, double eta);

	// Set the modulating parameter on the number of windows
	void setArbitraryModulation(const Vec &weights);

	// Set the alpha constant
	void setAlpha(double alpha) {
	    this->alpha = alpha;
	}

	// Execute the window selection
	virtual bool run();	

	// Returns the vector whose non-zer elements indicate the selection 
	// of windows in Omega
	void getWindowVector(Vec &out) {
	    assert(state == SUCCESS || state == FAILURE);
	    out = vResult;
	}
	
    private:
	// The number of windows in the system
	int m,n,k,nw;

	// The alpha and eps constants
	double alpha,eps;

	// The objective vector c
	Vec c;

	// The window length vector (sum of Omega's columns)
	Vec wl;

	// The weight vector that modulates v in the equation
	Vec wgt;

	// The resulting vectors
	Vec vResult;
};

#endif // _WINDOW_SELECTION_H_
