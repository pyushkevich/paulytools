#ifndef _ABSTRACT_LP_CLASSIFICATION_H_
#define _ABSTRACT_LP_CLASSIFICATION_H_

#include "LPWrapper.h"

class AbstractLPClassification {
public:
	AbstractLPClassification() {
		state = UNINIT;
	}

	virtual ~AbstractLPClassification() {
	}

	// Set the linear programming provider
	virtual void setLPProvider(LPWrapper *provider) {
		m_provider = provider;
		state = UNINIT;
	}

	// Return the separating plane and its cross-section
	double getSeparatingPlane(Vec &w) {
		assert(state == SUCCESS || state == FAILURE);
		w = wResult;
		return gResult;
	}

	// Get the objective function value
	double getObjectiveValue() {
		assert(state == SUCCESS || state == FAILURE);
		return objResult;
	}

	// Common run method
	virtual bool run() = 0;

protected:
	// Dimensions of A and B
	int m, k, n;

	// State of the process
	enum State {UNINIT,INIT,SUCCESS,FAILURE} state;

	// The LP provider
	LPWrapper *m_provider; 

	// The resulting vectors
	Vec wResult;
	double gResult,objResult;
};

#endif // _ABSTRACT_LP_CLASSIFICATION_H_

