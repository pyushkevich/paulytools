#ifndef _LINEAR_SEPARATION_H_
#define _LINEAR_SEPARATION_H_

#include "AbstractLPClassification.h"

class LinearSeparation : public AbstractLPClassification {
public:
    LinearSeparation();
    ~LinearSeparation();

    // Initialize the solver by passing two matrices
    void initialize(const Mat &A, const Mat &B);

    // Run the solver
    bool run();

    // Return the solution
    double getSeparatingPlane(Vec &w);

private:
    // The objective function vector
    Vec c;
};

#endif //_LINEAR_SEPARATION_H_
