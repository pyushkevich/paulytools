#include "FeatureSelection.h"


FeatureSelection::FeatureSelection()
    : WindowSelection() 
{
}

FeatureSelection::~FeatureSelection() 
{
}

void FeatureSelection::initialize(const Mat &A, const Mat &B) {
    Mat O(A.columns(),A.columns());
    O = 1;
    WindowSelection::initialize(A,B,O);
}

void FeatureSelection::initialize(const Mat &A, const Mat &B, const Mat &Omega) {
    initialize(A,B);
}


