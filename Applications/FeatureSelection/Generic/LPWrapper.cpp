#include "LPWrapper.h"
#include <limits>

LPWrapper::LPWrapper() {

}

LPWrapper::~LPWrapper() {

}

void LPWrapper::setProblem(const Vec &c,const Mat &M,const Vec &b)
{
    Vec upper(c.size()),lower(c.size());
    upper.fill(1e100);
    lower.fill(0.0);
    setProblem(c,M,b,lower,upper);
} 
