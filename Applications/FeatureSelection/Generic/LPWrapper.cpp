#include "LPWrapper.h"
#include <limits>

LPWrapper::LPWrapper() {

}

LPWrapper::~LPWrapper() {

}

void LPWrapper::setProblem(const Vec &c,const Mat &M,const Vec &b)
{
    Vec upper(c.rows()),lower(c.rows());
    upper.setAll(1e100);
    lower.setAll(0.0);
    setProblem(c,M,b,lower,upper);
} 
