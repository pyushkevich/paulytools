#include "SumBSplineCubUni.h"

SumBSplineCubUni::SumBSplineCubUni (const double coeff1[], const int np1,const double coeff2[], const int np2) : b1(coeff1,np1), BSplineCubUni(coeff2,np2) {
}
	
double SumBSplineCubUni::get (const double x) const {
  double tmp = b1.get(x);
  tmp += BSplineCubUni::get(x); 
  return tmp;
}

double SumBSplineCubUni::get1stDeriv (const double x) const {
  double tmp = b1.get1stDeriv(x);
  tmp += BSplineCubUni::get1stDeriv(x); 
  return tmp;
}
double SumBSplineCubUni::get2ndDeriv (const double x) const {
 double tmp = b1.get2ndDeriv(x);
  tmp += BSplineCubUni::get2ndDeriv(x); 
  return tmp;
}


