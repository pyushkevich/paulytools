
#ifndef _FunctionRep_H
#define _FunctionRep_H

class FunctionRep {
 public:
  FunctionRep() {}
  ~FunctionRep() {}
  virtual double get(const double x) const = 0;
  virtual double get1stDeriv(const double x) const = 0;
  virtual double get2ndDeriv(const double x) const = 0;
};

#endif
