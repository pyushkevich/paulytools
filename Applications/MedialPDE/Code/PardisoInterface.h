#ifndef __PardisoInterface_h_
#define __PardisoInterface_h_

#include <iostream>

class UnsymmetricRealPARDISO
{
public:
  // Initialize the solver 
  UnsymmetricRealPARDISO();
  
  // Release all memory used by the solver
  ~UnsymmetricRealPARDISO();
   
  // Factor the system for arbitrary right hand sides and matrices of the same
  // non-zer element structure
  void SymbolicFactorization(size_t n, int *idxRows, int *idxCols, double *xMatrix);

  // Factor the system for a specific matrix, but arbitrary right hand side
  void NumericFactorization(double *xMatrix);

  // Solve the system for the given right hand side, solution in xSoln
  void Solve(double *xRhs, double *xSoln);

private:
  /** Internal data for PARDISO */
  size_t PT[64];
  int MTYPE;
  int IPARM[64];

  // Storage for data in intermediate steps
  int n, *idxRows, *idxCols;
  double *xMatrix;
  bool flagPardisoCalled;
};

#endif 
