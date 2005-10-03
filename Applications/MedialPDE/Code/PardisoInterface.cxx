#include "PardisoInterface.h"
#include <iostream>

using namespace std;

// BLAS/PARDISO references
extern "C" {
  void pardisoinit_(size_t *, int *, int *);
  void pardiso_(size_t *, int *, int *, int *, int *, int *, double *, int *, int *, 
    int *, int *, int *, int *, double *, double *, int*);
}

UnsymmetricRealPARDISO::UnsymmetricRealPARDISO()
{
  // Set the type of matrix to unsymmetric real
  MTYPE = 11; 

  // Clear the parameter array
  memset(IPARM, 0, sizeof(int) * 64);

  // Initialize PARDISO to default values
  pardisoinit_(PT,&MTYPE,IPARM);

  // Specify the number of processors on the system (1)
  IPARM[2] = 1;

  flagPardisoCalled = false;
}

void 
UnsymmetricRealPARDISO
::SymbolicFactorization(size_t n, int *idxRows, int *idxCols, double *xMatrix)
{
  // Set the various parameters
  int MAXFCT = 1, MNUM = 1, PHASE = 11, N = n, NRHS = 1, MSGLVL = 0, ERROR = 0; 
  
  // Perform the symbolic factorization phase
  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, 
    xMatrix, idxRows, idxCols,
    NULL, &NRHS, IPARM, &MSGLVL, NULL, NULL, &ERROR);

  // Record the parameter for next phase
  this->idxCols = idxCols;
  this->idxRows = idxRows;
  this->n = n;

  // Set the flag so we know that pardiso was launched before
  flagPardisoCalled = true;
}

void 
UnsymmetricRealPARDISO
::NumericFactorization(double *xMatrix)
{
  // Set the various parameters
  int MAXFCT = 1, MNUM = 1, PHASE = 22, N = n, NRHS = 1, MSGLVL = 0, ERROR = 0; 
  
  // Perform the symbolic factorization phase
  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, 
    xMatrix, idxRows, idxCols,
    NULL, &NRHS, IPARM, &MSGLVL, NULL, NULL, &ERROR);

  // Record the parameter for next phase
  this->xMatrix = xMatrix;
}

void 
UnsymmetricRealPARDISO
::Solve(double *xRhs, double *xSoln)
{
  // Set the various parameters
  int MAXFCT = 1, MNUM = 1, PHASE = 33, N = n, NRHS = 1, MSGLVL = 0, ERROR = 0; 
  
  // Perform the symbolic factorization phase
  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, 
    xMatrix, idxRows, idxCols,
    NULL, &NRHS, IPARM, &MSGLVL, xRhs, xSoln, &ERROR);
}

UnsymmetricRealPARDISO::
~UnsymmetricRealPARDISO()
{
  // Set the various parameters
  int MAXFCT = 1, MNUM = 1, PHASE = -1, N = n, NRHS = 1, MSGLVL = 0, ERROR = 0; 
  
  // Perform the symbolic factorization phase
  if(flagPardisoCalled)
    pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, 
      xMatrix, idxRows, idxCols,
      NULL, &NRHS, IPARM, &MSGLVL, NULL, NULL, &ERROR);
}


