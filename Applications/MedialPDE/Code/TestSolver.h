#ifndef __TestSolver_h_
#define __TestSolver_h_

class MedialPDESolver;
class IMedialCoefficientMask;

int TestGradientComputation(
  MedialPDESolver *xSolver, IMedialCoefficientMask *xMask, int nRandVariations = 0);

int TestOptimizerGradientComputation(
  MedialOptimizationProblem &mop, 
  IMedialCoefficientMask &xMask,
  MedialPDESolver *xSolver);

#endif
