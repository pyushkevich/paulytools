#ifndef __TestSolver_h_
#define __TestSolver_h_

class MedialPDESolver;
class IMedialCoefficientMask;

void TestGradientComputation(
  MedialPDESolver *xSolver, IMedialCoefficientMask *xMask);

void TestOptimizerGradientComputation(
  MedialOptimizationProblem &mop, 
  IMedialCoefficientMask &xMask,
  MedialPDESolver *xSolver);

#endif
