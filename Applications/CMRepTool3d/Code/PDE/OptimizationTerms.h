#ifndef __OptimizationTerms_h_
#define __OptimizationTerms_h_

#include <smlmath.h>
#include "MedialPDESolver.h"
#include "BasisFunctions2D.h"
#include "ScriptInterface.h"
#include <optima.h>

using medialpde::FloatImage;

/** 
 * A structure to be passed in and out of the energy term functions
 */
struct SolutionData
{
  // Initialize the solution data using a grid
  SolutionData(MedialPDESolver *xSolver, bool flagCopyAtoms);
  ~SolutionData();

  MedialAtom *xAtoms;
  double *xCoeffArray;
  double *xBoundaryWeights, *xMedialWeights, xBoundaryArea;
  MedialAtomGrid *xAtomGrid;
  bool flagOwnAtoms;
};

class EnergyTerm 
{
public:
  // Compute the energy term 
  virtual double ComputeEnergy(SolutionData *data) = 0;

  // Compute the gradient of the energy term given the medial atoms
  virtual double ComputeGradient(SolutionData *S0, SolutionData **SGrad, 
    double xEpsilon, vnl_vector<double> &xOutGradient);

  // Print a verbose report
  virtual void PrintReport(ostream &sout, SolutionData *data) = 0;
};

class BoundaryImageMatchTerm : public EnergyTerm
{
public:
  // Constructor
  BoundaryImageMatchTerm(FloatImage *image)
    { this->xImage = image; }

  // Compute the image match
  double ComputeEnergy(SolutionData *data);

  // Compute the directional derivative of the image match function
  virtual double ComputeGradient(SolutionData *S0, SolutionData **SGrad, 
    double xEpsilon, vnl_vector<double> &xOutGradient);

  // Print a verbose report
  void PrintReport(ostream &sout, SolutionData *data);

private:
  FloatImage *xImage;
};

class BoundaryJacobianEnergyTerm : public EnergyTerm
{
public:
  // Compute the energy
  double ComputeEnergy(SolutionData *data);

private:
  // Maximum and minimum values of the Jacobian encountered during
  // optimization
  double xMinJacobian, xMaxJacobian;

  /** Penalty function applied to the squared jacobian */
  double PenaltyFunction(double x, double a, double b)
    { return exp( - a * x ) + exp( x - b ); }

  // Print a verbose report
  void PrintReport(ostream &sout, SolutionData *data);
};

class MedialOptimizationProblem : public DifferentiableFunction
{
public:
  /** Initialize the problem */
  MedialOptimizationProblem(MedialPDESolver *xSolver, IBasisRepresentation2D *xSurface)
    {
    this->xSurface = xSurface;
    this->xSolver = xSolver;
    }

  /** Add an image match term with a perscribed weight */
  void AddEnergyTerm(EnergyTerm *term, double xWeight);

  /** Evaluate the function (comes from optima library) */
  double evaluate(const Vector &X)
    { return Evaluate(X.getDataArray()); }

  /** Compute the function and the gradient (comes from optima library) */
  double computeOneJet(const Vector &X, Vector &XGrad)
    { return ComputeGradient( X.getDataArray(), XGrad.getDataArray() ); }

  /** Evaluate the function */
  double Evaluate(double *X);
  
  /** Compute the function and the gradient (comes from optima library) */
  double ComputeGradient(double *X, double *XGrad);

private:
  // Precision at which the problem is solved
  static const double xPrecision;

  // The epsilon used to compute finite difference derivatives
  static const double xEpsilon;

  MedialPDESolver *xSolver;
  IBasisRepresentation2D *xSurface;
  vector<double> xWeights;
  vector<EnergyTerm *> xTerms;
};

#endif
