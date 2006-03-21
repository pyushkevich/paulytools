#ifndef __OptimizationTerms_h_
#define __OptimizationTerms_h_

#include <smlmath.h>
#include "MedialPDESolver.h"
#include "BasisFunctions2D.h"
#include "ScriptInterface.h"
#include "CoefficientMask.h"
#include <optima.h>

using medialpde::FloatImage;

/** 
 * A base class for data associated with each solution as well as the 
 * data associated with the partial derivatives of solutions
 */
class SolutionDataBase
{
public:
  SolutionDataBase();
  virtual ~SolutionDataBase();
  
  // The grid on which atoms are sampled
  MedialAtomGrid *xAtomGrid;

  // The array of atoms
  MedialAtom *xAtoms;

  // Internal points of the m-rep
  SMLVec3d *xInternalPoints;

  // The array of boundary weights (area elements) and the total bnd. area
  double  *xBoundaryWeights, xBoundaryArea;

  // The array of internal weights (volume elts) and the total m-rep volume
  double  *xInternalWeights, xInternalVolume;

  // Weights associated with the chunks of profiles
  double *xInternalProfileWeights;

  // Whether the boundary/internal weights have been computed
  bool flagBoundaryWeights, flagInternalWeights, flagOwnAtoms;

  // Number of internal cuts at which the weights have been computed
  size_t nInternalCuts, nInternalPoints, nProfileIntervals;

  /** This method makes sure that the boundary area information is computed.
   * If the information has already been computed, it will not be computed
   * again. */
  virtual void UpdateBoundaryWeights() = 0;

  /** This method makes sure that the volume information is computed at the
   * resolution (nCuts = number of cuts in medial sail) given */
  virtual void UpdateInternalWeights(size_t nCuts) = 0;
};

/** 
 * A structure to be passed in and out of the energy term functions
 */
class SolutionData : public SolutionDataBase
{
public:

  // Initialize the solution data using a grid and an array of derivative
  // atoms
  SolutionData(MedialPDESolver *xSolver);
  virtual ~SolutionData() {}

  /** This method makes sure that the boundary area information is computed.
   * If the information has already been computed, it will not be computed
   * again. */
  void UpdateBoundaryWeights();

  /** This method makes sure that the volume information is computed at the
   * resolution (nCuts = number of cuts in medial sail) given */
  void UpdateInternalWeights(size_t nCuts);
  
  
};

class PartialDerivativeSolutionData : public SolutionDataBase
{
public:

  // Initialize the solution data using a grid
  PartialDerivativeSolutionData(SolutionData *xReference, MedialAtom *dAtoms);
  virtual ~PartialDerivativeSolutionData() {}

  /** This method makes sure that the boundary area information is computed.
   * If the information has already been computed, it will not be computed
   * again. */
  void UpdateBoundaryWeights();

  /** This method makes sure that the volume information is computed at the
   * resolution (nCuts = number of cuts in medial sail) given */
  void UpdateInternalWeights(size_t nCuts);

private:
  SolutionData *xReference;
};

class OffsetSolutionData : public SolutionDataBase
{
public:

  // Initialize the solution data using a grid
  OffsetSolutionData(
    SolutionData *sReference, 
    PartialDerivativeSolutionData *sDerivative, 
    double eps);
    
  virtual ~OffsetSolutionData() {}

  /** This method makes sure that the boundary area information is computed.
   * If the information has already been computed, it will not be computed
   * again. */
  void UpdateBoundaryWeights();

  /** This method makes sure that the volume information is computed at the
   * resolution (nCuts = number of cuts in medial sail) given */
  void UpdateInternalWeights(size_t nCuts);

private:
  SolutionData *xReference;
  PartialDerivativeSolutionData *xDerivative;
  double xEpsilon;
};



class EnergyTerm 
{
public:
  // Compute the energy term 
  virtual double ComputeEnergy(SolutionDataBase *data) = 0;

  // Initialize gradient computation and return the value of the solution
  // at the current state
  virtual double BeginGradientComputation(SolutionData *SCenter)
    { return ComputeEnergy(SCenter); }

  // Compute the partial derivative of the energy function
  virtual double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS) = 0;

  // Finish gradient computation, remove all temporary data
  virtual void EndGradientComputation() {};

  // Print a verbose report
  virtual void PrintReport(ostream &sout) = 0;
};

class NumericalGradientEnergyTerm : public EnergyTerm
{
public:
  NumericalGradientEnergyTerm()
    { xEpsilon = 0.0001; }
  
  void SetEpsilon(double eps) 
    { this->xEpsilon = eps; }
  
  double GetEpsilon()
    { return xEpsilon; }
  
  // Compute the finite difference gradient using epsilon
  virtual double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

private:
  double xEpsilon;
};

class BoundaryImageMatchTerm : public EnergyTerm
{
public:
  // Constructor
  BoundaryImageMatchTerm(FloatImage *image)
    { this->xImage = image; }

  // Compute the image match
  double ComputeEnergy(SolutionDataBase *data);

  // Initialize gradient computation and return the value of the solution
  // at the current state
  double BeginGradientComputation(SolutionData *SCenter);
  
  // Compute the partial derivative (must be called in the middle of Begin and
  // End of GradientComputation.
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Finish gradient computation, remove all temporary data
  void EndGradientComputation();

  // Print a verbose report
  void PrintReport(ostream &sout);

private:
  FloatImage *xImage;

  // Terms used in reporting details
  double xImageMatch, xBoundaryArea, xFinalMatch;

  // Temporaries for gradient computation
  SMLVec3d *xGradI;
  double *xImageVal;
};

/**
 * This term computes the amount of volume overlap between an image and the
 * m-rep. The image represents the object of interest with voxels of intensity
 * 1.0 and background with intensity 0.0. The image should have floating point
 * voxels, where values between 0 and 1 represent partial volume or
 * uncertainty
 */
class ProbabilisticEnergyTerm : public EnergyTerm
{
public:
  /** Initialize with an image and a number of sample points on each
   * medial sail vector */
  ProbabilisticEnergyTerm(FloatImage *image, size_t nCuts = 10);
  
  /** Compute the volume overlap fraction between image and m-rep */
  double ComputeEnergy(SolutionDataBase *data);

  // Initialize gradient computation and return the value of the solution
  // at the current state
  double BeginGradientComputation(SolutionData *SCenter);
  
  // Compute the partial derivative (must be called in the middle of Begin and
  // End of GradientComputation.
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Finish gradient computation, remove all temporary data
  void EndGradientComputation();

  // Print a verbose report
  void PrintReport(ostream &sout);
  
private:
  size_t nCuts;
  FloatImage *xImage;

  // Integral of object probability over the image
  double xImageIntegral;

  // Image values at the boundary
  double *xImageVal;

  // Scaled normal vectors along the boundary
  SMLVec3d *xBndNrm;

  // Cached Terms for reporting match details
  double xObjectIntegral, xRatio;
};

class BoundaryJacobianEnergyTerm : public EnergyTerm
{
public:
  // Compute the energy
  double ComputeEnergy(SolutionDataBase *data);

  // Print a verbose report
  void PrintReport(ostream &sout);

  // Compute the partial derivative term
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Get the worst Jacobian
  double GetMinJacobian() { return xMinJacobian; }
  
private:
  // Maximum and minimum values of the Jacobian encountered during
  // optimization
  double xMinJacobian, xMaxJacobian, xAvgJacobian, xTotalPenalty;

  // In addition to the total values, we keep track of per-quad values
  // computed in the course of calculating the penalty. This allows us
  // to compute derivatives more quickly
  struct QuadEntry {
    SMLVec3d XU, XV, YU[2], YV[2], NX, NY[2];
    double gX2, J[2], PenA[2], PenB[2];
  };
  
  typedef vector<QuadEntry> QuadVector;
  QuadVector xQuadEntries;

  /** Penalty function applied to the squared jacobian */
  double PenaltyFunction(double x, double a, double b)
    { return exp(-a * x) + exp(x - b); }

  /** The derivative of the penalty function */
  double PenaltyFunctionDerivative(double x, double a, double b)
    { return exp(x - b) - a * exp (-a * x); }
};


class MedialAnglesPenaltyTerm : public EnergyTerm
{
public:
  // Compute the energy
  double ComputeEnergy(SolutionDataBase *data);

  // Print a verbose report
  void PrintReport(ostream &sout);

  // Compute the partial derivative term
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);
  
private:
  double xTotalPenalty;
};

class CrestLaplacianEnergyTerm : public NumericalGradientEnergyTerm
{
public:
  double ComputeEnergy(SolutionDataBase *data);
  void PrintReport(ostream &sout);
private:
  double xMaxLaplacian, xAvgLaplacian, xTotalPenalty;
  size_t nCrestAtoms, nBadSites;
  double PenaltyFunction(double x, double a, double b)
    { return exp( a * x - b ); }
};

class AtomBadnessTerm : public EnergyTerm
{
public:
  double ComputeEnergy(SolutionDataBase *data);
  void PrintReport(ostream &sout);

  // Compute the partial derivative
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Get the number of bad atoms
  size_t GetNumberOfBadAtoms() 
    { return nBadAtoms; }
  
private:
  double xMaxBadness, xAvgBadness, xTotalPenalty;
  size_t nBadAtoms, nAtoms;
};

/**
 * Penalty for small values of the radius (phi / scale)^-8
 */
class RadiusPenaltyTerm : public EnergyTerm
{
public:
  RadiusPenaltyTerm(double xScale)
    { this->xScale = xScale; }

  double ComputeEnergy(SolutionDataBase *data);
  
  void PrintReport(ostream &sout);

  // Compute the partial derivative
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);
private:
  double xMinR2, xTotalPenalty, xScale;
};

class MedialRegularityTerm : public NumericalGradientEnergyTerm
{
public:
  /** Initialize with a sample solution */
  MedialRegularityTerm(MedialAtomGrid *xGrid, MedialAtom *xAtoms);
  
  /** Initialize the term with the template info */
  double ComputeEnergy(SolutionDataBase *data);

  // Compute partial derivative
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);
  
  void PrintReport(ostream &sout);
private:
  double xGradMagIntegral, xGradMagMaximum, xDomainArea;
  vnl_vector<double> xDomainWeights;
};

class MedialOptimizationProblem : public DifferentiableFunction
{
public:
  /** Initialize the problem */
  MedialOptimizationProblem(
    MedialPDESolver *xSolver, IMedialCoefficientMask *xCoeff);

  ~MedialOptimizationProblem();

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

  /** Compute the function and the partial derivative in one given direction */
  // double ComputePartialDerivative(double *xEvalPoint, double *xDirection, double &df);

  /** Print a detailed report */
  void PrintReport(ostream &sout);

private:
  // Precision at which the problem is solved
  static const double xPrecision;

  // The epsilon used to compute finite difference derivatives
  static const double xEpsilon;

  // The number of coefficients
  size_t nCoeff;

  // Value of the last match
  double xLastSolutionValue;

  MedialPDESolver *xSolver;
  IMedialCoefficientMask *xCoeff;
  vector<double> xWeights;
  vector<EnergyTerm *> xTerms;
  vector<CodeTimer> xTimers, xGradTimers;
  CodeTimer xSolveTimer, xSolveGradTimer, xWeightsTimer, xWeightsGradTimer;

  // Whether the gradient is available
  bool flagLastEvalAvailable, flagPhiGuessAvailable;
  
  // Last place where the function was evaluated
  vnl_vector<double> xLastEvalPoint;

  // The solution at the last evaluation point
  vnl_matrix<double> xLastPhiField;

  // vnl_vector<double> xLastGradEvalPoint, xLastGrad;

  // The phi / dPhi fields at the last evaluation point
  // vnl_matrix<double> xLastPhiField, xLastPhiDerivField; 

  // This method solves the MedialPDE, potentially using the last gradient
  // evaluation as the guess
  bool SolvePDE(double *xEvalPoint);

  // Legacy central difference solver
  void ComputeCentralDifferenceGradientPhi(double *x);

  // The array of derivative atoms
  // MedialAtom *dAtoms;
  vector<MedialAtom *> dAtomArray;
};

#endif
