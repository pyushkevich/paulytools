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
  bool flagBoundaryWeights, flagInternalWeights;

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
  ~SolutionData();

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
  ~PartialDerivativeSolutionData();

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

class EnergyTerm 
{
public:
  // Compute the energy term 
  virtual double ComputeEnergy(SolutionData *data) = 0;

  // Initialize gradient computation and return the value of the solution
  // at the current state
  virtual double BeginGradientComputation(SolutionData *SCenter);
  
  // Compute the partial derivative (must be called in the middle of Begin and
  // End of GradientComputation.
  virtual double ComputePartialDerivative(
    SolutionData *SCenter, SolutionData *SFwd, SolutionData *SBck,
    double xEpsilon); 

  // Compute the partial derivative given the partial derivative of the atoms
  // virtual double ComputePartialDerivative(
  //  SolutionData *S, SolutionData *dS) = 0;

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
  double ComputePartialDerivative(
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
  double ComputeEnergy(SolutionData *data);
/*
  // Initialize gradient computation and return the value of the solution
  // at the current state
  double BeginGradientComputation(SolutionData *SCenter);
  
  // Compute the partial derivative (must be called in the middle of Begin and
  // End of GradientComputation.
  double ComputePartialDerivative(
    SolutionData *SCenter, SolutionData *SFwd, SolutionData *SBck,
    double xEpsilon);

  // Finish gradient computation, remove all temporary data
  void EndGradientComputation();
*/
  // Print a verbose report
  void PrintReport(ostream &sout);

private:
  FloatImage *xImage;

  // Terms used in reporting details
  double xImageMatch, xBoundaryArea, xFinalMatch;

  // Temporaries for gradient computation
  SMLVec3d *xGradI;
  double *xImageVal, xMatch, xArea;
};

/**
 * This term computes the amount of volume overlap between an image and the
 * m-rep. The image represents the object of interest with voxels of intensity
 * 1.0 and background with intensity 0.0. The image should have floating point
 * voxels, where values between 0 and 1 represent partial volume or
 * uncertainty
 */
class VolumeOverlapEnergyTerm : public EnergyTerm
{
public:
  /** Initialize with an image and a number of sample points on each
   * medial sail vector */
  VolumeOverlapEnergyTerm(FloatImage *image, size_t nCuts = 10);
  
  /** Compute the volume overlap fraction between image and m-rep */
  double ComputeEnergy(SolutionData *data);
/*
  // Initialize gradient computation and return the value of the solution
  // at the current state
  double BeginGradientComputation(SolutionData *SCenter);
  
  // Compute the partial derivative (must be called in the middle of Begin and
  // End of GradientComputation.
  double ComputePartialDerivative(
    SolutionData *SCenter, SolutionData *SFwd, SolutionData *SBck,
    double xEpsilon);

  // Finish gradient computation, remove all temporary data
  void EndGradientComputation();
*/
  // Print a verbose report
  void PrintReport(ostream &sout);
  
private:
  size_t nCuts;
  FloatImage *xImage;
  double xImageVolume;

  // Temporary for gradient computation
  double *xImageVal;

  // Cached Terms for reporting match details
  double xIntersect, xObjectVolume, xOverlap;
};

class BoundaryJacobianEnergyTerm : public EnergyTerm
{
public:
  // Compute the energy
  double ComputeEnergy(SolutionData *data);

  // Print a verbose report
  void PrintReport(ostream &sout);

private:
  // Maximum and minimum values of the Jacobian encountered during
  // optimization
  double xMinJacobian, xMaxJacobian, xAvgJacobian, xTotalPenalty;

  /** Penalty function applied to the squared jacobian */
  double PenaltyFunction(double x, double a, double b)
    { return exp( - a * x ) + exp( x - b ); }
};

class CrestLaplacianEnergyTerm : public EnergyTerm
{
public:
  double ComputeEnergy(SolutionData *data);
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
  double ComputeEnergy(SolutionData *data);
  void PrintReport(ostream &sout);
private:
  double xMaxBadness, xAvgBadness, xTotalPenalty;
  size_t nBadAtoms, nAtoms;
};

class MedialRegularityTerm : public EnergyTerm
{
public:
  /** Initialize the term with the template info */
  MedialRegularityTerm(MedialAtom *xTempAtoms, MedialAtomGrid *xTempGrid);
  double ComputeEnergy(SolutionData *data);
  void PrintReport(ostream &sout);
private:
  double ComputeDistortionPenalty(double, double);
  vector<double> xEdgeLength;
  double xMaxDistortion, xMinDistortion, xMeanSquareDistortion, xTotalPenalty;
};

class MedialOptimizationProblem : public DifferentiableFunction
{
public:
  /** Initialize the problem */
  MedialOptimizationProblem(MedialPDESolver *xSolver, IMedialCoefficientMask *xCoeff)
    {
    this->xCoeff = xCoeff;
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

  /** Print a detailed report */
  void PrintReport(ostream &sout);

private:
  // Precision at which the problem is solved
  static const double xPrecision;

  // The epsilon used to compute finite difference derivatives
  static const double xEpsilon;

  // Value of the last match
  double xSolution;

  MedialPDESolver *xSolver;
  IMedialCoefficientMask *xCoeff;
  vector<double> xWeights;
  vector<EnergyTerm *> xTerms;
};

#endif
