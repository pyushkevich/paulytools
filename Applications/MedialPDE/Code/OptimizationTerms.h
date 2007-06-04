#ifndef __OptimizationTerms_h_
#define __OptimizationTerms_h_

#include <smlmath.h>
#include "CodeTimer.h"
#include "GenericMedialModel.h"
#include "BasisFunctions2D.h"
#include "ScriptInterface.h"
#include "CoefficientMapping.h"
#include <optima.h>
#include "Registry.h"

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
  
  // The medial model that specifies the grid/mesh on which atoms live
  MedialIterationContext *xAtomGrid;

  // The array of atoms in this solution. This array can either point to the xMedialModel's
  // array or it can be pointing to another set of atoms
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
  SolutionData(MedialIterationContext *xGrid, MedialAtom *xAtoms);

  // Destructor
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

/** 
 * Accumulator for statistics in energy terms
 */
class StatisticsAccumulator
{
public:
  StatisticsAccumulator() { Reset(); }

  void Reset()
    { n = 0; xMin = xMax = xSum = xSumSq = 0.0; }

  void Update(double x)
    {
    // Update Min and Max
    if(n == 0)
      { xMin = xMax = x; }
    else
      { xMin = std::min(xMin, x); xMax = std::max(xMax, x); }

    // Update the rest
    xSum += x;
    xSumSq += x * x;
    n++;
    }

  double GetMin() const { return xMin; }
  double GetMax() const { return xMax; }
  double GetMaxAbs() const { return std::max(fabs(xMin), fabs(xMax)); }
  double GetMean() const { return xSum / n; }
  double GetVariance() const { return (xSumSq - xSum * xSum / n) / (n-1); }
  double GetStdDev() const { return sqrt(GetVariance()); }
  double GetSum() const { return xSum; }
  size_t GetCount() const { return n; }

private:
  size_t n;
  double xMin, xMax, xSum, xSumSq;
};

inline ostream & operator << (ostream &sout, const StatisticsAccumulator &S)
{
  sout << "mean(sd): " << S.GetMean() << "(" << S.GetStdDev() << "); ";
  sout << "range: " << S.GetMin() << " to " << S.GetMax() << "; ";
  sout << "n: " << S.GetCount() << ";";
  return sout;
}

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

/**
 * This is a special energy term that involves integrating stuff
 * over the medial axis, ie. \int f(u,v) du dv. For such a term, we
 * want the du dv quantity, which is non-uniform over the mesh in many
 * cases (because of non-uniform grid spacing)
 */
class MedialIntegrationEnergyTerm : public EnergyTerm
{
public:
  MedialIntegrationEnergyTerm(GenericMedialModel *model);

protected:
  double xDomainArea;
  vnl_vector<double> xDomainWeights;
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
  SMLVec3d *xImageGrad;

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
  double GetMinJacobian() { return saJacobian.GetMin(); }
  
private:
  // Maximum and minimum values of the Jacobian encountered during
  // optimization
  StatisticsAccumulator saJacobian, saPenalty;

  // double xMinJacobian, xMaxJacobian, xAvgJacobian, xTotalPenalty;

  // In addition to the total values, we keep track of per-quad values
  // computed in the course of calculating the penalty. This allows us
  // to compute derivatives more quickly
  struct TriangleEntry {
    SMLVec3d XU, XV, YU[2], YV[2], NX, NY[2];
    double gX2, J[2], PenA[2], PenB[2];
  };
  
  typedef vector<TriangleEntry> TriangleVector;
  TriangleVector xTriangleEntries;

  /** Penalty function applied to the squared jacobian */
  double PenaltyFunction(double x, double a, double b)
    { return exp(-a * x) + exp(x - b); }

  /** The derivative of the penalty function */
  double PenaltyFunctionDerivative(double x, double a, double b)
    { return exp(x - b) - a * exp (-a * x); }

  // Constants used in the penalty computation
  const static double xPenaltyA, xPenaltyB;
};


class MedialAnglesPenaltyTerm : public MedialIntegrationEnergyTerm
{
public:
  // Initialize the term
  MedialAnglesPenaltyTerm(GenericMedialModel *model);

  // Compute the energy
  double ComputeEnergy(SolutionDataBase *data);

  // Print a verbose report
  void PrintReport(ostream &sout);

  // Compute the partial derivative term
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);
  
private:
  double xTotalPenalty, xMaxPenalty;
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
  vnl_vector<double> xPenalty;
  double xMinBadness, xAvgBadness, xTotalPenalty;
  size_t nBadAtoms, nAtoms;
};

/**
 * Term that penalizes |GradR| <> 1 on the boundary
 */
class BoundaryGradRPenaltyTerm : public EnergyTerm
{
public:
  // Compute the penalty
  double ComputeEnergy(SolutionDataBase *data);

  // Describe the terms of the penalty
  void PrintReport(ostream &sout);

  // Compute the partial derivative
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

private:
  // Accumulators for display and statistics calculation
  StatisticsAccumulator saGradR, saPenalty;

  // Scaling factor used to derive the penalty (why?)
  const static double xScale;
};


/**
 * Term that penalizes excessive curvature of the medial axis. This
 * penalty has the form r^2 * (kappa1^2 + kappa2^2). Thus it penalizes
 * situations where one of the radii of curvature is excessively greater
 * than the radius of the medial model. Of course the model can avoid 
 * this penalty by pushing the radius to zero, so this penalty should
 * always be used in conjunction with a radius penalty.
 */
class MedialCurvaturePenalty : public EnergyTerm
{
public:
  // Compute the penalty
  double ComputeEnergy(SolutionDataBase *data);

  // Describe the terms of the penalty
  void PrintReport(ostream &sout);

  // Compute the partial derivative
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

private:
  // Accumulators for display and statistics calculation
  StatisticsAccumulator 
    saMeanCurv, saGaussCurv, saSumSqKappa, saRad, saFeature, saPenalty;

  // Parameters of the penalty, which is of the form pow(f/scale, power)
  const static double xScale, xPower;
};


/**
 * This term is used to fit phi to an existing radius field. The fitting
 * is least-squares in phi, which may not be a great idea, but if we make it
 * least-squares in R, there may be unforseen problems with stability at R=0
 */
class DistanceToRadiusFieldEnergyTerm : public MedialIntegrationEnergyTerm
{
public:
  /** Constructor, takes the target radius field */
  DistanceToRadiusFieldEnergyTerm(GenericMedialModel *model, double *radius);
  
  double ComputeEnergy(SolutionDataBase *data);
  
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  void PrintReport(ostream &sout);

private:

  typedef vnl_vector<double> Vec;
  Vec xTargetPhi, xLastDelta;
  double xTotalDiff, xMaxDiff, xMinDiff, xTotalMatch, xTotalArea;
};


/**
 * This term is used to fit the medial surface to an existing set of points 
 * given at u and v coordinates. This should be used in conjunction with a
 * regularization prior
 */
class DistanceToPointSetEnergyTerm : public MedialIntegrationEnergyTerm
{
public:
  /** Constructor, takes the target XYZ field */
  DistanceToPointSetEnergyTerm(
    GenericMedialModel *model, double *x, double *y, double *z);
  
  double ComputeEnergy(SolutionDataBase *data);
  
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  void PrintReport(ostream &sout);

private:

  // Target points
  vector<SMLVec3d> target;

  // Match components
  double xTotalDist,  xTotalMatch, xTotalArea;
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

class MedialBendingEnergyTerm : public MedialIntegrationEnergyTerm
{
public:
  /** Initialize with a sample solution */
  MedialBendingEnergyTerm(GenericMedialModel *model);
  
  /** Initialize the term with the template info */
  double ComputeEnergy(SolutionDataBase *data);

  /** Compute partial derivative */
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);
  
  void PrintReport(ostream &sout);
private:
  double xMaxBending, xTotalBending, xMeanBending;
};

class MedialRegularityTerm : public MedialIntegrationEnergyTerm
{
public:
  /** Initialize with a sample solution */
  MedialRegularityTerm(GenericMedialModel *model);
  
  /** Initialize the term with the template info */
  double ComputeEnergy(SolutionDataBase *data);

  // Compute partial derivative
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);
  
  void PrintReport(ostream &sout);
private:
  double xGradMagIntegral, xGradMagMaximum;
};

class MedialOptimizationProblem : public DifferentiableFunction
{
public:
  /** Initialize the problem */
  MedialOptimizationProblem(
    GenericMedialModel *xMedialModel, CoefficientMapping *xCoeff);

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
  typedef CoefficientMapping::Vec Vec;

  // Precision at which the problem is solved
  static const double xPrecision;

  // The epsilon used to compute finite difference derivatives
  static const double xEpsilon;

  // The number of coefficients
  size_t nCoeff;

  // Value of the last match
  double xLastSolutionValue;

  // Value of the solution for each of the energy terms
  Vec xLastTermValues;

  // The value of the coefficients before optimization began
  Vec xInitialCoefficients;

  // The medial model being optimized 
  GenericMedialModel *xMedialModel;

  // The mapping from free parameters to the model's coefficients (e.g.,
  // affine transform)
  CoefficientMapping *xCoeff;

  // The weights of the different energy terms in the optimization
  vector<double> xWeights;

  // The pointers to the different energy terms
  vector<EnergyTerm *> xTerms;

  // Timers used to report code speeds
  vector<CodeTimer> xTimers, xGradTimers;
  CodeTimer xSolveTimer, xSolveGradTimer, xWeightsTimer, xWeightsGradTimer;

  // Whether the gradient is available
  bool flagLastEvalAvailable, flagPhiGuessAvailable;
  
  // Last place where the function was evaluated
  vnl_vector<double> xLastEvalPoint;

  // Last place where the gradient was evaluated and its value there
  vnl_vector<double> xLastGradPoint, xLastGradient, xLastGradHint;
  bool flagGradientComputed;

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
