#include "OptimizationTerms.h"
#include "ITKImageWrapper.h"
#include "MedialAtomGrid.h"
#include <iostream>

using namespace std;

namespace medialpde {

/**
 * A simple image adapter that samples the image without applying
 * any additional processing
 */
class FloatImageEuclideanFunctionAdapter : public EuclideanFunction
{
public:
  FloatImageEuclideanFunctionAdapter(FloatImage *image)
    { this->image = image; this->xBackground = image->xOutsideValue; }

  double Evaluate(const SMLVec3d &X)
    { return image->xImage->Interpolate(X[0], X[1], X[2], xBackground); }

  SMLVec3d ComputeGradient(const SMLVec3d &X)
    {
    return SMLVec3d(
      image->xGradient[0]->Interpolate(X[0], X[1], X[2], 0.0f),
      image->xGradient[1]->Interpolate(X[0], X[1], X[2], 0.0f),
      image->xGradient[2]->Interpolate(X[0], X[1], X[2], 0.0f));
    }
protected:
  FloatImage *image;
  float xBackground;
};

/** An image adapter used for volume computations. For images in the -1.0 to
 * 1.0 range, this adapter shifts it into a 0.0 to 1.0 range */
class FloatImageVolumeElementFunctionAdapter :
  public FloatImageEuclideanFunctionAdapter
{
public:
  FloatImageVolumeElementFunctionAdapter(FloatImage *image) :
    FloatImageEuclideanFunctionAdapter(image) {}

  // Y = 0.5 (X + 1)
  double Evaluate(const SMLVec3d &X)
    { return 0.5 * (1.0 + FloatImageEuclideanFunctionAdapter::Evaluate(X)); }

  SMLVec3d ComputeGradient(const SMLVec3d &X)
    { return 0.5 * FloatImageEuclideanFunctionAdapter::ComputeGradient(X); }
};

/** A function that computes the match as (I(x) - I0)^2, where I0 is some
 * level set in the image. Used for blurred binary images */
class FloatImageLevelSetFunctionAdapter : 
  public FloatImageEuclideanFunctionAdapter
{
public:
  FloatImageLevelSetFunctionAdapter(FloatImage *image, float xLevelSet = 0.0f) : 
    FloatImageEuclideanFunctionAdapter(image) 
    { this->xLevelSet = xLevelSet; }

  double Evaluate(const SMLVec3d &X)
    { 
    // Get the image match from the superclass
    double I = FloatImageEuclideanFunctionAdapter::Evaluate(X);
    return (I - xLevelSet) * (I - xLevelSet);
    }

  SMLVec3d ComputeGradient(const SMLVec3d &X)
    {
    // Get the gradient and match from the superclass
    SMLVec3d G = FloatImageEuclideanFunctionAdapter::ComputeGradient(X);
    double I = FloatImageEuclideanFunctionAdapter::Evaluate(X);

    // Scale the gradient by the match
    return (2.0 * (I - xLevelSet)) * G;
    }
private:
  float xLevelSet;
};

} // namespace
using namespace medialpde;



/*********************************************************************************
 * SOLUTION DATA                  
 ********************************************************************************/
SolutionData::SolutionData(MedialPDESolver *xSolver, bool flagCopyAtoms)
{
  // Point to the atom grid
  xAtomGrid = xSolver->GetAtomGrid();
  
  // Back up the atoms or point ot them
  flagOwnAtoms = flagCopyAtoms;
  if(flagOwnAtoms)
    {
    xAtoms = new MedialAtom[xAtomGrid->GetNumberOfAtoms()];
    for(size_t i = 0; i < xAtomGrid->GetNumberOfAtoms(); i++)
      xAtoms[i] = xSolver->GetAtomArray()[i];
    }
  else
    xAtoms = xSolver->GetAtomArray();

  // Set the flags
  flagInternalWeights = false;
  flagBoundaryWeights = false;
    
  // Create and compute boundary weights
  UpdateBoundaryWeights();
}

void SolutionData::UpdateBoundaryWeights()
{
  // Only run this method once
  if(flagBoundaryWeights) return;
  
  // Create and compute boundary weights
  xBoundaryWeights = new double[xAtomGrid->GetNumberOfBoundaryPoints()];
  xBoundaryArea = 
    ComputeMedialBoundaryAreaWeights(xAtomGrid, xAtoms, xBoundaryWeights);

  // Set the flag
  flagBoundaryWeights = true;
}

void SolutionData::UpdateInternalWeights(size_t nCuts)
{
  // Only run this method once
  if(flagInternalWeights && nCuts == nInternalCuts) return;

  // Delete the old data if necessary
  if(flagInternalWeights) 
    { delete xInternalWeights; delete xInternalPoints; }
  
  // Allocate the internal points and weights
  nInternalPoints = xAtomGrid->GetNumberOfInternalPoints(nCuts);
  xInternalPoints = new SMLVec3d[nInternalPoints];
  xInternalWeights = new double[nInternalPoints];

  // Interpolate the internal points
  ComputeMedialInternalPoints(xAtomGrid, xAtoms, nCuts, xInternalPoints);
  xInternalVolume = ComputeMedialInternalVolumeWeights(
    xAtomGrid, xInternalPoints, nCuts, xInternalWeights);

  // Set the flag and store the nCuts
  flagInternalWeights = true;
  nInternalCuts = nCuts;
}

SolutionData::~SolutionData()
{
  if(flagBoundaryWeights)
    delete xBoundaryWeights;
  if(flagInternalWeights)
    { delete xInternalWeights; delete xInternalPoints; }
  if(flagOwnAtoms)
    delete xAtoms;
}

/*********************************************************************************
 * ENERGY TERM
 ********************************************************************************/
double EnergyTerm::ComputeGradient(
  SolutionData *S0, SolutionData **SGrad, 
  double xEpsilon, vnl_vector<double> &xOutGradient)
{
  // Compute the solution
  double xSolution = ComputeEnergy(S0);
  
  // Take the finite differences
  for(size_t i = 0; i < xOutGradient.size(); i++)
    {
    double xDelta = ComputeEnergy(SGrad[i]) - xSolution;
    xOutGradient[i] = xDelta / xEpsilon;
    }

  return xSolution;
}

/*********************************************************************************
 * BOUNDARY IMAGE MATCH TERM
 ********************************************************************************/
double 
BoundaryImageMatchTerm
::ComputeEnergy(SolutionData *S) 
{
  // Create the image adapter for image / gradient interpolation
  FloatImageLevelSetFunctionAdapter fImage(xImage);

  // Integrate the image match
  xImageMatch = IntegrateFunctionOverBoundary(
    S->xAtomGrid, S->xAtoms, 
    S->xBoundaryWeights, &fImage);

  xBoundaryArea = S->xBoundaryArea;
  xFinalMatch = xImageMatch / xBoundaryArea;

  // Scale by area
  return xFinalMatch;
}
  
// Print a verbose report
void 
BoundaryImageMatchTerm
::PrintReport(ostream &sout)
{
  sout << "  Boundary Image Match Term " << endl;
  sout << "    image match   : " << xImageMatch << endl;
  sout << "    boundary area : " << xBoundaryArea << endl;
  sout << "    ratio         : " << xFinalMatch << endl;
}

double
BoundaryImageMatchTerm::ComputeGradient(
  SolutionData *S0, SolutionData **SGrad, 
  double xEpsilon, vnl_vector<double> &xOutGradient)
{
  // Create the image adapter for image / gradient interpolation
  FloatImageLevelSetFunctionAdapter fImage(xImage);

  // Compute the image and image gradient at each point in the image
  SMLVec3d *xGradI = new SMLVec3d[S0->xAtomGrid->GetNumberOfBoundaryPoints()];
  double *xImageVal = new double[S0->xAtomGrid->GetNumberOfBoundaryPoints()];
  double xMatch = 0.0;
  MedialBoundaryPointIterator *itBnd = S0->xAtomGrid->NewBoundaryPointIterator();
  while(!itBnd->IsAtEnd())
    {
    // Compute the image gradient
    SMLVec3d &X = GetBoundaryPoint(itBnd, S0->xAtoms).X;
    xGradI[itBnd->GetIndex()] = fImage.ComputeGradient(X);

    // Compute the image value
    xImageVal[itBnd->GetIndex()] = fImage.Evaluate(X);

    // Accumulate to get weighted match
    xMatch += xImageVal[itBnd->GetIndex()] * S0->xBoundaryWeights[itBnd->GetIndex()];
    
    ++(*itBnd);
    }
  
  // We will need the area in many calculations
  double A0 = S0->xBoundaryArea;

  // Clear the gradient vector
  xOutGradient.fill(0.0);

  // Repeat for each coefficient
  for(size_t iCoeff = 0; iCoeff < xOutGradient.size(); iCoeff++)
    {
    // Accumulator for the partial derivative of the weighted match function
    double dMatchdC = 0.0;
    
    // Compute the partial derivative for this coefficient
    for(itBnd->GoToBegin(); !itBnd->IsAtEnd(); ++(*itBnd))
      {
      // Get the index of this boundary point
      size_t iPoint = itBnd->GetIndex();
      
      // Get the boundary point before and after small deformation
      SMLVec3d X0 = GetBoundaryPoint(itBnd, S0->xAtoms).X;
      SMLVec3d X1 = GetBoundaryPoint(itBnd, SGrad[iCoeff]->xAtoms).X;

      // Get the area weights for this point
      double w0 = S0->xBoundaryWeights[ iPoint ];
      double w1 = SGrad[iCoeff]->xBoundaryWeights[ iPoint ];

      // Compute the change in intensity per change in coefficient
      double dIdC = dot_product( xGradI[iPoint] , X1 - X0 );
      dMatchdC += dIdC * w0 + xImageVal[iPoint] * (w1 - w0);
      }

    // Scale by epsilon to get the derivative
    dMatchdC /= xEpsilon;

    // Compute the derivative of M / A
    double dAdC = (SGrad[iCoeff]->xBoundaryArea - A0) / xEpsilon;
    xOutGradient[iCoeff] = (dMatchdC * A0 - dAdC * xMatch) / (A0 * A0);
    }

  // Clean up
  delete xGradI;
  delete xImageVal;
  delete itBnd;

  // Return the solution
  return xMatch / A0;
}

/*********************************************************************************
 * BOUNDARY JACOBIAN TERM
 ********************************************************************************/
double BoundaryJacobianEnergyTerm::ComputeEnergy(SolutionData *S)
{
  // Place to store the Jacobian
  xTotalPenalty = 0.0;
  xMaxJacobian = 1.0, xMinJacobian = 1.0, xAvgJacobian = 0.0;
  
  // Create a Quad-iterator through the atoms
  MedialQuadIterator *itQuad = S->xAtomGrid->NewQuadIterator();
  while(!itQuad->IsAtEnd())
    {
    // Get the four atoms in this quad
    MedialAtom &A00 = S->xAtoms[itQuad->GetAtomIndex(0, 0)];
    MedialAtom &A01 = S->xAtoms[itQuad->GetAtomIndex(0, 1)];
    MedialAtom &A10 = S->xAtoms[itQuad->GetAtomIndex(1, 0)];
    MedialAtom &A11 = S->xAtoms[itQuad->GetAtomIndex(1, 1)];

    // Compute area element vectors on the medial surface
    SMLVec3d NM00 = cross_3d(A01.X - A00.X, A10.X - A00.X);
    SMLVec3d NM11 = cross_3d(A11.X - A10.X, A11.X - A01.X);
    
    // Compute the jacobians for each boundary side
    for(size_t k = 0; k < 2; k++)
      {
      // Compute area element vectors on the boundary side
      SMLVec3d NB00 = 
        cross_3d(A01.xBnd[k].X - A00.xBnd[k].X, A10.xBnd[k].X - A00.xBnd[k].X);
      SMLVec3d NB11 = 
        cross_3d(A11.xBnd[k].X - A10.xBnd[k].X, A11.xBnd[k].X - A01.xBnd[k].X);

      // Compute the 'Jacobians'
      double xSign = (k == 0) ? 1 : -1;
      double J00 = dot_product(NM00, NB00) / dot_product(NM00, NM00);
      double J11 = dot_product(NM11, NB11) / dot_product(NM11, NM11);

      // Store the smallest and largest Jacobian values
      if(J00 < xMinJacobian) xMinJacobian = J00;
      if(J11 < xMinJacobian) xMinJacobian = J11;
      if(J00 > xMaxJacobian) xMaxJacobian = J00;
      if(J11 > xMaxJacobian) xMaxJacobian = J11;

      // Add to the average Jacobian
      xAvgJacobian += J00 + J11;
        
      // Compute the penalty function
      xTotalPenalty += PenaltyFunction(J00, 10, 40);
      xTotalPenalty += PenaltyFunction(J11, 10, 40);
      }

    ++(*itQuad);
    }

  // Scale the average jacobian
  xAvgJacobian /= 2.0 * S->xAtomGrid->GetNumberOfBoundaryQuads();

  // Return the total value
  return xTotalPenalty;
}

// Print a verbose report
void 
BoundaryJacobianEnergyTerm
::PrintReport(ostream &sout)
{
  sout << "  Boundary Jacobian Term " << endl;
  sout << "    total match  : " << xTotalPenalty << endl;
  sout << "    min jacobian : " << xMinJacobian << endl;  
  sout << "    max jacobian : " << xMaxJacobian << endl;  
  sout << "    avg jacobian : " << xAvgJacobian << endl;  
}

/*********************************************************************************
 * VOLUME OVERLAP IMAGE MATCH TERM
 ********************************************************************************/
VolumeOverlapEnergyTerm::VolumeOverlapEnergyTerm( FloatImage *xImage, size_t nCuts)
{
  // Store the inputs
  this->nCuts = nCuts;
  this->xImage = xImage;

  // Compute the image volume (it remains a constant throughout)
  xImageVolume = xImage->ComputeObjectVolume();
}

double VolumeOverlapEnergyTerm::ComputeEnergy(SolutionData *data)
{
  // Make sure that the solution data has internal weights and points
  data->UpdateInternalWeights(nCuts);

  // Compute the match with the image function
  FloatImageVolumeElementFunctionAdapter fImage(xImage);
  
  // Compute the intersection between the image and m-rep
  xIntersect = IntegrateFunctionOverInterior(
    data->xAtomGrid, data->xInternalPoints, 
    data->xInternalWeights, nCuts, &fImage);

  // Get (remember) the object volume 
  xObjectVolume = data->xInternalVolume;

  // Compute the volume overlap measure
  xOverlap = xIntersect / (xImageVolume + xObjectVolume - xIntersect); 

  // Return a value that can be minimized
  return 1.0 - xOverlap;
}

double
VolumeOverlapEnergyTerm::ComputeGradient(
  SolutionData *S0, SolutionData **SGrad, 
  double xEpsilon, vnl_vector<double> &xOutGradient)
{
  // Get info about the problem
  MedialAtomGrid *xGrid = S0->xAtomGrid;
  size_t nPoints = xGrid->GetNumberOfBoundaryPoints();

  // Compute the raw image match at this point (we will use the stored terms)
  ComputeEnergy(S0);

  // Create the image adapter for image interpolation at the boundary
  FloatImageVolumeElementFunctionAdapter fImage(xImage);

  // Interpolate the image at each sample point on the surface
  double *xImageVal = new double[xGrid->GetNumberOfBoundaryPoints()];
  MedialBoundaryPointIterator *itBnd = xGrid->NewBoundaryPointIterator();
  for( ; !itBnd->IsAtEnd() ; ++(*itBnd) )
    {
    // Compute and store the image value at the sample point
    SMLVec3d &X = GetBoundaryPoint(itBnd, S0->xAtoms).X;
    xImageVal[itBnd->GetIndex()] = fImage.Evaluate(X);
    }
  
  // Clear the gradient vector
  xOutGradient.fill(0.0);

  // Repeat for each coefficient
  for(size_t iCoeff = 0; iCoeff < xOutGradient.size(); iCoeff++)
    {
    // We will compute the partial derivative of the absolute overlap and the
    // partial derivative of the m-rep volume with respect to the coeff
    double dOverlap = 0.0, dVolume = 0.0;
    
    // Compute the partial derivative for this coefficient
    for(itBnd->GoToBegin(); !itBnd->IsAtEnd(); ++(*itBnd))
      {
      // Get the index of this boundary point
      size_t iPoint = itBnd->GetIndex();
      
      // Get the boundary point before and after small deformation
      BoundaryAtom &B = GetBoundaryPoint(itBnd, S0->xAtoms);
      SMLVec3d X1 = GetBoundaryPoint(itBnd, SGrad[iCoeff]->xAtoms).X;

      // Get the area weights for this point
      double w = S0->xBoundaryWeights[ iPoint ];

      // Compute the common term between the two
      double z = dot_product(X1 - B.X, B.N) * w;

      // Compute the two partials
      dOverlap += z * xImageVal[iPoint];
      dVolume += z;
      }

    // Scale the partials by epsilon to get the derivative
    dOverlap /= xEpsilon; dVolume /= xEpsilon;

    // Compute the partial derivative of relative overlap match
    double xUnion = (xImageVolume + xObjectVolume - xIntersect);        
    xOutGradient[iCoeff] = 
      (dVolume * xIntersect - dOverlap * (xImageVolume + xObjectVolume));
    xOutGradient[iCoeff] /= xUnion * xUnion;
    }

  // Clean up
  delete xImageVal;
  delete itBnd;

  // Return the solution
  return 1.0 - xOverlap;
}

void VolumeOverlapEnergyTerm::PrintReport(ostream &sout)
{
  sout << "  Volume Overlap Term: " << endl;
  sout << "    image volume   : " << xImageVolume << endl;
  sout << "    object volume  : " << xObjectVolume << endl;
  sout << "    intersection   : " << xIntersect << endl;
  sout << "    overlap ratio  : " << xOverlap << endl;
}
  
/*********************************************************************************
 * Crest Laplacian Penalty Term: Assure negative Rho on the crest
 ********************************************************************************/
double CrestLaplacianEnergyTerm::ComputeEnergy(SolutionData *data)
{
  // Initialize the accumulators
  xMaxLaplacian = -1e10;
  xTotalPenalty = xAvgLaplacian = 0.0;
  nCrestAtoms = nBadSites = 0;
  
  // Iterate over all crest atoms
  MedialAtomIterator *itAtom = data->xAtomGrid->NewAtomIterator();
  for( ; !itAtom->IsAtEnd(); ++(*itAtom) )
    {
    if(itAtom->IsEdgeAtom())
      {
      // Get the laplacian at this point
      double x = data->xAtoms[itAtom->GetIndex()].xLapR;

      // Update the minimum value
      if(x > xMaxLaplacian) xMaxLaplacian = x;

      // Update the bad atom counter 
      if(x > 0) nBadSites++;
      nCrestAtoms++;

      // Update the average
      xAvgLaplacian += x;
      
      // Update the total penalty
      xTotalPenalty += PenaltyFunction(x, 100, 0.0);
      }
    }
  delete itAtom;

  // Finish up
  xAvgLaplacian /= nCrestAtoms;
  xTotalPenalty;

  return xTotalPenalty;
}  

void CrestLaplacianEnergyTerm::PrintReport(ostream &sout)
{
  sout << "  Crest Laplacian Penalty Term: " << endl;
  sout << "    number of crest atoms    : " << nCrestAtoms << endl; 
  sout << "    number of bad atoms      : " << nBadSites << endl; 
  sout << "    largest laplacian value  : " << xMaxLaplacian << endl; 
  sout << "    average laplacian value  : " << xAvgLaplacian << endl; 
  sout << "    total penalty            : " << xTotalPenalty << endl;
}


/*********************************************************************************
 * MEDIAL OPTIMIZATION PROBLEM
 ********************************************************************************/
const double MedialOptimizationProblem::xPrecision = 1.0e-12;
const double MedialOptimizationProblem::xEpsilon = 1.0e-6;

// Add the energy term
void MedialOptimizationProblem::AddEnergyTerm(EnergyTerm *term, double xWeight)
{
  xWeights.push_back(xWeight);
  xTerms.push_back(term);
}


double MedialOptimizationProblem::Evaluate(double *X)
{
  // Update the surface with the new solution
  xCoeff->SetCoefficientArray(X);
  xSolver->Solve(xPrecision);

  // Create a solution data object representing current solution
  SolutionData S0(xSolver, false);

  // cout << "REPORT:" << endl;

  // Compute the solution for each term
  xSolution = 0.0;
  for(size_t iTerm = 0; iTerm < xTerms.size(); iTerm++)
    { xSolution += xTerms[iTerm]->ComputeEnergy(&S0) * xWeights[iTerm]; }

  // cout << "RESULT:" << xValue << endl;

  // Return the result
  return xSolution;
}

double 
MedialOptimizationProblem
::ComputeGradient(double *X, double *XGradient)
{
  // Copy the x-vector into the coefficients of the surface
  size_t nCoeff = xCoeff->GetNumberOfCoefficients();
  xCoeff->SetCoefficientArray(X);

  // Use current solution as the guess
  xSolver->Solve(xPrecision);
  xSolver->SetSolutionAsInitialGuess();

  // Create a solution data object representing current solution
  SolutionData *S0 = new SolutionData(xSolver, true);

  // Create an array of solution data objects
  SolutionData **SGrad = new SolutionData *[nCoeff];

  // Compute the time it takes to compute the derivative m-reps
  double tStart = clock();

  // Repeat for each coefficient
  for(size_t iCoeff = 0; iCoeff < nCoeff; iCoeff++)
    {
    // Increment the coefficient by epsilon
    double xOldValue = xCoeff->GetCoefficient(iCoeff);
    xCoeff->SetCoefficient(iCoeff, xOldValue + xEpsilon);

    // Solve the PDE for the new value and create a solution data object
    xSolver->Solve(xPrecision); cout << "." << flush;
    SGrad[iCoeff] = new SolutionData(xSolver, true);

    // Restore the old solution value
    xCoeff->SetCoefficient(iCoeff, xOldValue);
    }

  // See how much time elapsed for the medial computation
  double tMedial = clock() - tStart; 

  // At this point, we have computed the solution at current X and at a set of 
  // X + eps positions. We can use this to compute the gradient of each of the
  // terms involved in the optimization
  xSolution = 0.0; 
  vnl_vector<double> xGradTotal(nCoeff, 0.0), xGradTerm(nCoeff, 0.0);

  // Add up all the terms and all the gradients
  cout << endl;
  for(size_t iTerm = 0; iTerm < xTerms.size(); iTerm++)
    {
    double xValue = xTerms[iTerm]->ComputeGradient(S0, SGrad, xEpsilon, xGradTerm);
    // cout << xGradTerm << endl;
    xSolution += xWeights[iTerm] * xValue;
    xGradTotal += xWeights[iTerm] * xGradTerm;

    
    cout << "GRADIENT : " << xGradTerm << endl;
    }

  // Finish timing
  double tGradient = clock() - tMedial - tStart;
  cout << endl << "PDE Time: " << (tMedial / CLOCKS_PER_SEC) <<
   "; Gradient Time: " << (tGradient / CLOCKS_PER_SEC) << endl;

  // Increment the evaluation counter
  evaluationCost += nCoeff;

  // Clean up everything
  for(size_t iCoeff = 0; iCoeff < nCoeff; iCoeff++)
    delete SGrad[iCoeff];
  delete S0; delete SGrad;

  // Return the solution value
  for(size_t j = 0; j < nCoeff; j++) XGradient[j] = xGradTotal[j];
  return xSolution;
}

void MedialOptimizationProblem::PrintReport(ostream &sout)
{
  sout << "Optimization Summary: " << endl;
  sout << "  energy value   :  " << xSolution << endl; 
  sout << "  evaluation cost: " << this->getEvaluationCost() << endl;
  sout << "Per-Term Report:" << endl;

  for(size_t iTerm = 0; iTerm < xTerms.size(); iTerm++)
    xTerms[iTerm]->PrintReport(sout);
}

