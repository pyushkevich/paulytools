#include "OptimizationTerms.h"
#include "ITKImageWrapper.h"
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



/**
 * This function computes the area associated with each triangle on the
 * boundary of the object and assigns a third of this area to each of the 
 * vertices of the triangle. Thus for each vertex on the boundary a weight
 * is generated. The total of these weights, equal to the area of the boundary
 * surface is the return value 
 */
double ComputeMedialBoundaryAreaWeights( MedialAtomGrid *xGrid, 
  MedialAtom *xAtoms, double *xWeights)
{
  // A constant to hold 1/3
  const static double THIRD = 1.0f / 3.0f;
  double xTotalArea = 0.0f;

  // Clear the content of the weights
  memset(xWeights, 0, sizeof(double) * xGrid->GetNumberOfBoundaryPoints());
  
  // Create a quad-based iterator
  MedialBoundaryQuadIterator *itQuad = xGrid->NewBoundaryQuadIterator();

  // For each quad, compute the area associated with it
  while(!itQuad->IsAtEnd())
    {
    // Access the four medial points
    SMLVec3d X00 = GetBoundaryPoint(itQuad, xAtoms, 0, 0).X;
    SMLVec3d X01 = GetBoundaryPoint(itQuad, xAtoms, 0, 1).X;
    SMLVec3d X11 = GetBoundaryPoint(itQuad, xAtoms, 1, 1).X;
    SMLVec3d X10 = GetBoundaryPoint(itQuad, xAtoms, 1, 0).X;

    // Compute the areas of the two triangles involved
    double A1 = THIRD * TriangleArea( X00, X01, X11 );
    double A2 = THIRD * TriangleArea( X00, X11, X10 );
    xTotalArea += (A1 + A2);

    // Assign a third of each weight to each corner
    xWeights[itQuad->GetBoundaryIndex(0, 0)] += A1 + A2;
    xWeights[itQuad->GetBoundaryIndex(1, 1)] += A1 + A2;
    xWeights[itQuad->GetBoundaryIndex(0, 1)] += A1;
    xWeights[itQuad->GetBoundaryIndex(1, 0)] += A2;

    // Go to the next boundary quad
    ++(*itQuad);
    }

  // Delete the iterator object
  delete itQuad;

  // Return the total area
  return 3.0f * xTotalArea;
}

double IntegrateFunctionOverBoundary(
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, double *xWeights, 
  EuclideanFunction *fMatch) 
{
  // Match accumulator
  double xMatch = 0.0;
  
  // Create a boundary point iterator
  MedialBoundaryPointIterator *itBoundary = xGrid->NewBoundaryPointIterator();

  // Iterate through all boundary points
  while(!itBoundary->IsAtEnd())
    {
    // Evaluate the image at this location
    double xLocal = fMatch->Evaluate(GetBoundaryPoint(itBoundary, xAtoms).X);
    xMatch += xLocal * xWeights[itBoundary->GetIndex()];

    // On to the next point
    ++(*itBoundary); 
    }

  // Clean up
  delete itBoundary;

  // Return match scaled by total weight
  return xMatch;
}

/** Compute gradient of a function over the boundary */
double ComputeFunctionJetOverBoundary(
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, 
  double *xWeights, EuclideanFunction *fMatch, SMLVec3d *xOutGradient) 
{
  // Compute the image match
  double xMatch = IntegrateFunctionOverBoundary(xGrid, xAtoms, xWeights, fMatch);
  
  // Create a boundary point iterator and iterate over all points
  MedialBoundaryPointIterator *itBoundary = xGrid->NewBoundaryPointIterator();

  // Iterate through all boundary points
  while(!itBoundary->IsAtEnd())
    {
    // Evaluate the image gradient at this point
    unsigned int iBnd = itBoundary->GetIndex();
    SMLVec3d X = GetBoundaryPoint(itBoundary, xAtoms).X;
    xOutGradient[iBnd] = xWeights[iBnd] * fMatch->ComputeGradient(X);

    // On to the next point
    ++(*itBoundary);
    }

  // Return match scaled by total weight
  return xMatch;
}

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
    
  // Create and compute boundary weights
  xBoundaryWeights= new double[xAtomGrid->GetNumberOfBoundaryPoints()];
  xBoundaryArea = 
    ComputeMedialBoundaryAreaWeights(xAtomGrid, xAtoms, xBoundaryWeights);
}

SolutionData::~SolutionData()
{
  delete xBoundaryWeights;
  if(flagOwnAtoms)
    delete xAtoms;
}


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


double 
BoundaryImageMatchTerm
::ComputeEnergy(SolutionData *S) 
{
  // Create the image adapter for image / gradient interpolation
  FloatImageLevelSetFunctionAdapter fImage(xImage);

  // Integrate the image match
  double xMatch = IntegrateFunctionOverBoundary(
    S->xAtomGrid, S->xAtoms, 
    S->xBoundaryWeights, &fImage);

  // Scale by area
  return xMatch / S->xBoundaryArea;
}
  
// Print a verbose report
void 
BoundaryImageMatchTerm
::PrintReport(ostream &sout, SolutionData *data)
{
  double x = ComputeEnergy(data);
  sout << "  Boundary Image Match Term " << endl;
  sout << "    total match  : " << x << endl;
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

double BoundaryJacobianEnergyTerm::ComputeEnergy(SolutionData *S)
{
  // Place to store the Jacobian
  double xTotalPenalty = 0.0;
  xMaxJacobian = 1.0, xMinJacobian = 1.0;
  
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
        
      // Compute the penalty function
      xTotalPenalty += PenaltyFunction(J00, 10, 40);
      xTotalPenalty += PenaltyFunction(J11, 10, 40);
      }

    ++(*itQuad);
    }

  // Return the total value
  return xTotalPenalty;
}

// Print a verbose report
void 
BoundaryJacobianEnergyTerm
::PrintReport(ostream &sout, SolutionData *data)
{
  double x = ComputeEnergy(data);
  sout << "  Boundary Jacobian Term " << endl;
  sout << "    total match  : " << x << endl;
  sout << "    min jacobian : " << xMinJacobian << endl;  
  sout << "    max jacobian : " << xMaxJacobian << endl;  
}

// Constants
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
  xSurface->SetRawCoefficientArray(X);
  xSolver->Solve(xPrecision);

  // Create a solution data object representing current solution
  SolutionData S0(xSolver, false);

  // cout << "REPORT:" << endl;

  // Compute the solution for each term
  double xValue = 0.0;
  for(size_t iTerm = 0; iTerm < xTerms.size(); iTerm++)
    { 
    xValue += xTerms[iTerm]->ComputeEnergy(&S0) * xWeights[iTerm]; 
    // xTerms[iTerm]->PrintReport(cout, &S0);
    }

  // cout << "RESULT:" << xValue << endl;

  // Return the result
  return xValue;
}

double 
MedialOptimizationProblem
::ComputeGradient(double *X, double *XGradient)
{
  // Copy the x-vector into the coefficients of the surface
  size_t nCoeff = xSurface->GetNumberOfRawCoefficients();
  xSurface->SetRawCoefficientArray(X);

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
    double xOldValue = xSurface->GetRawCoefficient(iCoeff);
    xSurface->SetRawCoefficient(iCoeff, xOldValue + xEpsilon);

    // Solve the PDE for the new value and create a solution data object
    xSolver->Solve(xPrecision); cout << "." << flush;
    SGrad[iCoeff] = new SolutionData(xSolver, true);

    // Restore the old solution value
    xSurface->SetRawCoefficient(iCoeff, xOldValue);
    }

  // See how much time elapsed for the medial computation
  double tMedial = clock() - tStart; 

  // At this point, we have computed the solution at current X and at a set of 
  // X + eps positions. We can use this to compute the gradient of each of the
  // terms involved in the optimization
  double xSolution = 0.0; 
  vnl_vector<double> xGradTotal(nCoeff, 0.0), xGradTerm(nCoeff, 0.0);

  // Add up all the terms and all the gradients
  cout << endl << "REPORT" << endl;
  for(size_t iTerm = 0; iTerm < xTerms.size(); iTerm++)
    {
    double xValue = xTerms[iTerm]->ComputeGradient(S0, SGrad, xEpsilon, xGradTerm);
    cout << xGradTerm << endl;
    xSolution += xWeights[iTerm] * xValue;
    xGradTotal += xWeights[iTerm] * xGradTerm;
    
    xTerms[iTerm]->PrintReport(cout, S0);
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

