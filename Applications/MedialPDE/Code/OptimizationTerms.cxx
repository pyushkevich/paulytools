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
    { this->image = image; }

  double Evaluate(const SMLVec3d &X)
    { return (double) image->Interpolate(X); }
  
  void ComputeGradient(const SMLVec3d &X, SMLVec3d &G)
    { image->InterpolateImageGradient(X, G); }
  
protected:
  FloatImage *image;
};

/** A function that computes the match as (I(x) - I0)^2, where I0 is some
 * level set in the image. Used for blurred binary images */
class FloatImageSquareValueFunctionAdapter : public EuclideanFunction
{
public:
  FloatImageSquareValueFunctionAdapter(FloatImage *image, float xLevelSet = 0.0f) 
    { this->image = image; this->xLevelSet = xLevelSet; }

  double Evaluate(const SMLVec3d &X)
    { 
    // Get the image match from the superclass
    double I = (double) image->Interpolate(X);
    double d = I - xLevelSet; 
    return d * d;
    }

  void ComputeGradient(const SMLVec3d &X, SMLVec3d &G)
    {
    // Get the gradient and match from the superclass
    double I = (double) image->Interpolate(X);
    double d = I - xLevelSet;
    image->InterpolateImageGradient(X, G);
    G *= 2.0 * d;
    }

private:
  FloatImage *image;
  float xLevelSet;
};

} // namespace
using namespace medialpde;

/*********************************************************************************
 * SOLUTION DATA BASE CLASS
 ********************************************************************************/
SolutionDataBase::SolutionDataBase()
{
  // Set the flags
  flagInternalWeights = false;
  flagBoundaryWeights = false;
  flagOwnAtoms = false;
  
  xAtomGrid = NULL;
  xAtoms = NULL;
  xInternalPoints = NULL;
  xBoundaryWeights = NULL;
  xBoundaryArea = 0.0;
  xInternalWeights = NULL;
  xInternalVolume = 0.0;
  xInternalProfileWeights = NULL;  
}

SolutionDataBase::~SolutionDataBase()
{
  if(flagBoundaryWeights)
    { 
    delete xBoundaryWeights; 
    }
  if(flagInternalWeights)
    { 
    delete xInternalWeights; 
    delete xInternalPoints; 
    delete xInternalProfileWeights; 
    }
  if(flagOwnAtoms)
    {
    delete xAtoms;
    }
}
/*********************************************************************************
 * SOLUTION DATA                  
 ********************************************************************************/
SolutionData::SolutionData(MedialPDESolver *xSolver)
: SolutionDataBase()
{
  // Point to the atom grid
  xAtomGrid = xSolver->GetAtomGrid();
  
  // Back up the atoms or point ot them
  xAtoms = xSolver->GetAtomArray();

  // Set the flags
  flagInternalWeights = false;
  flagBoundaryWeights = false;
  flagOwnAtoms = false;
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
    { delete xInternalWeights; delete xInternalPoints; delete xInternalProfileWeights; }
  
  // Allocate the internal points and weights
  nInternalPoints = xAtomGrid->GetNumberOfInternalPoints(nCuts);
  nProfileIntervals = xAtomGrid->GetNumberOfProfileIntervals(nCuts);
  
  xInternalPoints = new SMLVec3d[nInternalPoints];
  xInternalWeights = new double[nInternalPoints];
  xInternalProfileWeights = new double[nProfileIntervals];

  // Interpolate the internal points - same method works for derivatives and
  // normal points because the interpolation is linear
  ComputeMedialInternalPoints(xAtomGrid, xAtoms, nCuts, xInternalPoints);

  // Compute the internal volumes - this requires the reference solution data
  // as well
  xInternalVolume = ComputeMedialInternalVolumeWeights(
    xAtomGrid, xInternalPoints, nCuts, xInternalWeights, xInternalProfileWeights);

  // Set the flag and store the nCuts
  flagInternalWeights = true;
  nInternalCuts = nCuts;
}

/*********************************************************************************
 * PARTIAL DERIVATIVE SOLUTION DATA
 ********************************************************************************/
PartialDerivativeSolutionData
::PartialDerivativeSolutionData(SolutionData *xReference, MedialAtom *dAtoms)
: SolutionDataBase()
{
  this->xReference = xReference;
  xAtomGrid = xReference->xAtomGrid;
  xAtoms = dAtoms;
}

void
PartialDerivativeSolutionData::UpdateBoundaryWeights()
{
  // Only run this method once
  if(flagBoundaryWeights) return;

  // Ensure that the reference weights exist
  xReference->UpdateBoundaryWeights();

  // Create boundary weights
  xBoundaryWeights = new double[xAtomGrid->GetNumberOfBoundaryPoints()];

  // We need parent information to run this method
  xBoundaryArea = ComputeMedialBoundaryAreaPartialDerivative(
    xAtomGrid, xReference->xAtoms, xAtoms, 
    xReference->xBoundaryWeights, xBoundaryWeights);
  
  // Set the flag
  flagBoundaryWeights = true;
}

void 
PartialDerivativeSolutionData::UpdateInternalWeights(size_t nCuts)
{
  // Only run this method once
  if(flagInternalWeights && nCuts == nInternalCuts) return;

  // Delete the old data if necessary
  if(flagInternalWeights) 
    { 
    delete xInternalWeights; 
    delete xInternalPoints; 
    delete xInternalProfileWeights; 
    }

  // Make sure the reference has the internal data available
  xReference->UpdateInternalWeights(nCuts);
  
  // Allocate the internal points and weights
  nInternalPoints = xAtomGrid->GetNumberOfInternalPoints(nCuts);
  nProfileIntervals = xAtomGrid->GetNumberOfProfileIntervals(nCuts);
  
  xInternalPoints = new SMLVec3d[nInternalPoints];
  xInternalWeights = new double[nInternalPoints];
  xInternalProfileWeights = new double[nProfileIntervals];

  // Interpolate the internal points
  ComputeMedialInternalPoints(xAtomGrid, xAtoms, nCuts, xInternalPoints);
  xInternalVolume = 
    ComputeMedialInternalVolumePartialDerivativeWeights(
      xAtomGrid, xReference->xInternalPoints, xInternalPoints, 
      nCuts, xInternalWeights, xInternalProfileWeights);

  // Set the flag and store the nCuts
  flagInternalWeights = true;
  nInternalCuts = nCuts;
}

/*********************************************************************************
 * OFFSET SOLUTION DATA
 ********************************************************************************/
 
OffsetSolutionData::OffsetSolutionData(
  SolutionData *sRef, 
  PartialDerivativeSolutionData *sDeriv,
  double eps)
{
  // Create pointers to reference objects
  xReference = sRef;
  xDerivative = sDeriv;
  xEpsilon = eps;
  
  // Point to the atom grid
  xAtomGrid = xReference->xAtomGrid;
  
  // Create a copy of the medial atoms
  xAtoms = new MedialAtom[xAtomGrid->GetNumberOfAtoms()];
  for(size_t iAtom = 0; iAtom < xAtomGrid->GetNumberOfAtoms(); iAtom++)
    AddScaleMedialAtoms(xReference->xAtoms[iAtom], 
      xDerivative->xAtoms[iAtom], xEpsilon, xAtoms[iAtom]);

  // Set the ownership flag
  flagOwnAtoms = true;    
}
  
void
OffsetSolutionData::UpdateBoundaryWeights()
{
  // Only run this method once
  if(flagBoundaryWeights) return;

  // Ensure that the reference weights exist
  xReference->UpdateBoundaryWeights();
  xDerivative->UpdateBoundaryWeights();

  // Create boundary weights
  xBoundaryWeights = new double[xAtomGrid->GetNumberOfBoundaryPoints()];

  // Compute the boundary weights by addition
  for(size_t i = 0; i < xAtomGrid->GetNumberOfBoundaryPoints(); i++)
    xBoundaryWeights[i] = xReference->xBoundaryWeights[i]
      + xEpsilon * xDerivative->xBoundaryWeights[i];
  
  // Compute the boundary area same way
  xBoundaryArea = xReference->xBoundaryArea
      + xEpsilon * xDerivative->xBoundaryArea;
  
  // Set the flag
  flagBoundaryWeights = true;
}

void 
OffsetSolutionData::UpdateInternalWeights(size_t nCuts)
{
  // Only run this method once
  if(flagInternalWeights && nCuts == nInternalCuts) return;

  // Delete the old data if necessary
  if(flagInternalWeights) 
    { 
    delete xInternalWeights; 
    delete xInternalPoints; 
    delete xInternalProfileWeights; 
    }

  // Make sure the reference has the internal data available
  xReference->UpdateInternalWeights(nCuts);
  xDerivative->UpdateInternalWeights(nCuts);
  
  // Allocate the internal points and weights
  nInternalPoints = xAtomGrid->GetNumberOfInternalPoints(nCuts);
  nProfileIntervals = xAtomGrid->GetNumberOfProfileIntervals(nCuts);
  
  xInternalPoints = new SMLVec3d[nInternalPoints];
  xInternalWeights = new double[nInternalPoints];
  xInternalProfileWeights = new double[nProfileIntervals];

  // Compute the internal weights by addition
  for(size_t i = 0; i < nInternalPoints; i++)
    xInternalWeights[i] = xReference->xInternalWeights[i]
      + xEpsilon * xDerivative->xInternalWeights[i];

  // Compute the profile interval weights
  for(size_t i = 0; i < nProfileIntervals; i++)
    xInternalProfileWeights[i] = xReference->xInternalProfileWeights[i]
      + xEpsilon * xDerivative->xInternalProfileWeights[i];

  // Compute the total weights
  xInternalVolume = xReference->xInternalVolume
      + xEpsilon * xDerivative->xInternalVolume;

  // Set the flag and store the nCuts
  flagInternalWeights = true;
  nInternalCuts = nCuts;
}


/*********************************************************************************
 * ENERGY TERM
 ********************************************************************************/

double 
NumericalGradientEnergyTerm
::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  cout << "GENERATING OFFSET DATA!!!" << endl;
  OffsetSolutionData s1(S, dS, xEpsilon), s2(S, dS, -xEpsilon);
  
  // Compute the two energy values
  double f1 = ComputeEnergy(&s1);
  double f2 = ComputeEnergy(&s2);
  
  // Return the central difference derivative
  return 0.5 * (f1 - f2) / xEpsilon;
}

/*********************************************************************************
 * BOUNDARY IMAGE MATCH TERM
 ********************************************************************************/
double 
BoundaryImageMatchTerm
::ComputeEnergy(SolutionDataBase *S) 
{
  // Create the image adapter for image / gradient interpolation
  FloatImageSquareValueFunctionAdapter fImage(xImage);

  // Make sure the weights exist in the solution
  S->UpdateBoundaryWeights();

  // Integrate the image match
  xImageMatch = IntegrateFunctionOverBoundary(
    S->xAtomGrid, S->xAtoms, S->xBoundaryWeights, &fImage);

  // Compute the adaptive match, making sure the interpolation works
  // double xMinArea = 0.5 * S->xBoundaryArea / S->xAtomGrid->GetNumberOfBoundaryQuads();
  // xImageMatch = AdaptivelyIntegrateFunctionOverBoundary(
  //   S->xAtomGrid, S->xAtoms, xMinArea, &fImage);

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
BoundaryImageMatchTerm::BeginGradientComputation(SolutionData *S)
{
  // Create the image adapter for image / gradient interpolation
  FloatImageSquareValueFunctionAdapter fImage(xImage);

  // Boundary weights will be needed
  S->UpdateBoundaryWeights();

  // Compute the image and image gradient at each point in the image
  xGradI = new SMLVec3d[S->xAtomGrid->GetNumberOfBoundaryPoints()];
  xImageVal = new double[S->xAtomGrid->GetNumberOfBoundaryPoints()];
  xImageMatch = 0.0;

  MedialBoundaryPointIterator *itBnd = S->xAtomGrid->NewBoundaryPointIterator();
  while(!itBnd->IsAtEnd())
    {
    // Compute the image gradient
    SMLVec3d &X = GetBoundaryPoint(itBnd, S->xAtoms).X;
    fImage.ComputeGradient(X, xGradI[itBnd->GetIndex()]);

    // Compute the image value
    xImageVal[itBnd->GetIndex()] = fImage.Evaluate(X);

    // Accumulate to get weighted match
    xImageMatch += 
      xImageVal[itBnd->GetIndex()] * S->xBoundaryWeights[itBnd->GetIndex()];
    
    ++(*itBnd);
    }
  
  // We will need the area in many calculations
  xBoundaryArea = S->xBoundaryArea;

  // Compute the final match
  xFinalMatch = xImageMatch / xBoundaryArea;
  
  // Clean up
  delete itBnd;

  // Return the solution
  return xFinalMatch;
}

double
BoundaryImageMatchTerm
::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *DS)
{
  // Create the image adapter for image / gradient interpolation
  FloatImageSquareValueFunctionAdapter fImage(xImage);

  // Accumulator for the partial derivative of the weighted match function
  double dMatchdC = 0.0;

  // We need the weights for the derivative terms
  S->UpdateBoundaryWeights();
  DS->UpdateBoundaryWeights();

  // Compute the partial derivative for this coefficient
  MedialBoundaryPointIterator *itBnd = S->xAtomGrid->NewBoundaryPointIterator();
  for( ; !itBnd->IsAtEnd(); ++(*itBnd))
    {
    // Get the index of this boundary point
    size_t iPoint = itBnd->GetIndex();

    // Get the change in the boundary point
    SMLVec3d dX = GetBoundaryPoint(itBnd, DS->xAtoms).X;

    // Get the area weights for this point
    double w = S->xBoundaryWeights[ iPoint ];
    double dw = DS->xBoundaryWeights[ iPoint ];

    // Compute the change in intensity per change in coefficient
    double I = xImageVal[iPoint];
    double dIdC = dot_product( xGradI[iPoint] , dX);
    
    // Increment the partial derivative of the weighted match
    dMatchdC += dw * I + w * dIdC;
    }
  
  delete itBnd;

  // Compute the derivative
  double dFinaldC = 
    ( dMatchdC * S->xBoundaryArea - xImageMatch * DS->xBoundaryArea ) / 
    ( S->xBoundaryArea * S->xBoundaryArea );

  // Return the final match derivative
  return dFinaldC;
}

void BoundaryImageMatchTerm::EndGradientComputation()
{
  // Clean up
  delete xGradI;
  delete xImageVal;
}

/*********************************************************************************
 * BOUNDARY JACOBIAN TERM
 ********************************************************************************/

inline double ComputeJacobian(const SMLVec3d &Xu, const SMLVec3d &Xv, 
  const SMLVec3d &Yu, const SMLVec3d &Yv)
{
  return 
    ( dot_product(Yu,Xu) * dot_product(Yv,Xv) - 
      dot_product(Yu,Xv) * dot_product(Yv,Xu) ) / 
    ( dot_product(Xu,Xu) * dot_product(Xv,Xv) - 
      dot_product(Xu,Xv) * dot_product(Xu,Xv) );
}

double BoundaryJacobianEnergyTerm::ComputeEnergy(SolutionDataBase *S)
{
  // Place to store the Jacobian
  xTotalPenalty = 0.0;
  xMaxJacobian = 1.0, xMinJacobian = 1.0, xAvgJacobian = 0.0;
  
  // Create a Quad-iterator through the atoms
  MedialQuadIterator *itQuad = S->xAtomGrid->NewQuadIterator();
  for(; !itQuad->IsAtEnd(); ++(*itQuad))
    {
    // Get the four atoms in this quad
    MedialAtom &A00 = S->xAtoms[itQuad->GetAtomIndex(0, 0)];
    MedialAtom &A01 = S->xAtoms[itQuad->GetAtomIndex(0, 1)];
    MedialAtom &A10 = S->xAtoms[itQuad->GetAtomIndex(1, 0)];
    MedialAtom &A11 = S->xAtoms[itQuad->GetAtomIndex(1, 1)];

    // Compute the average Xu and Xv vectors
    SMLVec3d XU = 0.5 * ((A11.X - A01.X) + (A10.X - A00.X));
    SMLVec3d XV = 0.5 * ((A11.X - A10.X) + (A01.X - A00.X));

    // Compute the same for the upper and lower boundaries
    SMLVec3d Y0U = 
      0.5 * ((A11.xBnd[0].X - A01.xBnd[0].X) + (A10.xBnd[0].X - A00.xBnd[0].X));
    SMLVec3d Y0V = 
      0.5 * ((A11.xBnd[0].X - A10.xBnd[0].X) + (A01.xBnd[0].X - A00.xBnd[0].X));
    SMLVec3d Y1U = 
      0.5 * ((A11.xBnd[1].X - A01.xBnd[1].X) + (A10.xBnd[1].X - A00.xBnd[1].X));
    SMLVec3d Y1V = 
      0.5 * ((A11.xBnd[1].X - A10.xBnd[1].X) + (A01.xBnd[1].X - A00.xBnd[1].X));

    // Compute the scaled normal vectors
    SMLVec3d NX = vnl_cross_3d(XU, XV);
    SMLVec3d NY0 = vnl_cross_3d(Y0U, Y0V);
    SMLVec3d NY1 = vnl_cross_3d(Y1U, Y1V);

    // Compute the Jacobians
    double J0 = dot_product(NY0, NX) / dot_product(NX, NX);
    double J1 = dot_product(NY1, NX) / dot_product(NX, NX);

    // Store the smallest and largest Jacobian values
    if(J0 < xMinJacobian) xMinJacobian = J0;
    if(J1 < xMinJacobian) xMinJacobian = J1;
    if(J0 > xMaxJacobian) xMaxJacobian = J0;
    if(J1 > xMaxJacobian) xMaxJacobian = J1;

    // Add to the average Jacobian
    xAvgJacobian += J0 + J1;

    // Compute the penalty function
    xTotalPenalty += PenaltyFunction(J0, 10, 400);
    xTotalPenalty += PenaltyFunction(J1, 10, 400);
    }
  delete itQuad;

  // Scale the average jacobian
  xAvgJacobian /= S->xAtomGrid->GetNumberOfBoundaryQuads();

  // Return the total value
  return xTotalPenalty;
}

double BoundaryJacobianEnergyTerm
::ComputePartialDerivative(
  SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // Place to store the Jacobian
  double dTotalPenalty = 0.0;
  
  // Create a Quad-iterator through the atoms
  MedialQuadIterator *itQuad = S->xAtomGrid->NewQuadIterator();
  for(; !itQuad->IsAtEnd(); ++(*itQuad))
    {
    // Get the four atoms in this quad
    MedialAtom &A00 = S->xAtoms[itQuad->GetAtomIndex(0, 0)];
    MedialAtom &A01 = S->xAtoms[itQuad->GetAtomIndex(0, 1)];
    MedialAtom &A10 = S->xAtoms[itQuad->GetAtomIndex(1, 0)];
    MedialAtom &A11 = S->xAtoms[itQuad->GetAtomIndex(1, 1)];

    // Get the derivative atoms too
    MedialAtom &dA00 = dS->xAtoms[itQuad->GetAtomIndex(0, 0)];
    MedialAtom &dA01 = dS->xAtoms[itQuad->GetAtomIndex(0, 1)];
    MedialAtom &dA10 = dS->xAtoms[itQuad->GetAtomIndex(1, 0)];
    MedialAtom &dA11 = dS->xAtoms[itQuad->GetAtomIndex(1, 1)];
    
    // Compute the average Xu and Xv vectors and derivatives
    SMLVec3d XU = 0.5 * ((A11.X - A01.X) + (A10.X - A00.X));
    SMLVec3d XV = 0.5 * ((A11.X - A10.X) + (A01.X - A00.X));
    SMLVec3d dXU = 0.5 * ((dA11.X - dA01.X) + (dA10.X - dA00.X));
    SMLVec3d dXV = 0.5 * ((dA11.X - dA10.X) + (dA01.X - dA00.X));

    // Compute the same for the upper and lower boundaries
    SMLVec3d Y0U = 
      0.5 * ((A11.xBnd[0].X - A01.xBnd[0].X) + (A10.xBnd[0].X - A00.xBnd[0].X));
    SMLVec3d Y0V = 
      0.5 * ((A11.xBnd[0].X - A10.xBnd[0].X) + (A01.xBnd[0].X - A00.xBnd[0].X));
    SMLVec3d Y1U = 
      0.5 * ((A11.xBnd[1].X - A01.xBnd[1].X) + (A10.xBnd[1].X - A00.xBnd[1].X));
    SMLVec3d Y1V = 
      0.5 * ((A11.xBnd[1].X - A10.xBnd[1].X) + (A01.xBnd[1].X - A00.xBnd[1].X));

    SMLVec3d dY0U = 
      0.5 * ((dA11.xBnd[0].X - dA01.xBnd[0].X) + (dA10.xBnd[0].X - dA00.xBnd[0].X));
    SMLVec3d dY0V = 
      0.5 * ((dA11.xBnd[0].X - dA10.xBnd[0].X) + (dA01.xBnd[0].X - dA00.xBnd[0].X));
    SMLVec3d dY1U = 
      0.5 * ((dA11.xBnd[1].X - dA01.xBnd[1].X) + (dA10.xBnd[1].X - dA00.xBnd[1].X));
    SMLVec3d dY1V = 
      0.5 * ((dA11.xBnd[1].X - dA10.xBnd[1].X) + (dA01.xBnd[1].X - dA00.xBnd[1].X));

    // Compute the scaled normal vectors and derivatives
    SMLVec3d NX  = vnl_cross_3d(XU,  XV);
    SMLVec3d NY0 = vnl_cross_3d(Y0U, Y0V);
    SMLVec3d NY1 = vnl_cross_3d(Y1U, Y1V);

    SMLVec3d dNX  = vnl_cross_3d(dXU,  XV)  + vnl_cross_3d(XU,  dXV);
    SMLVec3d dNY0 = vnl_cross_3d(dY0U, Y0V) + vnl_cross_3d(Y0U, dY0V);
    SMLVec3d dNY1 = vnl_cross_3d(dY1U, Y1V) + vnl_cross_3d(Y1U, dY1V);

    // Compute G and its derivative
    double GX = dot_product(NX, NX);
    double dGX = 2.0 * dot_product(dNX, NX);

    // Compute the Jacobians
    double J0 = dot_product(NY0, NX) / GX;
    double J1 = dot_product(NY1, NX) / GX;
    
    // Compute the derivatives of the Jacobians
    double dJ0 = (dot_product(dNY0, NX) + dot_product(NY0, dNX) - J0 *dGX) / GX;
    double dJ1 = (dot_product(dNY1, NX) + dot_product(NY1, dNX) - J1 *dGX) / GX;

    // Compute the penalty function
    dTotalPenalty += PenaltyFunctionDerivative(J0, 10, 400) * dJ0;
    dTotalPenalty += PenaltyFunctionDerivative(J1, 10, 400) * dJ1;
    }
  delete itQuad;

  // Return the total value
  return dTotalPenalty;
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

/* 
double BoundaryJacobianEnergyTerm::ComputeEnergy(SolutionDataBase *S)
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
    SMLVec3d NM00 = vnl_cross_3d(A01.X - A00.X, A10.X - A00.X);
    SMLVec3d NM11 = vnl_cross_3d(A11.X - A10.X, A11.X - A01.X);
    
    // Compute the jacobians for each boundary side
    for(size_t k = 0; k < 2; k++)
      {
      // Compute area element vectors on the boundary side
      SMLVec3d NB00 = 
        vnl_cross_3d(A01.xBnd[k].X - A00.xBnd[k].X, A10.xBnd[k].X - A00.xBnd[k].X);
      SMLVec3d NB11 = 
        vnl_cross_3d(A11.xBnd[k].X - A10.xBnd[k].X, A11.xBnd[k].X - A01.xBnd[k].X);

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
*/

/*********************************************************************************
 * VOLUME OVERLAP IMAGE MATCH TERM
 ********************************************************************************/
ProbabilisticEnergyTerm::ProbabilisticEnergyTerm(FloatImage *xImage, size_t nCuts)
{
  // Store the inputs
  this->nCuts = nCuts;
  this->xImage = xImage;

  // Compute the image volume (it remains a constant throughout)
  xImageIntegral = xImage->IntegratePositiveVoxels();
}

double ProbabilisticEnergyTerm::ComputeEnergy(SolutionDataBase *data)
{
  // Make sure that the solution data has internal weights and points
  data->UpdateInternalWeights(nCuts);

  // Construct the 'match with the image' function. We simply use the raw
  // image intensity here
  FloatImageEuclideanFunctionAdapter fImage(xImage);

  // Clear the intersection volume accumulator
  xObjectIntegral = 0;

  // double xVolume = 0.0;

  // Iterate over all the points inside the object and take the image value
  MedialInternalPointIterator *it = data->xAtomGrid->NewInternalPointIterator(nCuts);
  for(; !it->IsAtEnd(); ++(*it))
    {
    size_t i = it->GetIndex();
    xObjectIntegral += 
      data->xInternalWeights[i] * fImage.Evaluate(data->xInternalPoints[i]);
    // xVolume += data->xInternalWeights[i];
    }
  delete it;

  // cout << "Bonehead Volume : " << xVolume << endl;

  // Compute the ratio of image integral and maximum possible value
  return xRatio = 1.0 - xObjectIntegral / xImageIntegral;
}

double ProbabilisticEnergyTerm::BeginGradientComputation(SolutionData *S)
{
  // FloatImageTestFunction fImage;
  FloatImageEuclideanFunctionAdapter fImage(xImage);

  // Compute boundary weights
  S->UpdateBoundaryWeights();

  // double xVolume = 0.0;
  
  // Allocate the array of boundary normal vectors, scaled by area element
  xBndNrm = new SMLVec3d[S->xAtomGrid->GetNumberOfBoundaryPoints()];

  // Also, allocate an array of image intensities along the boundary
  xImageVal = new double[S->xAtomGrid->GetNumberOfBoundaryPoints()];

  // Compute the boundary normal vectors
  MedialBoundaryPointIterator *it = S->xAtomGrid->NewBoundaryPointIterator();
  for(; !it->IsAtEnd(); ++(*it))
    {
    size_t i = it->GetIndex();
    BoundaryAtom &bat = GetBoundaryPoint(it, S->xAtoms);
    xImageVal[i] = fImage.Evaluate(bat.X);
    xBndNrm[i] = S->xBoundaryWeights[i] * bat.N;

    // xVolume += dot_product(bat.X, xBndNrm[i]);
    }
  delete it;

  // cout << "Green's Theorem Volume: " << xVolume / 3.0 << endl;

  // Finally, compute the actual volume overlap measure
  return ComputeEnergy(S);
}
  
double ProbabilisticEnergyTerm::ComputePartialDerivative(
  SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // Main thing to compute is the derivative of the intersect region
  double dObjectIntegral = 0.0;

  // Iterate over boundary points, computing Green's Law derivative
  MedialBoundaryPointIterator *it = S->xAtomGrid->NewBoundaryPointIterator();
  for(; !it->IsAtEnd(); ++(*it))
    {
    size_t i = it->GetIndex();
    double NdX = dot_product(xBndNrm[i], GetBoundaryPoint(it, dS->xAtoms).X);
    dObjectIntegral += xImageVal[i] * NdX;
    }
  delete it;

  // Take the derivative of the relative overlap measure
  return - dObjectIntegral / xImageIntegral;
}

void ProbabilisticEnergyTerm::EndGradientComputation()
{
  delete xBndNrm;
  delete xImageVal;
}

void ProbabilisticEnergyTerm::PrintReport(ostream &sout)
{
  sout << "  Probabilistic Energy Term: " << endl;
  sout << "    object integral: " << xObjectIntegral << endl;
  sout << "    image integral : " << xImageIntegral << endl;
  sout << "    ratio          : " << xObjectIntegral / xImageIntegral << endl;
  sout << "    final value    : " << xRatio << endl;
}

/*********************************************************************************
 * Crest Laplacian Penalty Term: Assure negative Rho on the crest
 ********************************************************************************/
double CrestLaplacianEnergyTerm::ComputeEnergy(SolutionDataBase *data)
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
 * Atom Badness Penalty term
 ********************************************************************************/
double AtomBadnessTerm::ComputeEnergy(SolutionDataBase *data)
{
  // Initialize the accumulators
  xMaxBadness = 0;
  xTotalPenalty = xAvgBadness = 0.0;
  nBadAtoms = 0;
  
  // Iterate over all crest atoms
  MedialAtomIterator *itAtom = data->xAtomGrid->NewAtomIterator();
  for( ; !itAtom->IsAtEnd(); ++(*itAtom) )
    {
    // We only care about invalid atoms
    MedialAtom &A = data->xAtoms[itAtom->GetIndex()];
    if(!A.flagValid)
      {
      nBadAtoms++;
      xTotalPenalty += A.xGradRMagSqr - 1.0;
      }
    }
  delete itAtom;

  // Finish up
  nAtoms = data->xAtomGrid->GetNumberOfAtoms();
  xAvgBadness = xTotalPenalty / nBadAtoms;

  // Return the total penalty
  return xTotalPenalty;
}  

double AtomBadnessTerm::ComputePartialDerivative(
  SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // Initialize the accumulators
  double dTotalPenalty = 0.0;
  
  // Iterate over all crest atoms
  MedialAtomIterator *itAtom = S->xAtomGrid->NewAtomIterator();
  for( ; !itAtom->IsAtEnd(); ++(*itAtom) )
    {
    size_t i = itAtom->GetIndex();
    if(!S->xAtoms[i].flagValid)
      dTotalPenalty += dS->xAtoms[i].xGradRMagSqr;
    }
  delete itAtom;

  return dTotalPenalty;
}

void AtomBadnessTerm::PrintReport(ostream &sout)
{
  sout << "  Atom Penalty Term:           " << endl;
  sout << "    number of atoms          : " << nAtoms << endl; 
  sout << "    number of bad atoms      : " << nBadAtoms << endl; 
  sout << "    largest badness value    : " << xMaxBadness << endl; 
  sout << "    average badness value    : " << xAvgBadness << endl; 
  sout << "    total penalty            : " << xTotalPenalty << endl;
}

/*********************************************************************************
 * Angles penalty term
 ********************************************************************************/

// Compute squared cosines of four angles in a quad
double CosineSquareTuple(MedialAtom *A, MedialQuadIterator *it)
{
  // Get the four vectors
  const SMLVec3d &X00 = A[it->GetAtomIndex(0, 0)].X;
  const SMLVec3d &X01 = A[it->GetAtomIndex(0, 1)].X;
  const SMLVec3d &X10 = A[it->GetAtomIndex(1, 0)].X;
  const SMLVec3d &X11 = A[it->GetAtomIndex(1, 1)].X;

  // Compute the differences
  SMLVec3d DX0 = X10 - X00;
  SMLVec3d DX1 = X11 - X01;
  SMLVec3d DY0 = X01 - X00;
  SMLVec3d DY1 = X11 - X10;

  // Compute the lengths squared
  double LX0 = dot_product(DX0, DX0);
  double LX1 = dot_product(DX1, DX1);
  double LY0 = dot_product(DY0, DY0);
  double LY1 = dot_product(DY1, DY1);

  // Compute the cosines squared
  double A00 = dot_product(DX0, DY0);
  double A01 = dot_product(DX1, DY0);
  double A10 = dot_product(DX0, DY1);
  double A11 = dot_product(DX1, DY1);

  // Compute the actual values
  double C00 = (A00 * A00) / (LX0 * LY0);
  double C01 = (A01 * A01) / (LX1 * LY0);
  double C10 = (A10 * A10) / (LX0 * LY1);
  double C11 = (A11 * A11) / (LX1 * LY1);

  // Compute the weighted sum
  return C00+C01+C10+C11;
}

// Compute squared cosines of four angles in a quad
double CosineSquareTupleDerivative(MedialAtom *A, MedialAtom *dA, MedialQuadIterator *it)
{
  size_t i00 = it->GetAtomIndex(0, 0);
  size_t i01 = it->GetAtomIndex(0, 1);
  size_t i10 = it->GetAtomIndex(1, 0);
  size_t i11 = it->GetAtomIndex(1, 1);

  // Compute the differences and their derivatives
  SMLVec3d DX0 = A[i10].X - A[i00].X;
  SMLVec3d DX1 = A[i11].X - A[i01].X;
  SMLVec3d DY0 = A[i01].X - A[i00].X;
  SMLVec3d DY1 = A[i11].X - A[i10].X;

  SMLVec3d dDX0 = dA[i10].X - dA[i00].X;
  SMLVec3d dDX1 = dA[i11].X - dA[i01].X;
  SMLVec3d dDY0 = dA[i01].X - dA[i00].X;
  SMLVec3d dDY1 = dA[i11].X - dA[i10].X;

  // Compute the lengths squared
  double LX0 = dot_product(DX0, DX0);
  double LX1 = dot_product(DX1, DX1);
  double LY0 = dot_product(DY0, DY0);
  double LY1 = dot_product(DY1, DY1);

  double dLX0 = 2.0 * dot_product(DX0, dDX0);
  double dLX1 = 2.0 * dot_product(DX1, dDX1);
  double dLY0 = 2.0 * dot_product(DY0, dDY0);
  double dLY1 = 2.0 * dot_product(DY1, dDY1);

  // Compute the cosines squared and their derivatives
  double A00 = dot_product(DX0, DY0);
  double A01 = dot_product(DX1, DY0);
  double A10 = dot_product(DX0, DY1);
  double A11 = dot_product(DX1, DY1);

  double dA00 = dot_product(dDX0, DY0) + dot_product(DX0, dDY0);
  double dA01 = dot_product(dDX1, DY0) + dot_product(DX1, dDY0);
  double dA10 = dot_product(dDX0, DY1) + dot_product(DX0, dDY1);
  double dA11 = dot_product(dDX1, DY1) + dot_product(DX1, dDY1);

  // Compute the derivatives of the actual values
  double dC00 = ( 2.0 * (A00 * dA00) * (LX0 * LY0) - (A00 * A00) * (LX0 * dLY0 + dLX0 * LY0) ) / ( (LX0 * LY0) * (LX0 * LY0) );
  double dC01 = ( 2.0 * (A01 * dA01) * (LX1 * LY0) - (A01 * A01) * (LX1 * dLY0 + dLX1 * LY0) ) / ( (LX1 * LY0) * (LX1 * LY0) );
  double dC10 = ( 2.0 * (A10 * dA10) * (LX0 * LY1) - (A10 * A10) * (LX0 * dLY1 + dLX0 * LY1) ) / ( (LX0 * LY1) * (LX0 * LY1) );
  double dC11 = ( 2.0 * (A11 * dA11) * (LX1 * LY1) - (A11 * A11) * (LX1 * dLY1 + dLX1 * LY1) ) / ( (LX1 * LY1) * (LX1 * LY1) );

  // Compute the weighted sum
  return dC00 + dC01 + dC10 + dC11;
}

double MedialAnglesPenaltyTerm::ComputeEnergy(SolutionDataBase *data)
{
  xTotalPenalty = 0.0;

  // Iterate over all the quads in the medial mesh
  MedialQuadIterator *itQuad = data->xAtomGrid->NewQuadIterator();
  for( ; !itQuad->IsAtEnd(); ++(*itQuad))
    {
    // Compute the cosine square of each angle
    xTotalPenalty += CosineSquareTuple(data->xAtoms, itQuad);
    }
  delete itQuad;

  return 0.25 * xTotalPenalty / data->xAtomGrid->GetNumberOfAtoms();
}

double MedialAnglesPenaltyTerm::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  double dTotalPenalty = 0.0;

  // Iterate over all the quads in the medial mesh
  MedialQuadIterator *itQuad = S->xAtomGrid->NewQuadIterator();
  for( ; !itQuad->IsAtEnd(); ++(*itQuad))
    {
    // Compute the cosine square of each angle
    dTotalPenalty += CosineSquareTupleDerivative(S->xAtoms, dS->xAtoms, itQuad);
    }
  delete itQuad;

  return 0.25 * dTotalPenalty / S->xAtomGrid->GetNumberOfAtoms();
}

void MedialAnglesPenaltyTerm::PrintReport(ostream &sout)
{
  sout << "  Square Cosine Penalty Term : " << endl;
  sout << "    total penalty            : " << xTotalPenalty << endl;
}


/*********************************************************************************
 * Regularity term (used to maintain correspondence)
 ********************************************************************************/
MedialRegularityTerm::MedialRegularityTerm(
  MedialAtomGrid *inAtomGrid, MedialAtom *inAtomArray)
{
  // Generate an array of weights for irregular grids
  xDomainWeights.set_size(inAtomGrid->GetNumberOfAtoms());
  xDomainArea = ComputeMedialDomainAreaWeights(
    inAtomGrid, inAtomArray, xDomainWeights.data_block());
  cout << "Domain Area: " << xDomainArea << endl;
}


double MedialRegularityTerm::ComputeEnergy(SolutionDataBase *data)
{
  // Integral of (G2(1,11) + G(2,12))^2 + (G(1,21) + G(2,22)) du dv 
  // (no area elt. scaling)
  xGradMagIntegral = 0;
  xGradMagMaximum = 0;
  
  // Integrate over all atoms
  MedialAtomIterator *it = data->xAtomGrid->NewAtomIterator();
  for(; !it->IsAtEnd(); ++(*it))
    {
    // Get the medial atom
    MedialAtom &a = data->xAtoms[it->GetIndex()];

    // Only compute the penalty if the atom is internal
    if(a.flagCrest || !a.flagValid) 
      continue;

    // Compute the regularity term here
    double reg1 = a.G.xChristoffelSecond[0][0][0] + a.G.xChristoffelSecond[1][0][1];
    double reg2 = a.G.xChristoffelSecond[0][1][0] + a.G.xChristoffelSecond[1][1][1];
    double reg = reg1 * reg1 + reg2 * reg2;

    // Update the sum and the maximum
    if(reg > xGradMagMaximum)
      xGradMagMaximum = reg;

    // Update the integral
    xGradMagIntegral += xDomainWeights[it->GetIndex()] * reg;
    }

  delete it;

  return xGradMagIntegral / xDomainArea;
}

double MedialRegularityTerm
::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  double dIntegral = 0.0;
  
  // Integrate over all atoms
  MedialAtomIterator *it = S->xAtomGrid->NewAtomIterator();
  for(; !it->IsAtEnd(); ++(*it))
    {
    // Get the four atoms at the corner of the quad
    MedialAtom &a = S->xAtoms[it->GetIndex()];
    MedialAtom &da = dS->xAtoms[it->GetIndex()];

    // Only compute the penalty if the atom is internal
    if(a.flagCrest || !a.flagValid) 
      continue;

    // Compute the regularity term here
    double reg1 = a.G.xChristoffelSecond[0][0][0] + a.G.xChristoffelSecond[1][0][1];
    double reg2 = a.G.xChristoffelSecond[0][1][0] + a.G.xChristoffelSecond[1][1][1];
    double dreg1 = da.G.xChristoffelSecond[0][0][0] + da.G.xChristoffelSecond[1][0][1];
    double dreg2 = da.G.xChristoffelSecond[0][1][0] + da.G.xChristoffelSecond[1][1][1];
    double dreg = 2.0 * (reg1 * dreg1 + reg2 * dreg2);

    // Update the integral
    dIntegral += xDomainWeights[it->GetIndex()] * dreg;
    }

  delete it;

  // Return the error
  return dIntegral / xDomainArea; 
}

void MedialRegularityTerm::PrintReport(ostream &sout)
{
  sout << "  Medial Regularization Term: " << endl;
  sout << "    maximum grad mag sqr:     " << xGradMagMaximum << endl;
  sout << "    grad mag sqr integral: " << xGradMagIntegral << endl;
  sout << "    domain area: " << xDomainArea << endl;
}


/*********************************************************************************
 * MEDIAL RADIUS PENALTY TERM
 ********************************************************************************/
double RadiusPenaltyTerm::ComputeEnergy(SolutionDataBase *S)
{ 
  xTotalPenalty = 0.0;
  xMinR2 = 1e100;

  for(size_t i = 0; i < S->xAtomGrid->GetNumberOfAtoms(); i++)
    {
    // Get the square of R
    double phi = S->xAtoms[i].F;
    if(phi < xMinR2)
      xMinR2 = phi;

    // Apply the penalty function
    double x = phi / xScale;
    double x2 = x * x;
    double x4 = x2 * x2;
    
    xTotalPenalty += 1.0 / x4;
    }

  xTotalPenalty /= S->xAtomGrid->GetNumberOfAtoms();

  return xTotalPenalty;
}

double RadiusPenaltyTerm
::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  double dTotalPenalty = 0.0;

  for(size_t i = 0; i < S->xAtomGrid->GetNumberOfAtoms(); i++)
    {
    // Get the square of R
    double phi = S->xAtoms[i].F;
    double dphi = dS->xAtoms[i].F;

    // Apply the penalty function
    double x = phi / xScale, dx = dphi / xScale;
    double x2 = x * x;
    double x4 = x2 * x2;

    dTotalPenalty += -4.0 * dx / (x * x4);
    }

  dTotalPenalty /= S->xAtomGrid->GetNumberOfAtoms();

  return dTotalPenalty;
}  

void RadiusPenaltyTerm::PrintReport(ostream &sout)
{
  sout << "  Radius Penalty Term : " << endl;
  sout << "    total penalty            : " << xTotalPenalty << endl;
  sout << "    smallest R^2             : " << xMinR2 << endl;
}


/*********************************************************************************
 * MEDIAL OPTIMIZATION PROBLEM
 ********************************************************************************/
const double MedialOptimizationProblem::xPrecision = 1.0e-10;
const double MedialOptimizationProblem::xEpsilon = 1.0e-6;

MedialOptimizationProblem
::MedialOptimizationProblem(
  MedialPDESolver *xSolver, IMedialCoefficientMask *xCoeff)
{
  this->xCoeff = xCoeff;
  this->xSolver = xSolver;
  this->nCoeff = xCoeff->GetNumberOfCoefficients();

  flagLastEvalAvailable = false;
  flagPhiGuessAvailable = false;
  
  // Prepare medial atoms for gradient computation
  for(size_t i = 0; i < nCoeff; i++)
    {
    MedialAtom *a = new MedialAtom[xSolver->GetAtomGrid()->GetNumberOfAtoms()];
    IHyperSurface2D *xVariation = xCoeff->GetComponentSurface(i);
    xSolver->PrepareAtomsForVariationalDerivative(xVariation, a);
    xCoeff->ReleaseComponentSurface(xVariation);
    dAtomArray.push_back(a);
    }
}

MedialOptimizationProblem
::~MedialOptimizationProblem()
{
  for(size_t i = 0; i < nCoeff; i++)
    delete dAtomArray[i];
}

// Add the energy term
void MedialOptimizationProblem::AddEnergyTerm(EnergyTerm *term, double xWeight)
{
  xWeights.push_back(xWeight);
  xTerms.push_back(term);
  xTimers.push_back(CodeTimer());
  xGradTimers.push_back(CodeTimer());
}


bool MedialOptimizationProblem::SolvePDE(double *xEvalPoint)
{
  // Make a vector from the eval point
  vnl_vector<double> X(xEvalPoint, nCoeff);

  // This flag will tell us if the last evaluation failed
  bool flagSolveOk;

  // Check if we are re-evaluating at a previously tried location
  if(flagLastEvalAvailable && X == xLastEvalPoint)
    return false;

  // Solve the equation
  xSolveTimer.Start();

  // Update the surface with the new solution
  xCoeff->SetCoefficientArray(xEvalPoint);

  // Solve the PDE using the last phi field if we have it
  // flagSolveOk = (flagPhiGuessAvailable)
  //   ? xSolver->Solve(xLastPhiField, xPrecision)
  //   : xSolver->Solve(xPrecision);
  flagSolveOk = xSolver->Solve(xPrecision);

  // Stop timing
  xSolveTimer.Stop();

  // Set the last evaluation point (this is to ensure we never double-evaluate)
  xLastEvalPoint = X;
  flagLastEvalAvailable = true;

  // Return true : solution updated
  return true;
}

double MedialOptimizationProblem::Evaluate(double *xEvalPoint)
{
  // Solve the PDE - if there is no update, return the last solution value
  if( SolvePDE(xEvalPoint) )
    {
    // Create a solution data object representing current solution
    SolutionData S0(xSolver);

    // Compute the solution for each term
    xLastSolutionValue = 0.0;
    for(size_t iTerm = 0; iTerm < xTerms.size(); iTerm++)
      { 
      xTimers[iTerm].Start();
      xLastSolutionValue += xTerms[iTerm]->ComputeEnergy(&S0) * xWeights[iTerm]; 
      xTimers[iTerm].Stop();
      }
    }

  // Return the result
  return xLastSolutionValue;
}

/*
double MedialOptimizationProblem
::ComputePartialDerivative(double *xEvalPoint, double *xDirection, double &df)
{
  // Solve the PDE at the current point
  SolvePDE(xEvalPoint);

  // Initialize the problem with the current solution
  xSolver->SetSolutionAsInitialGuess();

  // Create a solution data 
  SolutionData S(xSolver);

  // Compute the variation in the gradient direction and the var. deriv.
  IHyperSurface2D *xVariation = xCoeff->GetVariationSurface(xDirection);
  xSolveGradTimer.Start();
  xSolver->ComputeVariationalDerivative(xVariation, dAtoms); 
  xSolveGradTimer.Stop();
  PartialDerivativeSolutionData dS(&S, dAtoms);

  // Compute the solution at the current point
  xLastSolutionValue = 0.0; df = 0.0;
  for(size_t iTerm = 0; iTerm < xTerms.size(); iTerm++)
    {
    // Compute the solution
    xTimers[iTerm].Start();
    xLastSolutionValue +=
      xWeights[iTerm] * xTerms[iTerm]->BeginGradientComputation(&S);
    xTimers[iTerm].Stop();

    // Compute the partial derivative
    xSolveGradTimer.Start();
    df += xWeights[iTerm] * xTerms[iTerm]->ComputePartialDerivative(&S, &dS);
    xSolveGradTimer.Stop();

    // Clean up
    xTerms[iTerm]->EndGradientComputation();
    }

  // Release the variation surface
  xCoeff->ReleaseVariationSurface(xVariation);

  // Return the function value
  return xLastSolutionValue;
}
*/

double 
MedialOptimizationProblem
::ComputeGradient(double *xEvalPoint, double *xGradient)
{
  size_t iTerm;

  // Solve the PDE
  SolvePDE(xEvalPoint);

  // Compute the gradient
  xSolveGradTimer.Start();
  xSolver->ComputeVariationalGradient(dAtomArray); 
  xSolveGradTimer.Stop();

  // Store this solution as the new initialization for the PDE solver
  xSolver->SetSolutionAsInitialGuess();

  // Create a solution data object representing current solution
  SolutionData S(xSolver);

  // Compute the value at the solution and init gradient computation
  xLastSolutionValue = 0.0;
  for(iTerm = 0; iTerm < xTerms.size(); iTerm++)
    {
    xTimers[iTerm].Start();
    xLastSolutionValue += 
      xWeights[iTerm] * xTerms[iTerm]->BeginGradientComputation(&S);
    xTimers[iTerm].Stop();
    }
    
  // Repeat for each coefficient
  for(size_t iCoeff = 0; iCoeff < nCoeff; iCoeff++)
    {
    // Create a solution partial derivative object
    PartialDerivativeSolutionData dS(&S, dAtomArray[iCoeff]);
    
    // Compute the partial derivatives for each term
    xGradient[iCoeff] = 0.0;
    for(iTerm = 0; iTerm < xTerms.size(); iTerm++)
      {
      xGradTimers[iTerm].Start();
      xGradient[iCoeff] += xWeights[iTerm] *
        xTerms[iTerm]->ComputePartialDerivative(&S, &dS);
      xGradTimers[iTerm].Stop();
      }

    // Dump the gradient
    // cout << iCoeff << "; " << XGradient[iCoeff] << endl;
    cout << "." << flush;
    }

  cout << endl;

  // Clear up gradient computation
  for(iTerm = 0; iTerm < xTerms.size(); iTerm++)
    xTerms[iTerm]->EndGradientComputation();

  // Increment the evaluation counter
  evaluationCost += nCoeff;

  // Return the solution value
  return xLastSolutionValue;
}

void MedialOptimizationProblem::PrintReport(ostream &sout)
{
  sout << "Optimization Summary: " << endl;
  sout << "  energy value   :  " << xLastSolutionValue << endl; 
  sout << "  evaluation cost: " << this->getEvaluationCost() << endl;
  sout << "  solver time : " << xSolveTimer.Read() << endl;
  sout << "  solver gradient time : " << xSolveGradTimer.Read() << endl;
  sout << "  weights time : " << xWeightsTimer.Read() << endl;
  sout << "  weights gradient time : " << xWeightsGradTimer.Read() << endl;
  sout << "Per-Term Report:" << endl;

  for(size_t iTerm = 0; iTerm < xTerms.size(); iTerm++)
    {
    sout << "Term " << iTerm << endl;
    xTerms[iTerm]->PrintReport(sout);
    sout << "    elapsed time: " << xTimers[iTerm].Read() << endl;
    sout << "    gradient time: " << xGradTimers[iTerm].Read() << endl;
    }
}

