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
  xFinalMatch = sqrt( xImageMatch / xBoundaryArea );

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
  xFinalMatch = sqrt( xImageMatch / xBoundaryArea );
  
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

  // Account for the square root
  dFinaldC *= 0.5 / xFinalMatch;

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
    xTotalPenalty += PenaltyFunction(J0, 10, 40);
    xTotalPenalty += PenaltyFunction(J1, 10, 40);
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
    dTotalPenalty += PenaltyFunctionDerivative(J0, 10, 40) * dJ0;
    dTotalPenalty += PenaltyFunctionDerivative(J1, 10, 40) * dJ1;
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


/*
double VolumeOverlapEnergyTerm::BeginGradientComputation(SolutionData *S)
{
  // Clear the intersection volume accumulator
  xIntersect = 0;

  // Make sure that the solution data has internal weights and points
  S->UpdateInternalWeights(nCuts);

  // Construct the 'match with the image' function
  // FloatImageVolumeElementFunctionAdapter fImage(xImage);
  FloatImageTestFunction fImage;

  // Initialize the array of image values and gradients
  xImageVal = new double[S->xAtomGrid->GetNumberOfInternalPoints(nCuts)];
  xImageGrad = new SMLVec3d[S->xAtomGrid->GetNumberOfInternalPoints(nCuts)];

  // Iterate over all the points inside the object and take the image value
  MedialInternalPointIterator *it = S->xAtomGrid->NewInternalPointIterator(nCuts);
  for(; !it->IsAtEnd(); ++(*it))
    {
    size_t iPoint = it->GetIndex();

    // Take the image value and store it for later
    xImageVal[iPoint] = fImage.Evaluate(S->xInternalPoints[iPoint]);

    // Compute the image gradient as well
    xImageGrad[iPoint] = fImage.ComputeGradient(S->xInternalPoints[iPoint]);

    // Compute the weight and the weighted image value
    double w = S->xInternalWeights[iPoint];
    xIntersect += xImageVal[iPoint] * w; 
    }
  delete it;

  // Get (remember) the object volume 
  xObjectVolume = S->xInternalVolume;

  // Compute the volume overlap measure
  xUnion = xImageVolume + xObjectVolume - xIntersect;
  xOverlap = xIntersect / xUnion;

  // Return a value that can be minimized
  return 1.0 - xOverlap;
}

double VolumeOverlapEnergyTerm::ComputePartialDerivative(
  SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // Make sure the weights are computed in the derivative data
  dS->UpdateInternalWeights(nCuts);

  // Compute the derivative of the intersection volume
  double dIntersect = 0.0;

  for(size_t i = 0; i < S->xAtomGrid->GetNumberOfInternalPoints(nCuts); i++)
    {
    // Simple chain rule applies
    dIntersect += dS->xInternalWeights[i] * xImageVal[i] + 
      S->xInternalWeights[i] * dot_product(dS->xInternalPoints[i], xImageGrad[i]);
    }

  // Now compute the derivative of the overlap volume
  double dObjectVolume = dS->xInternalVolume;
  double dUnion = dObjectVolume - dIntersect;
  double dOverlap = 
    ( dIntersect * xUnion - xIntersect * dUnion ) / (xUnion * xUnion );

  // Done!
  return -dOverlap;
}

void VolumeOverlapEnergyTerm::EndGradientComputation()
{
  delete xImageVal;
  delete xImageGrad;
}
*/
  
  
  /*
double VolumeOverlapEnergyTerm::BeginGradientComputation(SolutionData *S)
{
  // Get info about the problem
  MedialAtomGrid *xGrid = S->xAtomGrid;
  size_t nPoints = xGrid->GetNumberOfBoundaryPoints();

  // Create the image adapter for image interpolation at the boundary
  FloatImageVolumeElementFunctionAdapter fImage(xImage);

  // We require the weights from the central solution
  S->UpdateBoundaryWeights();

  // Interpolate the image at each sample point on the surface
  xImageVal = new double[nPoints];
  MedialBoundaryPointIterator *itBnd = xGrid->NewBoundaryPointIterator();
  for( ; !itBnd->IsAtEnd() ; ++(*itBnd) )
    {
    // Compute and store the image value at the sample point
    SMLVec3d &X = GetBoundaryPoint(itBnd, S->xAtoms).X;
    xImageVal[itBnd->GetIndex()] = 1.0 + 0.5 * fImage.Evaluate(X);
    }
  delete itBnd;
  
  // Return the raw solution
  return ComputeEnergy(S);
}

double
VolumeOverlapEnergyTerm
::ComputePartialDerivative(
    SolutionData *SCenter, SolutionData *SFwd, SolutionData *SBck,
    double xEpsilon)
{
  // Get info about the problem
  MedialAtomGrid *xGrid = SCenter->xAtomGrid;
  size_t nPoints = xGrid->GetNumberOfBoundaryPoints();

  // We will compute the partial derivative of the absolute overlap and the
  // partial derivative of the m-rep volume with respect to the coeff
  double dOverlap = 0.0, dVolume = 0.0;

  // Compute the partial derivative for this coefficient
  MedialBoundaryPointIterator *itBnd = xGrid->NewBoundaryPointIterator();
  
  for( ; !itBnd->IsAtEnd(); ++(*itBnd))
    {
    // Get the index of this boundary point
    size_t iPoint = itBnd->GetIndex();

    // Get the boundary point before and after small deformation
    BoundaryAtom &B = GetBoundaryPoint(itBnd, SCenter->xAtoms);
    SMLVec3d X1 = GetBoundaryPoint(itBnd, SBck->xAtoms).X;
    SMLVec3d X2 = GetBoundaryPoint(itBnd, SFwd->xAtoms).X;
    SMLVec3d dX = 0.5 * (X2 - X1);

    // Get the area weights for this point
    double w = SCenter->xBoundaryWeights[ iPoint ];

    // Compute the common term between the two
    double z = dot_product( dX, B.N ) * w;

    // Compute the two partials
    dOverlap += z * xImageVal[iPoint];
    dVolume += z;
    }
  
  delete itBnd;

  // Scale the partials by epsilon to get the derivative
  dOverlap /= xEpsilon; dVolume /= xEpsilon;

  // Compute the partial derivative of relative overlap match
  double xUnion = (xImageVolume + xObjectVolume - xIntersect);        
  double xPartial = 
    (dVolume * xIntersect - dOverlap * (xImageVolume + xObjectVolume));
  xPartial /= xUnion * xUnion;

  // Return the partial derivative
  return xPartial;
}

void VolumeOverlapEnergyTerm::EndGradientComputation()
{
  delete xImageVal;
}
*/
/*
void VolumeOverlapEnergyTerm::PrintReport(ostream &sout)
{
  sout << "  Volume Overlap Term: " << endl;
  sout << "    image volume   : " << xImageVolume << endl;
  sout << "    object volume  : " << xObjectVolume << endl;
  sout << "    intersection   : " << xIntersect << endl;
  sout << "    overlap ratio  : " << xOverlap << endl;
}
*/
  
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
 * Regularity term (used to maintain correspondence)
 ********************************************************************************/
MedialRegularityTerm::MedialRegularityTerm(
  MedialAtom *xTempAtoms, MedialAtomGrid *xTempGrid)
{
  // Allocate appropriate space
  xEdgeLength.reserve(xTempGrid->GetNumberOfQuads() * 6);

  // Compute all the distances between adjacent atoms. This involves quad diagonals 
  // as well as sides, because a quad can be deformed without changing edge lengths, 
  // while a triangle is defined by the edge lengths.
  MedialQuadIterator *itQuad = xTempGrid->NewQuadIterator();
  for(; !itQuad->IsAtEnd(); ++(*itQuad))
    {
    // Get the four atoms at the corner of the quad
    MedialAtom &a00 = xTempAtoms[itQuad->GetAtomIndex(0, 0)];
    MedialAtom &a01 = xTempAtoms[itQuad->GetAtomIndex(0, 1)];
    MedialAtom &a10 = xTempAtoms[itQuad->GetAtomIndex(1, 0)];
    MedialAtom &a11 = xTempAtoms[itQuad->GetAtomIndex(1, 1)];

    // Record the distances (repeat in this order later)
    xEdgeLength.push_back( (a00.X - a01.X).two_norm() );
    xEdgeLength.push_back( (a01.X - a11.X).two_norm() );
    xEdgeLength.push_back( (a11.X - a10.X).two_norm() );
    xEdgeLength.push_back( (a10.X - a00.X).two_norm() );
    xEdgeLength.push_back( (a00.X - a11.X).two_norm() );
    xEdgeLength.push_back( (a01.X - a10.X).two_norm() );
    }

  delete itQuad;
}

double MedialRegularityTerm::DistortionPenalty(
  double dTemplate, double dSubject)
{
  double xDistortion = dSubject - dTemplate;
  double xPenalty = xDistortion * xDistortion;
  
  // Record the maximum and minimum distortion
  if(xDistortion < xMinDistortion) xMinDistortion = xDistortion;
  if(xDistortion < xMaxDistortion) xMaxDistortion = xDistortion;

  // Return the squared distortion
  return xPenalty;
}

double MedialRegularityTerm::DistortionPenaltyAndDerivative(
  double dTemplate, const SMLVec3d &A, const SMLVec3d &B, double &xGradTerm)
{
  double dSubject = (A - B).two_norm();
  double xDistortion = dSubject - dTemplate;
  double xPenalty = xDistortion * xDistortion;
  
  // Record the maximum and minimum distortion
  if(xDistortion < xMinDistortion) xMinDistortion = xDistortion;
  if(xDistortion < xMaxDistortion) xMaxDistortion = xDistortion;

  // Compute the derivative term
  xGradTerm = 2.0 * xDistortion / dSubject;
  
  // Return the squared distortion
  return xPenalty;
}



double MedialRegularityTerm::ComputeEnergy(SolutionDataBase *data)
{
  // The prior is based on the distortion in terms of distances between adjacent 
  // atoms as compared to the template. We iterate over the edges on the medial 
  // surface and compute these distortions
  xTotalPenalty = 0.0;
  xMaxDistortion = xMinDistortion = 0.0;
  
  // Iterate over the quads
  MedialQuadIterator *itQuad = data->xAtomGrid->NewQuadIterator();
  for(size_t i = 0; !itQuad->IsAtEnd(); ++(*itQuad), i+=6)
    {
    // Get the four atoms at the corner of the quad
    MedialAtom &a00 = data->xAtoms[itQuad->GetAtomIndex(0, 0)];
    MedialAtom &a01 = data->xAtoms[itQuad->GetAtomIndex(0, 1)];
    MedialAtom &a10 = data->xAtoms[itQuad->GetAtomIndex(1, 0)];
    MedialAtom &a11 = data->xAtoms[itQuad->GetAtomIndex(1, 1)];

    // Record the distances (repeat in this order later)
    xTotalPenalty += DistortionPenalty(xEdgeLength[i]  , (a00.X - a01.X).two_norm());
    xTotalPenalty += DistortionPenalty(xEdgeLength[i+1], (a01.X - a11.X).two_norm());
    xTotalPenalty += DistortionPenalty(xEdgeLength[i+2], (a11.X - a10.X).two_norm());
    xTotalPenalty += DistortionPenalty(xEdgeLength[i+3], (a10.X - a00.X).two_norm());
    xTotalPenalty += DistortionPenalty(xEdgeLength[i+4], (a00.X - a11.X).two_norm());
    xTotalPenalty += DistortionPenalty(xEdgeLength[i+5], (a01.X - a10.X).two_norm());
    }

  // Compute the average error
  xTotalPenalty /= xEdgeLength.size();
  xMeanSquareDistortion = xTotalPenalty;

  // Return the error
  return xTotalPenalty; 
}

double MedialRegularityTerm::BeginGradientComputation(SolutionData *data)
{
  // The prior is based on the distortion in terms of distances between adjacent 
  // atoms as compared to the template. We iterate over the edges on the medial 
  // surface and compute these distortions
  xTotalPenalty = 0.0;
  xMaxDistortion = xMinDistortion = 0.0;

  // Allocate an array to hold the derivative of the distortion function
  xGradientCommon.resize( xEdgeLength.size(), 0.0 );
  
  // Iterate over the quads
  MedialQuadIterator *itQuad = data->xAtomGrid->NewQuadIterator();
  for(size_t i = 0; !itQuad->IsAtEnd(); ++(*itQuad), i+=6)
    {
    // Get the four atoms at the corner of the quad
    MedialAtom &a00 = data->xAtoms[itQuad->GetAtomIndex(0, 0)];
    MedialAtom &a01 = data->xAtoms[itQuad->GetAtomIndex(0, 1)];
    MedialAtom &a10 = data->xAtoms[itQuad->GetAtomIndex(1, 0)];
    MedialAtom &a11 = data->xAtoms[itQuad->GetAtomIndex(1, 1)];

    // Record the distances (repeat in this order later)
    xTotalPenalty += DistortionPenaltyAndDerivative(
      xEdgeLength[i+0], a00.X, a01.X, xGradientCommon[i+0]);
    xTotalPenalty += DistortionPenaltyAndDerivative(
      xEdgeLength[i+1], a01.X, a11.X, xGradientCommon[i+1]);
    xTotalPenalty += DistortionPenaltyAndDerivative(
      xEdgeLength[i+2], a11.X, a10.X, xGradientCommon[i+2]);
    xTotalPenalty += DistortionPenaltyAndDerivative(
      xEdgeLength[i+3], a10.X, a00.X, xGradientCommon[i+3]);
    xTotalPenalty += DistortionPenaltyAndDerivative(
      xEdgeLength[i+4], a00.X, a11.X, xGradientCommon[i+4]);
    xTotalPenalty += DistortionPenaltyAndDerivative(
      xEdgeLength[i+5], a01.X, a10.X, xGradientCommon[i+5]);
    }

  // Compute the average error
  xTotalPenalty /= xEdgeLength.size();
  xMeanSquareDistortion = xTotalPenalty;

  // Return the error
  return xTotalPenalty; 
}

double MedialRegularityTerm
::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // The prior is based on the distortion in terms of distances between adjacent 
  // atoms as compared to the template. We iterate over the edges on the medial 
  // surface and compute these distortions
  double dTotalPenalty = 0.0;
  
  // Iterate over the quads
  MedialQuadIterator *itQuad = S->xAtomGrid->NewQuadIterator();
  for(size_t i = 0; !itQuad->IsAtEnd(); ++(*itQuad), i+=6)
    {
    // Get the four atoms at the corner of the quad
    MedialAtom &a00 = S->xAtoms[itQuad->GetAtomIndex(0, 0)];
    MedialAtom &a01 = S->xAtoms[itQuad->GetAtomIndex(0, 1)];
    MedialAtom &a10 = S->xAtoms[itQuad->GetAtomIndex(1, 0)];
    MedialAtom &a11 = S->xAtoms[itQuad->GetAtomIndex(1, 1)];

    // Get the four atom derivatives at the corner of the quad
    MedialAtom &d00 = dS->xAtoms[itQuad->GetAtomIndex(0, 0)];
    MedialAtom &d01 = dS->xAtoms[itQuad->GetAtomIndex(0, 1)];
    MedialAtom &d10 = dS->xAtoms[itQuad->GetAtomIndex(1, 0)];
    MedialAtom &d11 = dS->xAtoms[itQuad->GetAtomIndex(1, 1)];

    // Record the distances (repeat in this order later)
    dTotalPenalty += xGradientCommon[i+0] * dot_product(a00.X - a01.X, d00.X - d01.X);
    dTotalPenalty += xGradientCommon[i+1] * dot_product(a01.X - a11.X, d01.X - d11.X);
    dTotalPenalty += xGradientCommon[i+2] * dot_product(a11.X - a10.X, d11.X - d10.X);
    dTotalPenalty += xGradientCommon[i+3] * dot_product(a10.X - a00.X, d10.X - d00.X);
    dTotalPenalty += xGradientCommon[i+4] * dot_product(a00.X - a11.X, d00.X - d11.X);
    dTotalPenalty += xGradientCommon[i+5] * dot_product(a01.X - a10.X, d01.X - d10.X);
    }

  // Compute the average error
  dTotalPenalty /= xEdgeLength.size();

  // Return the error
  return dTotalPenalty; 
}

void MedialRegularityTerm::PrintReport(ostream &sout)
{
  sout << "  Medial Regularization Term: " << endl;
  sout << "    minimum distortion:     " << xMinDistortion << endl;
  sout << "    maximum distortion:     " << xMaxDistortion << endl;
  sout << "    mean square distortion: " << xMeanSquareDistortion << endl;
  sout << "    total penalty:          " << xTotalPenalty << endl;
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
  SolutionData S0(xSolver);

  // cout << "REPORT:" << endl;

  // Compute the solution for each term
  xSolution = 0.0;
  for(size_t iTerm = 0; iTerm < xTerms.size(); iTerm++)
    { xSolution += xTerms[iTerm]->ComputeEnergy(&S0) * xWeights[iTerm]; }

  // cout << "RESULT:" << xSolution << endl;
  // PrintReport(cout);

  // Return the result
  return xSolution;
}

double 
MedialOptimizationProblem
::ComputeGradient(double *X, double *XGradient)
{
  size_t iTerm;

  // Copy the x-vector into the coefficients of the surface
  size_t nCoeff = xCoeff->GetNumberOfCoefficients();
  xCoeff->SetCoefficientArray(X);

  // Solve using default initial guess
  // xSolver->SetDefaultInitialGuess();
  xSolver->Solve(xPrecision);

  // Create a solution data object representing current solution
  SolutionData S(xSolver);

  // Compute the value at the solution and init gradient computation
  xSolution = 0.0;
  for(iTerm = 0; iTerm < xTerms.size(); iTerm++)
    xSolution += 
      xWeights[iTerm] * xTerms[iTerm]->BeginGradientComputation(&S);

  // Allocate an array of derivative atoms
  MedialAtom *dAtoms = new MedialAtom[xSolver->GetAtomGrid()->GetNumberOfAtoms()];

  // Repeat for each coefficient
  for(size_t iCoeff = 0; iCoeff < nCoeff; iCoeff++)
    {
    // Get an adapter for evaluating the function
    IHyperSurface2D *xVariation = xCoeff->GetComponentSurface(iCoeff);

    // Compute the partial derivatives with respect to the coefficients
    xSolver->ComputeVariationalDerivative(xVariation, dAtoms); 

    // Create a solution partial derivative object
    PartialDerivativeSolutionData dS(&S, dAtoms);
    
    // Compute the partial derivatives for each term
    XGradient[iCoeff] = 0.0;
    for(iTerm = 0; iTerm < xTerms.size(); iTerm++)
      XGradient[iCoeff] += xWeights[iTerm] *
        xTerms[iTerm]->ComputePartialDerivative(&S, &dS);

    // Dump the gradient
    // cout << iCoeff << "; " << XGradient[iCoeff] << endl;
    cout << "." << flush;

    // Release the component surface
    xCoeff->ReleaseComponentSurface(xVariation);
    }

  // Clear up gradient computation
  for(iTerm = 0; iTerm < xTerms.size(); iTerm++)
    xTerms[iTerm]->EndGradientComputation();

  // Delete the derivative atoms array
  delete dAtoms;

  // Increment the evaluation counter
  cout << endl;
  evaluationCost += nCoeff;

  // Return the solution value
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

