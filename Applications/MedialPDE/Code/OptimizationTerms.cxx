#include "OptimizationTerms.h"
#include "MedialAtomGrid.h"
#include "MedialAtom.h"
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
SolutionData::SolutionData( MedialIterationContext *xGrid, MedialAtom *xAtoms) :
  SolutionDataBase()
{
  this->xAtomGrid = xGrid;
  this->xAtoms = xAtoms;

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

  // Loop over all boundary points
  for(MedialBoundaryPointIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    // Compute the image gradient
    SMLVec3d &X = GetBoundaryPoint(it, S->xAtoms).X;
    fImage.ComputeGradient(X, xGradI[it.GetIndex()]);

    // Compute the image value
    xImageVal[it.GetIndex()] = fImage.Evaluate(X);

    // Accumulate to get weighted match
    xImageMatch += xImageVal[it.GetIndex()] * S->xBoundaryWeights[it.GetIndex()];
    }
  
  // We will need the area in many calculations
  xBoundaryArea = S->xBoundaryArea;

  // Compute the final match
  xFinalMatch = xImageMatch / xBoundaryArea;
  
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
  
  for(MedialBoundaryPointIterator it(S->xAtomGrid) ; !it.IsAtEnd(); ++it)
    {
    // Get the index of this boundary point
    size_t iPoint = it.GetIndex();

    // Get the change in the boundary point
    SMLVec3d dX = GetBoundaryPoint(it, DS->xAtoms).X;

    // Get the area weights for this point
    double w = S->xBoundaryWeights[ iPoint ];
    double dw = DS->xBoundaryWeights[ iPoint ];

    // Compute the change in intensity per change in coefficient
    double I = xImageVal[iPoint];
    double dIdC = dot_product( xGradI[iPoint] , dX);
    
    // Increment the partial derivative of the weighted match
    dMatchdC += dw * I + w * dIdC;
    }
  
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

  // Reset the TriangleEntry vector
  if(xTriangleEntries.size() != S->xAtomGrid->GetNumberOfTriangles())
    xTriangleEntries.resize(S->xAtomGrid->GetNumberOfTriangles());

  // Keep track of the entry index
  TriangleVector::iterator eit = xTriangleEntries.begin();
  
  // Create a Triangle-iterator through the atoms
  for(MedialTriangleIterator itt(S->xAtomGrid); !itt.IsAtEnd(); ++itt, ++eit)
    {
    // Get the four atoms in this quad
    MedialAtom &A0 = S->xAtoms[itt.GetAtomIndex(0)];
    MedialAtom &A1 = S->xAtoms[itt.GetAtomIndex(1)];
    MedialAtom &A2 = S->xAtoms[itt.GetAtomIndex(2)];

    // Compute the Xu and Xv vectors and the normal
    eit->XU = A1.X - A0.X; eit->XV = A2.X - A0.X;
    eit->NX = vnl_cross_3d(eit->XU, eit->XV);
    eit->gX2 = dot_product(eit->NX, eit->NX);

    // Compute side-wise entries
    for(size_t z = 0; z < 2; z++)
      {
      // Compute the same for the upper and lower boundaries
      eit->YU[z] = A1.xBnd[z].X - A0.xBnd[z].X;
      eit->YV[z] = A2.xBnd[z].X - A0.xBnd[z].X;
      eit->NY[z] = vnl_cross_3d(eit->YU[z], eit->YV[z]);
      
      // Compute the Jacobian
      eit->J[z] = dot_product(eit->NY[z], eit->NX) / eit->gX2;

      // Store the smallest and largest Jacobian values
      if(eit->J[z] < xMinJacobian) xMinJacobian = eit->J[z];
      if(eit->J[z] > xMaxJacobian) xMaxJacobian = eit->J[z];

      // Add to the average Jacobian
      xAvgJacobian += eit->J[z];

      // Compute the penalty terms
      // return exp(-a * x) + exp(x - b); 
      eit->PenA[z] = exp(-10 * eit->J[z]);
      eit->PenB[z] = exp(eit->J[z] - 400);
      xTotalPenalty += eit->PenA[z] + eit->PenB[z];
      }
    }

  // Scale the average jacobian
  xAvgJacobian /= S->xAtomGrid->GetNumberOfBoundaryTriangles();

  // Return the total value
  return xTotalPenalty;
}

double BoundaryJacobianEnergyTerm
::ComputePartialDerivative(
  SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // Place to store the Jacobian
  double dTotalPenalty = 0.0;

  // Make sure that the quad entry array has been initialized
  assert(xTriangleEntries.size() == S->xAtomGrid->GetNumberOfTriangles());
  
  // Keep track of the entry index
  TriangleVector::iterator eit = xTriangleEntries.begin();
  
  // Create a Triangle-iterator through the atoms
  
  for(MedialTriangleIterator itt(S->xAtomGrid); !itt.IsAtEnd(); ++itt, ++eit)
    {
    // Get the derivative atoms too
    MedialAtom &dA0 = dS->xAtoms[itt.GetAtomIndex(0)];
    MedialAtom &dA1 = dS->xAtoms[itt.GetAtomIndex(1)];
    MedialAtom &dA2 = dS->xAtoms[itt.GetAtomIndex(2)];

    // Compute the average Xu and Xv vectors and derivatives
    SMLVec3d dXU = dA1.X - dA0.X; SMLVec3d dXV = dA2.X - dA0.X;
    SMLVec3d dNX = vnl_cross_3d(dXU,  eit->XV) + vnl_cross_3d(eit->XU,  dXV);

    // Compute G and its derivative
    double dgX2 = 2.0 * dot_product(dNX, eit->NX);

    // Compute side-wise derivatives
    for(size_t z = 0; z < 2; z++)
      {
      // Compute boundary vector derivatives
      SMLVec3d dYU = dA1.xBnd[z].X - dA0.xBnd[z].X;
      SMLVec3d dYV = dA2.xBnd[z].X - dA0.xBnd[z].X;
      SMLVec3d dNY = vnl_cross_3d(dYU, eit->YV[z]) + vnl_cross_3d(eit->YU[z], dYV);

      // Compute the Jacobian derivative
      double dJ = (
        dot_product(dNY, eit->NX) + 
        dot_product(eit->NY[z], dNX) - eit->J[z] * dgX2) / eit->gX2;

      // Compute the penalty terms
      dTotalPenalty += 
        (-10 * eit->PenA[z] + eit->PenB[z]) * dJ;
      }
    }

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

double ProbabilisticEnergyTerm::ComputeEnergy(SolutionDataBase *S)
{
  // Make sure that the solution data has internal weights and points
  S->UpdateInternalWeights(nCuts);

  // Construct the 'match with the image' function. We simply use the raw
  // image intensity here
  FloatImageEuclideanFunctionAdapter fImage(xImage);

  // Clear the intersection volume accumulator
  xObjectIntegral = 0;

  // double xVolume = 0.0;

  // Iterate over all the points inside the object and take the image value
  for(MedialInternalPointIterator it(S->xAtomGrid, nCuts); !it.IsAtEnd(); ++it)
    {
    size_t i = it.GetIndex();
    xObjectIntegral += 
      S->xInternalWeights[i] * fImage.Evaluate(S->xInternalPoints[i]);
    // xVolume += S->xInternalWeights[i];
    }

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
  for(MedialBoundaryPointIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    size_t i = it.GetIndex();
    BoundaryAtom &bat = GetBoundaryPoint(it, S->xAtoms);
    xImageVal[i] = fImage.Evaluate(bat.X);
    xBndNrm[i] = S->xBoundaryWeights[i] * bat.N;

    // xVolume += dot_product(bat.X, xBndNrm[i]);
    }

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
  for(MedialBoundaryPointIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    size_t i = it.GetIndex();
    double NdX = dot_product(xBndNrm[i], GetBoundaryPoint(it, dS->xAtoms).X);
    dObjectIntegral += xImageVal[i] * NdX;
    }

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
double CrestLaplacianEnergyTerm::ComputeEnergy(SolutionDataBase *S)
{
  // Initialize the accumulators
  xMaxLaplacian = -1e10;
  xTotalPenalty = xAvgLaplacian = 0.0;
  nCrestAtoms = nBadSites = 0;
  
  // Iterate over all crest atoms
  for(MedialAtomIterator it(S->xAtomGrid) ; !it.IsAtEnd(); ++it)
    {
    if(it.IsEdgeAtom())
      {
      // Get the laplacian at this point
      double x = S->xAtoms[it.GetIndex()].xLapR;

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
double AtomBadnessTerm::ComputeEnergy(SolutionDataBase *S)
{
  // Initialize the accumulators
  xMaxBadness = 0;
  xTotalPenalty = xAvgBadness = 0.0;
  nBadAtoms = 0;
  
  // Iterate over all crest atoms
  for(MedialAtomIterator itAtom(S->xAtomGrid) ; !itAtom.IsAtEnd(); ++itAtom)
    {
    // We only care about invalid atoms
    MedialAtom &A = S->xAtoms[itAtom.GetIndex()];
    if(!A.flagValid)
      {
      nBadAtoms++;
      xTotalPenalty += A.xGradRMagSqr - 1.0;
      }
    }

  // Finish up
  nAtoms = S->xAtomGrid->GetNumberOfAtoms();
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
  for(MedialAtomIterator itAtom(S->xAtomGrid) ; !itAtom.IsAtEnd(); ++itAtom)
    {
    size_t i = itAtom.GetIndex();
    if(!S->xAtoms[i].flagValid)
      dTotalPenalty += dS->xAtoms[i].xGradRMagSqr;
    }

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
double CosineSquareTuple(MedialAtom *A, MedialTriangleIterator &it)
{
  // Get the four vectors
  const SMLVec3d &X0 = A[it.GetAtomIndex(0)].X;
  const SMLVec3d &X1 = A[it.GetAtomIndex(1)].X;
  const SMLVec3d &X2 = A[it.GetAtomIndex(2)].X;

  // Compute the differences
  SMLVec3d D10 = X1 - X0;
  SMLVec3d D21 = X2 - X1;
  SMLVec3d D02 = X0 - X2;

  // Compute the lengths squared
  double L10 = dot_product(D10, D10);
  double L21 = dot_product(D21, D21);
  double L02 = dot_product(D02, D02);

  // Compute the cosines squared
  double A0 = dot_product(D10, D02);
  double A1 = dot_product(D21, D10);
  double A2 = dot_product(D02, D21);

  // Compute the actual values
  double C0 = (A0 * A0) / (L10 * L02);
  double C1 = (A1 * A1) / (L21 * L10);
  double C2 = (A2 * A2) / (L02 * L21);

  // Compute the weighted sum
  return C0 + C1 + C2;
}

// Compute squared cosines of four angles in a quad
double CosineSquareTupleDerivative(MedialAtom *A, MedialAtom *dA, MedialTriangleIterator &it)
{
  size_t i0 = it.GetAtomIndex(0);
  size_t i1 = it.GetAtomIndex(1);
  size_t i2 = it.GetAtomIndex(2);

  // Compute the differences and their derivatives
  SMLVec3d D10 = A[i1].X - A[i0].X;
  SMLVec3d D21 = A[i2].X - A[i1].X;
  SMLVec3d D02 = A[i0].X - A[i2].X;

  SMLVec3d dD10 = dA[i1].X - dA[i0].X;
  SMLVec3d dD21 = dA[i2].X - dA[i1].X;
  SMLVec3d dD02 = dA[i0].X - dA[i2].X;

  // Compute the lengths squared
  double L10 = dot_product(D10, D10);
  double L21 = dot_product(D21, D21);
  double L02 = dot_product(D02, D02);

  double dL10 = 2.0 * dot_product(D10, dD10);
  double dL21 = 2.0 * dot_product(D21, dD21);
  double dL02 = 2.0 * dot_product(D02, dD02);

  // Compute the cosines squared
  double A0 = dot_product(D10, D02);
  double A1 = dot_product(D21, D10);
  double A2 = dot_product(D02, D21);

  double dA0 = dot_product(D10, dD02) + dot_product(dD10, D02);
  double dA1 = dot_product(D21, dD10) + dot_product(dD21, D10);
  double dA2 = dot_product(D02, dD21) + dot_product(dD02, D21);

  // Compute the derivatives of the actual values
  double C0 = (A0 * A0) / (L10 * L02);
  double C1 = (A1 * A1) / (L21 * L10);
  double C2 = (A2 * A2) / (L02 * L21);

  // (a/b)' = (a'b - ab')/b^2 = (a' - (a/b)b')/b
  double dC0 = (2.0 * A0 * dA0 - C0 * (L10 * dL02 + dL10 * L02)) / (L10 * L02);
  double dC1 = (2.0 * A1 * dA1 - C1 * (L21 * dL10 + dL21 * L10)) / (L21 * L10);
  double dC2 = (2.0 * A2 * dA2 - C2 * (L02 * dL21 + dL02 * L21)) / (L02 * L21);

  // Compute the weighted sum
  return dC0 + dC1 + dC2;
}

double MedialAnglesPenaltyTerm::ComputeEnergy(SolutionDataBase *S)
{
  xTotalPenalty = 0.0;

  // Sum the penalty over all triangles in the mesh
  for(MedialTriangleIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    xTotalPenalty += CosineSquareTuple(S->xAtoms, it);

  return 0.25 * xTotalPenalty / S->xAtomGrid->GetNumberOfAtoms();
}

double MedialAnglesPenaltyTerm::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  double dTotalPenalty = 0.0;

  // Sum the penalty over all triangles in the mesh
  for(MedialTriangleIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    dTotalPenalty += CosineSquareTupleDerivative(S->xAtoms, dS->xAtoms, it);

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
  MedialIterationContext *inAtomGrid, MedialAtom *inAtomArray)
{
  // Generate an array of weights for irregular grids
  xDomainWeights.set_size(inAtomGrid->GetNumberOfAtoms());
  xDomainArea = ComputeMedialDomainAreaWeights(
    inAtomGrid, inAtomArray, xDomainWeights.data_block());
  cout << "Domain Area: " << xDomainArea << endl;
}

double MedialRegularityTerm::ComputeEnergy(SolutionDataBase *S)
{
  // Integral of (G2(1,11) + G(2,12))^2 + (G(1,21) + G(2,22)) du dv 
  // (no area elt. scaling)
  xGradMagIntegral = 0;
  xGradMagMaximum = 0;
  
  // Integrate over all atoms
  for(MedialAtomIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    // Get the medial atom
    MedialAtom &a = S->xAtoms[it.GetIndex()];

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
    xGradMagIntegral += xDomainWeights[it.GetIndex()] * reg;
    }

  return xGradMagIntegral / xDomainArea;
}

double MedialRegularityTerm
::ComputePartialDerivative(SolutionData *S, PartialDerivativeSolutionData *dS)
{
  double dIntegral = 0.0;
  
  // Integrate over all atoms
  for(MedialAtomIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    // Get the four atoms at the corner of the quad
    MedialAtom &a = S->xAtoms[it.GetIndex()];
    MedialAtom &da = dS->xAtoms[it.GetIndex()];

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
    dIntegral += xDomainWeights[it.GetIndex()] * dreg;
    }

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
  GenericMedialModel *xMedialModel, CoefficientMapping *xCoeff)
{
  this->xCoeff = xCoeff;
  this->xMedialModel = xMedialModel;
  this->nCoeff = xCoeff->GetNumberOfParameters();

  flagLastEvalAvailable = false;
  flagPhiGuessAvailable = false;
  flagGradientComputed = false;

  // Save the initial coefficients currently in the model
  xInitialCoefficients = xMedialModel->GetCoefficientArray();
  
  // Prepare medial atoms for gradient computation (this can only be done for 
  // linear coefficient mappings; for non-linear ones, the variation will change
  // depending on where in parameter space we are).
  for(size_t i = 0; i < nCoeff; i++)
    {
    // Create an array of atoms for this directional derivative
    MedialAtom *a = new MedialAtom[xMedialModel->GetNumberOfAtoms()];

    // Precompute the variation (direction) corresponding to i-th component
    if(xCoeff->IsLinear()) 
      {
      Vec xVariation = xCoeff->GetVariationForParameter(
        xInitialCoefficients, Vec(xCoeff->GetNumberOfParameters(), 0.0), i);

      // Precompute the 'simple' terms in the variational derivative
      xMedialModel->PrepareAtomsForVariationalDerivative(xVariation, a);
      }
    
    // Save the atoms
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

  // Check if we are re-evaluating at a previously tried location
  if(flagLastEvalAvailable && X == xLastEvalPoint)
    return false;

  // Turn of the last eval flag in case that exception is thrown below
  flagLastEvalAvailable = false;

  // Solve the equation
  xSolveTimer.Start();

  // Update the medial model with the new coefficients
  xMedialModel->SetCoefficientArray(xCoeff->Apply(xInitialCoefficients, X));

  // Solve the PDE using the last phi field if we have it
  if(flagGradientComputed)
    xMedialModel->ComputeAtoms(xLastGradHint.data_block());
  else
    xMedialModel->ComputeAtoms();

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
  double t0 = clock();

  bool flagUpdate = SolvePDE(xEvalPoint);

  cout << "[SLV: " << (clock() - t0)/CLOCKS_PER_SEC << " s] " << flush;

  // If the solution changed compute the terms
  if(flagUpdate)
    {
    t0 = clock();

    // Create a solution data object representing current solution
    SolutionData S0(xMedialModel->GetIterationContext(), xMedialModel->GetAtomArray());

    // Compute the solution for each term
    xLastSolutionValue = 0.0;
    xLastTermValues.set_size(xTerms.size());

    for(size_t iTerm = 0; iTerm < xTerms.size(); iTerm++)
      { 
      xTimers[iTerm].Start();
      xLastTermValues[iTerm] = xTerms[iTerm]->ComputeEnergy(&S0);
      xLastSolutionValue += xLastTermValues[iTerm] * xWeights[iTerm]; 
      xTimers[iTerm].Stop();
      }

    cout << "[MAP: " << (clock() - t0)/CLOCKS_PER_SEC << " s] " << endl;
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
  xMedialModel->SetSolutionAsInitialGuess();

  // Create a solution data 
  SolutionData S(xMedialModel);

  // Compute the variation in the gradient direction and the var. deriv.
  IHyperSurface2D *xVariation = xCoeff->GetVariationSurface(xDirection);
  xSolveGradTimer.Start();
  xMedialModel->ComputeVariationalDerivative(xVariation, dAtoms); 
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

void
MedialOptimizationProblem
::ComputeCentralDifferenceGradientPhi(double *x)
{
  // Epsilon used for central differences
  double eps = 0.0001;
  size_t i, n = xMedialModel->GetNumberOfAtoms();

  // Allocate medial atom arrays for central gradient computation
  MedialAtom *a1 = new MedialAtom[n];
  MedialAtom *a2 = new MedialAtom[n];

  // Make a copy of the atoms at the central location
  MedialAtom *a0 = new MedialAtom[n];
  std::copy(xMedialModel->GetAtomArray(), xMedialModel->GetAtomArray()+n, a0);

  // Get a hint array to improve computation of c.d. gradients
  GenericMedialModel::Vec xHint = xMedialModel->GetHintArray(); 

  // Create a timer for this procedure
  CodeTimer tCDG;

  // Compute the central difference gradient for each direction
  for(size_t iCoeff = 0; iCoeff < nCoeff; iCoeff++)
    {
    // Get the current value of the parameters
    double c = xLastEvalPoint[iCoeff];

    // Compute the forward differences
    xLastEvalPoint[iCoeff] = c + eps;
    xMedialModel->SetCoefficientArray(xCoeff->Apply(xInitialCoefficients, xLastEvalPoint));
    xMedialModel->ComputeAtoms(xHint.data_block());

    // Save the atoms
    std::copy(xMedialModel->GetAtomArray(), xMedialModel->GetAtomArray()+n, a1);

    // Compute the backward differences
    xLastEvalPoint[iCoeff] = c - eps;
    xMedialModel->SetCoefficientArray(xCoeff->Apply(xInitialCoefficients, xLastEvalPoint));
    xMedialModel->ComputeAtoms(xHint.data_block());

    // Save the atoms
    std::copy(xMedialModel->GetAtomArray(), xMedialModel->GetAtomArray()+n, a2);

    // Reset the coefficient
    xLastEvalPoint[iCoeff] = c;

    // Compute the atom derivative arrays
    for(i = 0; i< n; i++)
      MedialAtomCentralDifference(a1[i], a2[i], eps, dAtomArray[iCoeff][i]);
    }

  // Print timing information
  cout << "[CDG: " << tCDG.StopAndRead() << " ms] " << flush;

  // Restore the atoms in the solution
  std::copy(a0, a0 + n, xMedialModel->GetAtomArray());

  // Clean up
  delete[] a1; delete[] a2; delete[] a0;
}

double 
MedialOptimizationProblem
::ComputeGradient(double *xEvalPoint, double *xGradient)
{
  size_t iTerm;

  double t00 = clock();

  // Solve the PDE
  double t0 = clock();
  SolvePDE(xEvalPoint);
  cout << " [SLV: " << (clock() - t0) / CLOCKS_PER_SEC << " s] " << flush;

  // Compute the gradient
  xSolveGradTimer.Start();
  
  // If the coefficient mapping is non-linear, we must compute all the variations again
  if(!xCoeff->IsLinear())
    {
    for(size_t i = 0; i < xCoeff->GetNumberOfParameters(); i++) 
      {
      Vec xVariation = xCoeff->GetVariationForParameter(
        xInitialCoefficients, Vec(xCoeff->GetNumberOfParameters(), 0.0), i);

      // Precompute the 'simple' terms in the variational derivative
      xMedialModel->PrepareAtomsForVariationalDerivative(
        xVariation, dAtomArray[i]);
      }
    }

  // Now, compute the actual gradient
  t0 = clock();
  xMedialModel->ComputeAtomGradient(dAtomArray); 
  cout << " [GRAD: " << (clock() - t0) / CLOCKS_PER_SEC << " s] " << flush;
  
  // Stop the timer
  xSolveGradTimer.Stop();

  // Create a solution data object representing current solution
  SolutionData S(xMedialModel->GetIterationContext(), xMedialModel->GetAtomArray());

  // Compute the value at the solution and init gradient computation
  t0 = clock();
  xLastSolutionValue = 0.0;
  xLastTermValues.set_size(xTerms.size());

  for(iTerm = 0; iTerm < xTerms.size(); iTerm++)
    {
    xTimers[iTerm].Start();
    xLastTermValues[iTerm] = xTerms[iTerm]->BeginGradientComputation(&S);
    xLastSolutionValue += xWeights[iTerm] * xLastTermValues[iTerm];
    xTimers[iTerm].Stop();
    }
  cout << " [" << (clock() - t0) / CLOCKS_PER_SEC << " s] " << flush;
    
  // Repeat for each coefficient
  t0 = clock();
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
    }
  cout << " [" << (clock() - t0) / CLOCKS_PER_SEC << " s] " << flush;

  // Clear up gradient computation
  for(iTerm = 0; iTerm < xTerms.size(); iTerm++)
    xTerms[iTerm]->EndGradientComputation();

  // Increment the evaluation counter
  evaluationCost += nCoeff;

  cout << " [[" << (clock() - t00) / CLOCKS_PER_SEC << " s]] " << flush;

  cout << endl;

  // Store the information about the gradient
  xLastGradPoint = vnl_vector<double>(xEvalPoint, nCoeff);
  xLastGradient = vnl_vector<double>(xGradient, nCoeff);
  xLastGradHint = xMedialModel->GetHintArray();
  flagGradientComputed = true;

  // Store the information about the atoms, to guide PDE solving afterwards
  

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
    sout << "  Contribution: " << 
      xWeights[iTerm] << " * " << xLastTermValues[iTerm] <<
      " = " << xWeights[iTerm] * xLastTermValues[iTerm] << endl;
    sout << "  Elapsed time: " << xTimers[iTerm].Read() << endl;
    sout << "  Gradient time: " << xGradTimers[iTerm].Read() << endl;
    }
}

