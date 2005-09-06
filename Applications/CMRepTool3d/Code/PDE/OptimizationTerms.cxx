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
    { 
    // return FloatImageEuclideanFunctionAdapter::Evaluate(X); 
    return FloatImageEuclideanFunctionAdapter::Evaluate(X);
    }

  SMLVec3d ComputeGradient(const SMLVec3d &X)
    { return FloatImageEuclideanFunctionAdapter::ComputeGradient(X); }
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
    double d = I - xLevelSet;
    return d * d * d * d;
    }

  SMLVec3d ComputeGradient(const SMLVec3d &X)
    {
    // Get the gradient and match from the superclass
    SMLVec3d G = FloatImageEuclideanFunctionAdapter::ComputeGradient(X);
    double I = FloatImageEuclideanFunctionAdapter::Evaluate(X);
    double d = I - xLevelSet;

    // Scale the gradient by the match
    return (4.0 * d * d * d) * G;
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

  // Interpolate the internal points
  ComputeMedialInternalPoints(xAtomGrid, xAtoms, nCuts, xInternalPoints);
  xInternalVolume = ComputeMedialInternalVolumeWeights(
    xAtomGrid, xInternalPoints, nCuts, xInternalWeights, xInternalProfileWeights);

  // Set the flag and store the nCuts
  flagInternalWeights = true;
  nInternalCuts = nCuts;
}

SolutionData::~SolutionData()
{
  if(flagBoundaryWeights)
    { delete xBoundaryWeights; }
  if(flagInternalWeights)
    { delete xInternalWeights; delete xInternalPoints; delete xInternalProfileWeights; }
  if(flagOwnAtoms)
    delete xAtoms;
}

/*********************************************************************************
 * ENERGY TERM
 ********************************************************************************/
double EnergyTerm::BeginGradientComputation(SolutionData *SCenter)
{
  return ComputeEnergy(SCenter);
}

double EnergyTerm::ComputePartialDerivative(
  SolutionData *SCenter, SolutionData *SFwd, SolutionData *SBck,
  double xEpsilon)
{
  // Take the finite difference
  double xDelta = ComputeEnergy(SFwd) - ComputeEnergy(SBck);
  return 0.5 * xDelta / xEpsilon;
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

  // Make sure the weights exist in the solution
  S->UpdateBoundaryWeights();

  // Integrate the image match
  // xImageMatch = IntegrateFunctionOverBoundary(
  //  S->xAtomGrid, S->xAtoms, 
  //  S->xBoundaryWeights, &fImage);

  // Compute the adaptive match, making sure the interpolation works
  double xMinArea = 0.5 * S->xBoundaryArea / S->xAtomGrid->GetNumberOfBoundaryQuads();
  xImageMatch = AdaptivelyIntegrateFunctionOverBoundary(
    S->xAtomGrid, S->xAtoms, xMinArea, &fImage);

  xBoundaryArea = S->xBoundaryArea;
  xFinalMatch = xImageMatch / xBoundaryArea;
  xFinalMatch = sqrt(sqrt(xFinalMatch));

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
/*
double
BoundaryImageMatchTerm::BeginGradientComputation(SolutionData *S)
{
  // Create the image adapter for image / gradient interpolation
  FloatImageLevelSetFunctionAdapter fImage(xImage);

  // Boundary weights will be needed
  S->UpdateBoundaryWeights();

  // Compute the image and image gradient at each point in the image
  xGradI = new SMLVec3d[S->xAtomGrid->GetNumberOfBoundaryPoints()];
  xImageVal = new double[S->xAtomGrid->GetNumberOfBoundaryPoints()];
  xMatch = 0.0;

  MedialBoundaryPointIterator *itBnd = S->xAtomGrid->NewBoundaryPointIterator();
  while(!itBnd->IsAtEnd())
    {
    // Compute the image gradient
    SMLVec3d &X = GetBoundaryPoint(itBnd, S->xAtoms).X;
    xGradI[itBnd->GetIndex()] = fImage.ComputeGradient(X);

    // Compute the image value
    xImageVal[itBnd->GetIndex()] = fImage.Evaluate(X);

    // Accumulate to get weighted match
    xMatch += xImageVal[itBnd->GetIndex()] * S->xBoundaryWeights[itBnd->GetIndex()];
    
    ++(*itBnd);
    }
  
  // We will need the area in many calculations
  xArea = S->xBoundaryArea;
  
  // Clean up
  delete itBnd;

  // Return the solution
  return xMatch / xArea;
}

double
BoundaryImageMatchTerm::ComputePartialDerivative(
    SolutionData *SCenter, SolutionData *SFwd, SolutionData *SBck,
    double xEpsilon)
{
  // Create the image adapter for image / gradient interpolation
  FloatImageLevelSetFunctionAdapter fImage(xImage);

  // Accumulator for the partial derivative of the weighted match function
  double dMatchdC = 0.0;

  // We need the weights for the derivative terms
  SFwd->UpdateBoundaryWeights();
  SBck->UpdateBoundaryWeights();

  // Compute the partial derivative for this coefficient
  MedialBoundaryPointIterator *itBnd = 
    SCenter->xAtomGrid->NewBoundaryPointIterator();
  
  for( ; !itBnd->IsAtEnd(); ++(*itBnd))
    {
    // Get the index of this boundary point
    size_t iPoint = itBnd->GetIndex();

    // Get the boundary point before and after small deformation
    SMLVec3d X2 = GetBoundaryPoint(itBnd, SFwd->xAtoms).X;
    SMLVec3d X1 = GetBoundaryPoint(itBnd, SBck->xAtoms).X;
    SMLVec3d dX = 0.5 * (X2 - X1);

    // Get the area weights for this point
    double w0 = SCenter->xBoundaryWeights[ iPoint ];
    double w1 = SBck->xBoundaryWeights[ iPoint ];
    double w2 = SFwd->xBoundaryWeights[ iPoint ];
    double dw = 0.5 * (w2 - w1);

    // Compute the change in intensity per change in coefficient
    double dIdC = dot_product( xGradI[iPoint] , dX);
    dMatchdC += dIdC * w0 + xImageVal[iPoint] * dw;
    }
  
  delete itBnd;

  // Scale by epsilon to get the derivative
  dMatchdC /= xEpsilon;

  // Compute the derivative of M / A
  double dAdC = 0.5 * (SFwd->xBoundaryArea - SBck->xBoundaryArea) / xEpsilon;
  return (dMatchdC * xArea - dAdC * xMatch) / (xArea * xArea);
}

void BoundaryImageMatchTerm::EndGradientComputation()
{
  // Clean up
  delete xGradI;
  delete xImageVal;
}
*/
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
  // Clear the intersection volume accumulator
  xIntersect = 0;
  
  // Make sure that the solution data has internal weights and points
  data->UpdateInternalWeights(nCuts);

  // Construct the 'match with the image' function
  FloatImageVolumeElementFunctionAdapter fImage(xImage);

  // Interpolate the image at each internal point
  size_t nInternal = data->xAtomGrid->GetNumberOfInternalPoints(nCuts);
  double *xImageValue = new double[nInternal];
  for(size_t i = 0; i < nInternal; i++)
    xImageValue[i] = fImage.Evaluate(data->xInternalPoints[i]);

  // Iterate over all the profiles (intervals of the sail vectors)
  /* 
  MedialProfileIntervalIterator *itProfile = 
    data->xAtomGrid->NewProfileIntervalIterator(nCuts);
  
  for( ; !itProfile->IsAtEnd(); ++(*itProfile) )
    {
    // Get the two points on this profile stick
    size_t i1 = itProfile->GetInnerPointIndex();
    size_t i2 = itProfile->GetOuterPointIndex();
    double m1 = xImageValue[i1], m2 = xImageValue[i2];
    
    // Compute the fraction of the profile that lies inside the positive
    // range of the image. 
    double xLocal;    
    if(m1 > 0 && m2 > 0)
      xLocal = 1.0;
    else if(m1 <= 0 && m2 <= 0)
      xLocal = 0.0;
    else if(m1 > 0)
      xLocal = m1 / (m1 - m2);
    else
      xLocal = m2 / (m2 - m1);

    // Get the weight assigned to this profile
    double w = data->xInternalProfileWeights[itProfile->GetIndex()];

    // Add the weighted sum to the match
    xIntersect += w * xLocal;
    }
  
  delete itProfile; 
  */
  
  MedialInternalPointIterator *it = data->xAtomGrid->NewInternalPointIterator(nCuts);
  for(; !it->IsAtEnd(); ++(*it))
    {
    double xLocal = 0.5 * (1.0 + xImageValue[it->GetIndex()]);
    double w = data->xInternalWeights[it->GetIndex()];
    xIntersect += xLocal * w;
    }
  delete it;
  
  // Get (remember) the object volume 
  xObjectVolume = data->xInternalVolume;

  // Compute the volume overlap measure
  xOverlap = xIntersect / (xImageVolume + xObjectVolume - xIntersect); 

  // Clean up
  delete xImageValue;

  // Return a value that can be minimized
  return 1.0 - xOverlap;
}

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
 * Atom Badness Penalty term
 ********************************************************************************/
double AtomBadnessTerm::ComputeEnergy(SolutionData *data)
{
  // Initialize the accumulators
  xMaxBadness = 0;
  xTotalPenalty = xAvgBadness = 0.0;
  nBadAtoms = nAtoms = 0;
  
  // Iterate over all crest atoms
  MedialAtomIterator *itAtom = data->xAtomGrid->NewAtomIterator();
  for( ; !itAtom->IsAtEnd(); ++(*itAtom) )
    {
    double x = data->xAtoms[itAtom->GetIndex()].xBadness;

    // Update the minimum value
    if(x > xMaxBadness) xMaxBadness = x;

    // Update the bad atom counter 
    if(x > 0) nBadAtoms++;
    nAtoms++;

    // Update the average
    xAvgBadness += x;

    // Update the total penalty
    xTotalPenalty += 100 * x;
    }
  delete itAtom;

  // Finish up
  xAvgBadness /= nBadAtoms;

  return xTotalPenalty;
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
    xEdgeLength.push_back( vnl_vector_ssd(a00.X, a01.X) );
    xEdgeLength.push_back( vnl_vector_ssd(a01.X, a11.X) );
    xEdgeLength.push_back( vnl_vector_ssd(a11.X, a10.X) );
    xEdgeLength.push_back( vnl_vector_ssd(a10.X, a00.X) );
    xEdgeLength.push_back( vnl_vector_ssd(a00.X, a11.X) );
    xEdgeLength.push_back( vnl_vector_ssd(a01.X, a10.X) );
    }

  delete itQuad;
}

double MedialRegularityTerm::ComputeDistortionPenalty(
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

double MedialRegularityTerm::ComputeEnergy(SolutionData *data)
{
  // The prior is based on the distortion in terms of distances between adjacent 
  // atoms as compared to the template. We iterate over the edges on the medial 
  // surface and compute these distortions
  xTotalPenalty = 0.0;
  xMaxDistortion = xMinDistortion = 0.0;
  
  // Iterate over the quads
  MedialQuadIterator *itQuad = data->xAtomGrid->NewQuadIterator();
  vector<double>::iterator itDist = xEdgeLength.begin();
  for(; !itQuad->IsAtEnd(); ++(*itQuad), ++itDist)
    {
    // Get the four atoms at the corner of the quad
    MedialAtom &a00 = data->xAtoms[itQuad->GetAtomIndex(0, 0)];
    MedialAtom &a01 = data->xAtoms[itQuad->GetAtomIndex(0, 1)];
    MedialAtom &a10 = data->xAtoms[itQuad->GetAtomIndex(1, 0)];
    MedialAtom &a11 = data->xAtoms[itQuad->GetAtomIndex(1, 1)];

    // Record the distances (repeat in this order later)
    xTotalPenalty += ComputeDistortionPenalty( *itDist, vnl_vector_ssd(a00.X, a01.X) );
    xTotalPenalty += ComputeDistortionPenalty( *itDist, vnl_vector_ssd(a01.X, a11.X) );
    xTotalPenalty += ComputeDistortionPenalty( *itDist, vnl_vector_ssd(a11.X, a10.X) );
    xTotalPenalty += ComputeDistortionPenalty( *itDist, vnl_vector_ssd(a10.X, a00.X) );
    xTotalPenalty += ComputeDistortionPenalty( *itDist, vnl_vector_ssd(a00.X, a11.X) );
    xTotalPenalty += ComputeDistortionPenalty( *itDist, vnl_vector_ssd(a01.X, a10.X) );
    }

  // Compute the average error
  xTotalPenalty /= xEdgeLength.size();
  xMeanSquareDistortion = xTotalPenalty;

  // Return the error
  return xTotalPenalty; 
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
  SolutionData S0(xSolver, false);

  // cout << "REPORT:" << endl;

  // Compute the solution for each term
  xSolution = 0.0;
  for(size_t iTerm = 0; iTerm < xTerms.size(); iTerm++)
    { xSolution += xTerms[iTerm]->ComputeEnergy(&S0) * xWeights[iTerm]; }

  cout << "RESULT:" << xSolution << endl;
  PrintReport(cout);

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

  // Use current solution as the guess
  xSolver->SetSolutionAsInitialGuess();

  // Create a solution data object representing current solution
  SolutionData *SCenter = new SolutionData(xSolver, true);
  
  // Compute the value at the solution and init gradient computation
  xSolution = 0.0;
  for(iTerm = 0; iTerm < xTerms.size(); iTerm++)
    xSolution += 
      xWeights[iTerm] * xTerms[iTerm]->BeginGradientComputation(SCenter);

  // Repeat for each coefficient
  for(size_t iCoeff = 0; iCoeff < nCoeff; iCoeff++)
    {
    // Increment the coefficient by epsilon
    double xOldValue = xCoeff->GetCoefficient(iCoeff);

    // Compute the forward difference
    xCoeff->SetCoefficient(iCoeff, xOldValue + xEpsilon);
    xSolver->Solve(xPrecision); 
    SolutionData *SFwd = new SolutionData(xSolver, true);

    // Compute the backward difference
    xCoeff->SetCoefficient(iCoeff, xOldValue - xEpsilon);
    xSolver->Solve(xPrecision); 
    SolutionData *SBck = new SolutionData(xSolver, true);

    // Restore the old solution value
    xCoeff->SetCoefficient(iCoeff, xOldValue);

    // Compute the partial derivatives for each term
    XGradient[iCoeff] = 0.0;
    for(iTerm = 0; iTerm < xTerms.size(); iTerm++)
      XGradient[iCoeff] += xWeights[iTerm] *
        xTerms[iTerm]->ComputePartialDerivative(SCenter, SFwd, SBck, xEpsilon);

    // Dump the gradient
    // cout << iCoeff << "; " << XGradient[iCoeff] << endl;
    cout << "." << flush;

    // Delete the solutions
    delete SFwd; delete SBck;
    }

  // Increment the evaluation counter
  cout << endl;
  evaluationCost += nCoeff;

  // Clean up everything
  delete SCenter;

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

