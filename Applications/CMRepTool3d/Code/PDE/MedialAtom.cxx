#include "MedialAtom.h"

void
MedialAtom
::ComputeDifferentialGeometry()
{
  G.SetJet( X.data_block(), Xu.data_block(), Xv.data_block(),
    Xuu.data_block(), Xuv.data_block(), Xvv.data_block());
}

void MedialAtom ::ComputeNormalVector()
{
  N = cross_3d(Xu, Xv) * sqrt(G.gInv);
}

bool MedialAtom::ComputeBoundaryAtoms()
{
  // Terms going into gradR
  double CXu = 
    Ru * G.xContravariantTensor[0][0] + Rv * G.xContravariantTensor[0][1]; 
  double CXv = 
    Ru * G.xContravariantTensor[0][1] + Rv * G.xContravariantTensor[1][1];
    
  // Compute Grad R
  xGradR[0] = Xv[0] * CXv + Xu[0] * CXu;
  xGradR[1] = Xv[1] * CXv + Xu[1] * CXu;
  xGradR[2] = Xv[2] * CXv + Xu[2] * CXu;

  // Compute squared length of GradR
  double xMagGradR2 = Ru * CXu + Rv * CXv;
  double sinTheta2 = 1.0f - xMagGradR2;
  
  // Correct the floating point / solver error
  if(fabs(sinTheta2) < 1e-8) sinTheta2 = 0.0;
  
  // Compute the boundary sites
  if (sinTheta2 > 0.0)
    {
    // We are not at a crest, valid atom
    flagCrest = false; flagValid = true;
    
    // Compute the normal and tangent components of the boundaries
    SMLVec3d CN = N * sqrt(sinTheta2);

    // Compute the position and normals of boundary points
    xBnd[0].N = - xGradR - CN;
    xBnd[1].N = - xGradR + CN;
    xBnd[0].X = X + xBnd[0].N * R;
    xBnd[1].X = X + xBnd[1].N * R;
    }
  else if (sinTheta2 < 0.0)
    { 
    flagValid = false; 
    // cerr << "Bad atom ("<< R << ", " << sinTheta2 << ") at " << u << ", " << v << endl;
    cout << "x" << flush;
    }
  else
    { 
    // We are at a crest, valid atom
    flagValid = true; flagCrest = true;

    // Simpler geometry, save a square root!
    xBnd[0].N = xBnd[1].N = -xGradR;
    xBnd[0].X = xBnd[1].X = X - xGradR * R;
    }

  return flagValid;
}

double CellVolume(
  const SMLVec3d &X000, const SMLVec3d &X001, 
  const SMLVec3d &X010, const SMLVec3d &X011, 
  const SMLVec3d &X100, const SMLVec3d &X101, 
  const SMLVec3d &X110, const SMLVec3d &X111)
{
  // Compute the centerpoint inside the cell
  SMLVec3d C = (X000 + X001 + X010 + X011 + X100 + X101 + X110 + X111) * 0.125;

  // Compute each of the twelve tetrahedra inside the cuboid (group by the
  // face of the cell)
  double v00 = dot_product(cross_3d(X001 - X000, X100 - X000), C - X000);
  double v01 = dot_product(cross_3d(X100 - X101, X001 - X101), C - X101);

  double v10 = dot_product(cross_3d(X101 - X100, X110 - X100), C - X100);
  double v11 = dot_product(cross_3d(X110 - X111, X101 - X111), C - X111);
  
  double v20 = dot_product(cross_3d(X111 - X110, X010 - X110), C - X110);
  double v21 = dot_product(cross_3d(X010 - X011, X111 - X011), C - X011);
  
  double v30 = dot_product(cross_3d(X011 - X010, X000 - X010), C - X010);
  double v31 = dot_product(cross_3d(X000 - X001, X011 - X001), C - X001);
  
  double v40 = dot_product(cross_3d(X111 - X011, X001 - X011), C - X011);
  double v41 = dot_product(cross_3d(X001 - X101, X111 - X101), C - X101);
  
  double v50 = dot_product(cross_3d(X000 - X010, X110 - X010), C - X010);
  double v51 = dot_product(cross_3d(X110 - X100, X000 - X100), C - X100);

  // Return a sixth of these sums
  static const double SIXTH = 1.0 / 6.0;
  return SIXTH *
    (v00 + v01 + v10 + v11 + v20 + v21 + v30 + v31 + v40 + v41 + v50 + v51);
}

/** Helper function to compute the coordinate of an internal medial point
 * that is referenced by an iterator */
SMLVec3d GetInternalPoint(MedialInternalPointIterator *it, MedialAtom *xAtoms)
{
  // If a medial atom, return the atom
  size_t iAtom = it->GetAtomIndex();
  size_t iDepth = it->GetDepth();

  // If depth is zero (we are on the medial axis) return the medial point
  if(iDepth == 0) return xAtoms[iAtom].X;

  // If depth is equal to max-depth, we can return the boundary site instead
  SMLVec3d Y = xAtoms[iAtom].xBnd[it->GetBoundarySide()].X;
  if(iDepth == it->GetMaxDepth()) return Y;

  // In the remaining case, the point is internal and should be interpolated
  SMLVec3d X = xAtoms[iAtom].X;
  return X + (Y - X) * it->GetRelativeDistanceToMedialAxis();
}

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

/** 
 * Interpolate a list of internal medial points from the m-rep.
 */
void ComputeMedialInternalPoints(
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, size_t nCuts, SMLVec3d *xPoints)
{
  // Create an internal point iterator
  MedialInternalPointIterator *it = xGrid->NewInternalPointIterator(nCuts);
  while(!it->IsAtEnd())
    {
    // Compute the internal point
    xPoints[it->GetIndex()] = GetInternalPoint(it, xAtoms);
    ++(*it);
    }
  delete it;
}

/**
 * This function computes the volume associated with each cell in the medial
 * interior and assigns a portion of this cell to each internal point that is
 * the vertex of the cell. The sum of these resulting internal point weights
 * is equal to the volume of the m-rep
 */
double ComputeMedialInternalVolumeWeights(
  MedialAtomGrid *xGrid, SMLVec3d *xPoints, size_t nCuts, 
  double *xWeights, double *xInternalProfileWeights)
{
  // Set all the weights to zero to begin with
  size_t nPoints = xGrid->GetNumberOfInternalPoints(nCuts);
  memset(xWeights, 0, sizeof(double) * nPoints);
  memset(xInternalProfileWeights, 0, sizeof(double) 
    * xGrid->GetNumberOfProfileIntervals(nCuts));
  
  // Create an internal cell iterator
  MedialInternalCellIterator *itCell = xGrid->NewInternalCellIterator(nCuts);

  // For each cell, compute the internal volume
  double xTotalVolume = 0.0;
  while(!itCell->IsAtEnd())
    {
    // Get the indices of the cell corners
    size_t i000 = itCell->GetInternalPointIndex(0, 0, 0);
    size_t i001 = itCell->GetInternalPointIndex(0, 0, 1);
    size_t i010 = itCell->GetInternalPointIndex(0, 1, 0);
    size_t i100 = itCell->GetInternalPointIndex(1, 0, 0);
    size_t i011 = itCell->GetInternalPointIndex(0, 1, 1);
    size_t i101 = itCell->GetInternalPointIndex(1, 0, 1);
    size_t i110 = itCell->GetInternalPointIndex(1, 1, 0);
    size_t i111 = itCell->GetInternalPointIndex(1, 1, 1);

    // Compute the cell's volume
    // cout.width(16);
    // cout << "Cell : " << itCell->GetAtomIndex(0,0) << " ";
    double V = -CellVolume(
      xPoints[i000], xPoints[i001], xPoints[i010], xPoints[i011], 
      xPoints[i100], xPoints[i101], xPoints[i110], xPoints[i111]); 

    // Add to total volume of the structure
    xTotalVolume += V;

    // Scale by an eigth (since each cube has eight corners)
    V *= 0.125;

    // Assign the fractions of the volume to the corresponding internal points
    xWeights[i000] += V; xWeights[i110] += V;
    xWeights[i010] += V; xWeights[i100] += V;
    xWeights[i001] += V; xWeights[i111] += V;
    xWeights[i011] += V; xWeights[i101] += V;

    // Also assign fractions of the weights to the profile intervals
    if(xInternalProfileWeights)
      {
      size_t j00 = itCell->GetProfileIntervalIndex(0, 0);
      size_t j01 = itCell->GetProfileIntervalIndex(0, 1);
      size_t j10 = itCell->GetProfileIntervalIndex(1, 0);
      size_t j11 = itCell->GetProfileIntervalIndex(1, 1);
      V *= 2.0;
      xInternalProfileWeights[j00] += V;
      xInternalProfileWeights[j01] += V;
      xInternalProfileWeights[j10] += V;
      xInternalProfileWeights[j11] += V;
      }

    // Go to the next boundary quad
    ++(*itCell);
    }

  // Delete the iterator object
  delete itCell;

  // Return the total area
  return xTotalVolume;
}

/** This method integrates any given function over the interior of mrep */
double IntegrateFunctionOverInterior (
  MedialAtomGrid *xGrid, SMLVec3d *xPoints, 
  double *xWeights, size_t nCuts, EuclideanFunction *fMatch) 
{
  // Match accumulator
  double xMatch = 0.0;
  
  // Create an internal point iterator
  MedialInternalPointIterator *it = xGrid->NewInternalPointIterator(nCuts);
  while(!it->IsAtEnd())
    {
    // Evaluate the image at this location
    size_t iPoint = it->GetIndex();
    xMatch += xWeights[iPoint] * fMatch->Evaluate( xPoints[iPoint] );

    // On to the next point
    ++(*it); 
    }

  // Clean up
  delete it;

  // Return match scaled by total weight
  return xMatch;
}

/** Integrate any function over the interior of a medial atom */
double IntegrateFunctionOverBoundary (
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, 
  double *xWeights, EuclideanFunction *fMatch)
{
  // Match accumulator
  double xMatch = 0.0;

  // Iterate through all boundary points
  MedialBoundaryPointIterator *itBoundary = xGrid->NewBoundaryPointIterator();
  while(!itBoundary->IsAtEnd())
    {
    double xLocal = fMatch->Evaluate(GetBoundaryPoint(itBoundary, xAtoms).X);
    xMatch += xLocal * xWeights[itBoundary->GetIndex()];
    ++(*itBoundary);
    }
  
  delete itBoundary;
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




