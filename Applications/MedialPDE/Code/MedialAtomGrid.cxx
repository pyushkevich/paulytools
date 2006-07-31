#include "MedialAtomGrid.h"


/*****************************************************************************
 * CELL VOLUME AND DERIVATIVE COMPUTATION
 ****************************************************************************/

inline double ScalarTripleProduct(
  const SMLVec3d &A, const SMLVec3d &B, const SMLVec3d &C)
{
  return
    A[0] * (B[1] * C[2] - B[2] * C[1]) + 
    A[1] * (B[2] * C[0] - B[0] * C[2]) + 
    A[2] * (B[0] * C[1] - B[1] * C[0]);
}

inline double ScalarTripleProductDerivative(
  const SMLVec3d &A, const SMLVec3d &B, const SMLVec3d &C,
  const SMLVec3d &DA, const SMLVec3d &DB, const SMLVec3d &DC)
{
  return
    DA[0] * (B[1] * C[2] - B[2] * C[1]) + 
    DA[1] * (B[2] * C[0] - B[0] * C[2]) + 
    DA[2] * (B[0] * C[1] - B[1] * C[0]) +
    A[0] * (DB[1] * C[2] - DB[2] * C[1] + B[1] * DC[2] - B[2] * DC[1]) + 
    A[1] * (DB[2] * C[0] - DB[0] * C[2] + B[2] * DC[0] - B[0] * DC[2]) + 
    A[2] * (DB[0] * C[1] - DB[1] * C[0] + B[0] * DC[1] - B[1] * DC[0]);
}


// This returns 6 x the volume of the tetrahedron
inline double TetrahedronVolume(
  const SMLVec3d &X1, const SMLVec3d &X2, const SMLVec3d &X3, const SMLVec3d X4)
{
  return dot_product(vnl_cross_3d(X2 - X1, X3 - X1), X4 - X1);
}

// This returns 6 x the derivative
inline double TetrahedronVolumePartialDerivative(
  const SMLVec3d &X1, const SMLVec3d &X2, const SMLVec3d &X3, const SMLVec3d X4,
  const SMLVec3d &D1, const SMLVec3d &D2, const SMLVec3d &D3, const SMLVec3d D4)
{
  return 
    dot_product(
      vnl_cross_3d( X1 - X3, D2 ) + 
      vnl_cross_3d( X2 - X1, D3 ) + 
      vnl_cross_3d( X3 - X2, D1 ) , X4 - X1 ) + 
    dot_product( 
      vnl_cross_3d( X2 - X1, X3 - X1 ) , D4 - D1 );
}

/**
 * 
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
  double v00 = TetrahedronVolume(X000, X001, X100, C);
  double v01 = TetrahedronVolume(X101, X100, X001, C);

  double v10 = TetrahedronVolume(X100, X101, X110, C);
  double v11 = - TetrahedronVolume(X111, X110, X101, C);
  
  double v20 = - TetrahedronVolume(X110, X111, X010, C);
  double v21 = TetrahedronVolume(X011, X010, X111, C);
  
  double v30 = TetrahedronVolume(X010, X011, X000, C);
  double v31 = TetrahedronVolume(X001, X000, X011, C);
  
  double v40 = TetrahedronVolume(X011, X111, X001, C);
  double v41 = - TetrahedronVolume(X101, X001, X111, C);
  
  double v50 = - TetrahedronVolume(X010, X000, X110, C);
  double v51 = TetrahedronVolume(X100, X110, X000, C);

   cout << "VOLS: " 
    << (v00 > 0 ? "+1" : "-1") << ", " 
    << (v01 > 0 ? "+1" : "-1") << ", " 
    << (v10 > 0 ? "+1" : "-1") << ", " 
    << (v11 > 0 ? "+1" : "-1") << ", " 
    << (v20 > 0 ? "+1" : "-1") << ", " 
    << (v21 > 0 ? "+1" : "-1") << ", " 
    << (v30 > 0 ? "+1" : "-1") << ", " 
    << (v31 > 0 ? "+1" : "-1") << ", " 
    << (v40 > 0 ? "+1" : "-1") << ", " 
    << (v41 > 0 ? "+1" : "-1") << ", " 
    << (v50 > 0 ? "+1" : "-1") << ", " 
    << (v51 > 0 ? "+1" : "-1") << endl;
  
  // Return a sixth of these sums
  static const double SIXTH = 1.0 / 6.0;
  return - SIXTH *
    (v00 + v01 + v10 + v11 + v20 + v21 + v30 + v31 + v40 + v41 + v50 + v51);
}
*/

/**
 * This in another, less precise, formula for the volume element. We just
 * model the solid as a parallellepiped with sides computed as averages
 */
double CellVolume(
  const SMLVec3d &X000, const SMLVec3d &X001, 
  const SMLVec3d &X010, const SMLVec3d &X011, 
  const SMLVec3d &X100, const SMLVec3d &X101, 
  const SMLVec3d &X110, const SMLVec3d &X111)
{
  SMLVec3d U = X111 + X110 + X101 + X100 - X011 - X010 - X001 - X000;
  SMLVec3d V = X111 + X110 - X101 - X100 + X011 + X010 - X001 - X000;
  SMLVec3d W = X111 - X110 + X101 - X100 + X011 - X010 + X001 - X000;
  return - 0.015625 * ScalarTripleProduct(U, V, W);
} 

double CellVolumePartialDerivative(
  const SMLVec3d &X000, const SMLVec3d &X001, 
  const SMLVec3d &X010, const SMLVec3d &X011, 
  const SMLVec3d &X100, const SMLVec3d &X101, 
  const SMLVec3d &X110, const SMLVec3d &X111,
  const SMLVec3d &D000, const SMLVec3d &D001, 
  const SMLVec3d &D010, const SMLVec3d &D011, 
  const SMLVec3d &D100, const SMLVec3d &D101, 
  const SMLVec3d &D110, const SMLVec3d &D111)
{
  SMLVec3d U  = X111 + X110 + X101 + X100 - X011 - X010 - X001 - X000;
  SMLVec3d V  = X111 + X110 - X101 - X100 + X011 + X010 - X001 - X000;
  SMLVec3d W  = X111 - X110 + X101 - X100 + X011 - X010 + X001 - X000;
  SMLVec3d dU = D111 + D110 + D101 + D100 - D011 - D010 - D001 - D000;
  SMLVec3d dV = D111 + D110 - D101 - D100 + D011 + D010 - D001 - D000;
  SMLVec3d dW = D111 - D110 + D101 - D100 + D011 - D010 + D001 - D000;
  return - 0.015625 * ScalarTripleProductDerivative(U, V, W, dU, dV, dW);
}

/**
 * Here is a similar formula for a wedge volume. The wedge here is treated as
 * a half of a parallellepiped. The vectors forming its edges are computed by 
 * averaging the contributing points
 */
double WedgeVolume(
  const SMLVec3d &A0, const SMLVec3d &B0, const SMLVec3d &C0,
  const SMLVec3d &A1, const SMLVec3d &B1, const SMLVec3d &C1)
{
  const static double FACTOR = 1.0 / 24.0;
  SMLVec3d U = (B1 + B0) - (A1 + A0);
  SMLVec3d V = (C1 + C0) - (A1 + A0);
  SMLVec3d W = (A1 + B1 + C1) - (A0 + B0 + C0);
  
  return FACTOR * ScalarTripleProduct(U, V, W);
}

double WedgeVolumePartialDerivative(
  const SMLVec3d &A0, const SMLVec3d &B0, const SMLVec3d &C0,
  const SMLVec3d &A1, const SMLVec3d &B1, const SMLVec3d &C1,
  const SMLVec3d &dA0, const SMLVec3d &dB0, const SMLVec3d &dC0,
  const SMLVec3d &dA1, const SMLVec3d &dB1, const SMLVec3d &dC1)
{
  const static double FACTOR = 1.0 / 24.0;
  SMLVec3d U = (B1 + B0) - (A1 + A0);
  SMLVec3d V = (C1 + C0) - (A1 + A0);
  SMLVec3d W = (A1 + B1 + C1) - (A0 + B0 + C0);
  SMLVec3d dU = (dB1 + dB0) - (dA1 + dA0);
  SMLVec3d dV = (dC1 + dC0) - (dA1 + dA0);
  SMLVec3d dW = (dA1 + dB1 + dC1) - (dA0 + dB0 + dC0);
  return FACTOR * ScalarTripleProductDerivative(U, V, W, dU, dV, dW);
}

/** Helper function to compute the coordinate of an internal medial point
 * that is referenced by an iterator */
SMLVec3d GetInternalPoint(MedialInternalPointIterator &it, MedialAtom *xAtoms)
{
  // If a medial atom, return the atom
  size_t iAtom = it.GetAtomIndex();
  size_t iDepth = it.GetDepth();

  // If depth is zero (we are on the medial axis) return the medial point
  if(iDepth == 0) return xAtoms[iAtom].X;

  // If depth is equal to max-depth, we can return the boundary site instead
  SMLVec3d Y = xAtoms[iAtom].xBnd[it.GetBoundarySide()].X;
  if(iDepth == it.GetMaxDepth()) return Y;

  // In the remaining case, the point is internal and should be interpolated
  SMLVec3d X = xAtoms[iAtom].X;
  return X + (Y - X) * it.GetRelativeDistanceToMedialAxis();
}


/**
 * This function computes the area associated with each triangle on the
 * boundary of the object and assigns a third of this area to each of the 
 * vertices of the triangle. Thus for each vertex on the boundary a weight
 * is generated. The total of these weights, equal to the area of the boundary
 * surface is the return value 
 */
/* 
double ComputeMedialBoundaryAreaWeights( MedialAtomGrid *xGrid, 
  MedialAtom *xAtoms, double *xWeights)
{
  // Clear the content of the weights
  memset(xWeights, 0, sizeof(double) * xGrid->GetNumberOfBoundaryPoints());
  double xTotalArea = 0.0f;
  
  // Create a quad-based iterator
  MedialBoundaryQuadIterator *itQuad = xGrid->NewBoundaryQuadIterator();

  // For each quad, compute the area associated with it
  for(; !itQuad->IsAtEnd(); ++(*itQuad))
    {
    // Access the four medial points
    const BoundaryAtom& A00 = GetBoundaryPoint(itQuad, xAtoms, 0, 0);
    const BoundaryAtom& A01 = GetBoundaryPoint(itQuad, xAtoms, 0, 1);
    const BoundaryAtom& A11 = GetBoundaryPoint(itQuad, xAtoms, 1, 1);
    const BoundaryAtom& A10 = GetBoundaryPoint(itQuad, xAtoms, 1, 0);

    // Compute the four vectors along the edges of this rectangle
    SMLVec3d U0 = A10.X - A00.X; SMLVec3d V0 = A01.X - A00.X;
    SMLVec3d U1 = A11.X - A01.X; SMLVec3d V1 = A11.X - A10.X;

    // Compute the contribution wo the weight elements of each node
    double W00 = -0.25 * ScalarTripleProduct(U0, V0, A00.N);
    double W10 = -0.25 * ScalarTripleProduct(U0, V1, A10.N);
    double W01 = -0.25 * ScalarTripleProduct(U1, V0, A01.N);
    double W11 = -0.25 * ScalarTripleProduct(U1, V1, A11.N);

    // Update the total area
    xTotalArea += W00 + W10 + W01 + W11;
    
    // Assign a third of each weight to each corner
    xWeights[itQuad->GetBoundaryIndex(0, 0)] += W00;
    xWeights[itQuad->GetBoundaryIndex(1, 1)] += W11;
    xWeights[itQuad->GetBoundaryIndex(0, 1)] += W01;
    xWeights[itQuad->GetBoundaryIndex(1, 0)] += W10;
    }

  // Delete the iterator object
  delete itQuad;

  // Return the total area
  return xTotalArea;
}
*/


// This method computes the weights for integration over the domain of the medial
// surface. Because the domain may be non-uniform, we must scale all integrals in
// du dv by these weights
double ComputeMedialDomainAreaWeights(
  MedialIterationContext *xGrid, MedialAtom *xAtoms, double *xWeights)
{
  // A constant to hold 1/3
  const static double THIRD = 1.0f / 3.0f;
  double xTotalArea = 0.0f;

  // Clear the content of the weights
  memset(xWeights, 0, sizeof(double) * xGrid->GetNumberOfAtoms());
  
  // Iterate over the triangles
  for(MedialTriangleIterator itt(xGrid); !itt.IsAtEnd() ; ++itt)
    {
    // Access the four medial atoms
    size_t i0 = itt.GetAtomIndex(0);
    size_t i1 = itt.GetAtomIndex(1);
    size_t i2 = itt.GetAtomIndex(2);

    // Get the size of the quad. 
    SMLVec3d X0(xAtoms[i0].u, xAtoms[i0].v, 0.0);
    SMLVec3d X1(xAtoms[i1].u, xAtoms[i1].v, 0.0);
    SMLVec3d X2(xAtoms[i2].u, xAtoms[i2].v, 0.0);

    // Compute the area of the triangle, divided by three
    double A = THIRD * TriangleArea( X0, X1, X2 );

    // Add to the total area
    xTotalArea += A;

    // Assign a third of each weight to each corner
    xWeights[i0] += A; xWeights[i1] += A; xWeights[i2] += A;
    }

  // Return the total area
  return 3.0f * xTotalArea;
}


// Old Method: always returns positive areas - dangerous!
double ComputeMedialBoundaryAreaWeights(
  MedialIterationContext *xGrid, MedialAtom *xAtoms, double *xWeights)
{
  // A constant to hold 1/3
  const static double THIRD = 1.0f / 3.0f;
  double xTotalArea = 0.0f;

  // Clear the content of the weights
  memset(xWeights, 0, sizeof(double) * xGrid->GetNumberOfBoundaryPoints());
  
  // Iterate over the triangles
  for(MedialBoundaryTriangleIterator itt(xGrid); !itt.IsAtEnd() ; ++itt)
    {
    // Access the four medial atoms
    size_t i0 = itt.GetBoundaryIndex(0);
    size_t i1 = itt.GetBoundaryIndex(1);
    size_t i2 = itt.GetBoundaryIndex(2);

    // Access the four medial points
    SMLVec3d X0 = GetBoundaryPoint(itt, xAtoms, 0).X;
    SMLVec3d X1 = GetBoundaryPoint(itt, xAtoms, 1).X;
    SMLVec3d X2 = GetBoundaryPoint(itt, xAtoms, 2).X;

    // Compute the area of triangle
    double A = THIRD * TriangleArea( X0, X1, X2 );
    
    // Add to the total area
    xTotalArea += A;

    // Assign a third of each weight to each corner
    xWeights[i0] += A; xWeights[i1] += A; xWeights[i2] += A;
    }

  // Return the total area
  return 3.0f * xTotalArea;
}

double ComputeMedialBoundaryAreaPartialDerivative(
  MedialIterationContext *xGrid, MedialAtom *xAtoms, MedialAtom *dAtoms, 
  double *xWeights, double *dWeights)
{
  // A constant to hold 1/3
  const static double THIRD = 1.0f / 3.0f;
  double dTotalArea = 0.0f;

  // Clear the content of the weights
  memset(dWeights, 0, sizeof(double) * xGrid->GetNumberOfBoundaryPoints());
  
  // For each quad, compute the area associated with it
  for(MedialBoundaryTriangleIterator it(xGrid); !it.IsAtEnd(); ++it)
    {
    // Access the four medial points
    SMLVec3d X0 = GetBoundaryPoint(it, xAtoms, 0).X;
    SMLVec3d X1 = GetBoundaryPoint(it, xAtoms, 1).X;
    SMLVec3d X2 = GetBoundaryPoint(it, xAtoms, 2).X;

    // Access the derivatives of the points
    SMLVec3d D0 = GetBoundaryPoint(it, dAtoms, 0).X;
    SMLVec3d D1 = GetBoundaryPoint(it, dAtoms, 1).X;
    SMLVec3d D2 = GetBoundaryPoint(it, dAtoms, 2).X;

    // Compute the derivative of the triangle area
    double dA = THIRD * 
      TriangleAreaPartialDerivative(X0, X1, X2, D0, D1, D2, TriangleArea(X0, X1, X2));
    
    // Update the total area
    dTotalArea += dA;

    // Assign a third of each weight to each corner
    dWeights[it.GetBoundaryIndex(0)] += dA;
    dWeights[it.GetBoundaryIndex(1)] += dA;
    dWeights[it.GetBoundaryIndex(2)] += dA;
    }

  // Return the total area
  return 3.0f * dTotalArea;
}

/** 
 * Interpolate a list of internal medial points from the m-rep.
 */
void ComputeMedialInternalPoints(
  MedialIterationContext *xGrid, MedialAtom *xAtoms, size_t nCuts, SMLVec3d *xPoints)
{
  // Create an internal point iterator
  for(MedialInternalPointIterator it(xGrid, nCuts); !it.IsAtEnd(); ++it)
    xPoints[it.GetIndex()] = GetInternalPoint(it, xAtoms);
}

class CellVolumeWeightComputer {
public:
  
  CellVolumeWeightComputer(SMLVec3d *xPoints)
    { this->xPoints = xPoints; }
  
  double ComputeWeight (
    size_t i000, size_t i001, size_t i010, size_t i011, 
    size_t i100, size_t i101, size_t i110, size_t i111)
    {
    return CellVolume(
      xPoints[i000], xPoints[i001], xPoints[i010], xPoints[i011], 
      xPoints[i100], xPoints[i101], xPoints[i110], xPoints[i111]); 
    }

private:
  SMLVec3d *xPoints;
};

class CellVolumePartialDerivativeWeightComputer {
public:
  
  CellVolumePartialDerivativeWeightComputer(SMLVec3d *xPoints, SMLVec3d *dPoints)
    { this->xPoints = xPoints; this->dPoints = dPoints; }
  
  double ComputeWeight (
    size_t i000, size_t i001, size_t i010, size_t i011, 
    size_t i100, size_t i101, size_t i110, size_t i111)
    {
    return CellVolumePartialDerivative (
      xPoints[i000], xPoints[i001], xPoints[i010], xPoints[i011], 
      xPoints[i100], xPoints[i101], xPoints[i110], xPoints[i111],
      dPoints[i000], dPoints[i001], dPoints[i010], dPoints[i011], 
      dPoints[i100], dPoints[i101], dPoints[i110], dPoints[i111]); 
    }

private:
  SMLVec3d *xPoints, *dPoints;
};

class WedgeVolumeWeightComputer {
public:
  
  WedgeVolumeWeightComputer(SMLVec3d *xPoints)
    { this->xPoints = xPoints; }
  
  double ComputeWeight(
    size_t i00, size_t i01, size_t i02,
    size_t i10, size_t i11, size_t i12)
    {
    return WedgeVolume(
      xPoints[i00], xPoints[i01], xPoints[i02], 
      xPoints[i10], xPoints[i11], xPoints[i12]);
    }

private:
  SMLVec3d *xPoints;
};

class WedgeVolumePartialDerivativeWeightComputer {
public:
  
  WedgeVolumePartialDerivativeWeightComputer(SMLVec3d *xPoints, SMLVec3d *dPoints)
    { this->xPoints = xPoints; this->dPoints = dPoints; }
  
  double ComputeWeight (
    size_t i00, size_t i01, size_t i02,
    size_t i10, size_t i11, size_t i12)
    {
    return WedgeVolumePartialDerivative (
      xPoints[i00], xPoints[i01], xPoints[i02], 
      xPoints[i10], xPoints[i11], xPoints[i12],
      dPoints[i00], dPoints[i01], dPoints[i02], 
      dPoints[i10], dPoints[i11], dPoints[i12]);
    }

private:
  SMLVec3d *xPoints, *dPoints;
};

/**
 * This generic function reuses code to compute both the cell volume-based
 * weights and their derivatives.
 */
template<class TWeightComputer>
double ComputeMedialInteralWeights(
  MedialIterationContext *xGrid, TWeightComputer *xComputer, size_t nCuts, 
  double *xWeights, double *xInternalProfileWeights)
{
  const static double SIXTH = 1.0 / 6.0;

  // Set all the weights to zero to begin with
  double xTotalVolume = 0.0;
  size_t nPoints = xGrid->GetNumberOfInternalPoints(nCuts);
  memset(xWeights, 0, sizeof(double) * nPoints);
  memset(xInternalProfileWeights, 0, sizeof(double) 
    * xGrid->GetNumberOfProfileIntervals(nCuts));
  
  // Iterate over all the wedges
  for(MedialInternalCellIterator it(xGrid, nCuts); !it.IsAtEnd(); ++it)
    {
    // Get the indices of the six points
    size_t i00 = it.GetInternalPointIndex(0, 0);
    size_t i01 = it.GetInternalPointIndex(1, 0);
    size_t i02 = it.GetInternalPointIndex(2, 0);
    size_t i10 = it.GetInternalPointIndex(0, 1);
    size_t i11 = it.GetInternalPointIndex(1, 1);
    size_t i12 = it.GetInternalPointIndex(2, 1);

    // Add to total volume of the structure
    double V = xComputer->ComputeWeight(i00, i01, i02, i10, i11, i12); 
    xTotalVolume += V;

    // Scale by a sixth (since each wedge has six corners)
    V *= SIXTH;

    // Assign the fractions of the volume to the corresponding internal points
    xWeights[i00] += V; xWeights[i10] += V;
    xWeights[i01] += V; xWeights[i11] += V;
    xWeights[i02] += V; xWeights[i12] += V;

    // Also assign fractions of the weights to the profile intervals
    if(xInternalProfileWeights)
      {
      size_t j0 = it.GetProfileIntervalIndex(0);
      size_t j1 = it.GetProfileIntervalIndex(1);
      size_t j2 = it.GetProfileIntervalIndex(2);
      V = V + V;
      xInternalProfileWeights[j0] += V;
      xInternalProfileWeights[j1] += V;
      xInternalProfileWeights[j2] += V;
      }
    }

  // Return the total area
  cout << "TOTAL VOLUME " << xTotalVolume << endl;
  return xTotalVolume; 
}


/**
 * This function computes the volume associated with each cell in the medial
 * interior and assigns a portion of this cell to each internal point that is
 * the vertex of the cell. The sum of these resulting internal point weights
 * is equal to the volume of the m-rep
 */
double ComputeMedialInternalVolumeWeights(
  MedialIterationContext *xGrid, SMLVec3d *xPoints, size_t nCuts, 
  double *xWeights, double *xInternalProfileWeights)
{
  // Create a specialization object
  WedgeVolumeWeightComputer xComputer(xPoints);

  // Compute the weights
  return ComputeMedialInteralWeights(xGrid, &xComputer, nCuts, 
    xWeights, xInternalProfileWeights);
}

double ComputeMedialInternalVolumePartialDerivativeWeights(
  MedialIterationContext *xGrid, SMLVec3d *xPoints, SMLVec3d *dPoints, size_t nCuts, 
  double *dWeights, double *dInternalProfileWeights)
{
  // Create a specialization object
  WedgeVolumePartialDerivativeWeightComputer xComputer(xPoints, dPoints);

  // Compute the weights
  return ComputeMedialInteralWeights(xGrid, &xComputer, nCuts, 
    dWeights, dInternalProfileWeights);
}

/** This method integrates any given function over the interior of mrep */
double IntegrateFunctionOverInterior (
  MedialIterationContext *xGrid, SMLVec3d *xPoints, 
  double *xWeights, size_t nCuts, EuclideanFunction *fMatch) 
{
  // Match accumulator
  double xMatch = 0.0;
  
  // Create an internal point iterator
  for(MedialInternalPointIterator it(xGrid, nCuts); !it.IsAtEnd(); ++it)
    {
    // Evaluate the image at this location
    size_t iPoint = it.GetIndex();
    xMatch += xWeights[iPoint] * fMatch->Evaluate( xPoints[iPoint] );
    }

  // Return match scaled by total weight
  return xMatch;
}

/** Integrate any function over the interior of a medial atom */
double IntegrateFunctionOverBoundary (
  MedialIterationContext *xGrid, MedialAtom *xAtoms, 
  double *xWeights, EuclideanFunction *fMatch)
{
  // Match accumulator
  double xMatch = 0.0;

  // Iterate through all boundary points
  for(MedialBoundaryPointIterator it(xGrid); !it.IsAtEnd(); ++it)
    {
    double xLocal = fMatch->Evaluate(GetBoundaryPoint(it, xAtoms).X);
    xMatch += xLocal * xWeights[it.GetIndex()];
    }
  
  return xMatch;
}


/********************************************************************************
 * DEPRECATED CODE
 * ******************************************************************************

double RecursiveGetTriangleMatch(
  SMLVec3d &A, SMLVec3d &B, SMLVec3d &C, 
  double xMinArea, EuclideanFunction *fMatch, double &xArea)
{
  // Compute the area of this triangle
  xArea = TriangleArea(A, B, C);

  // If the area is smaller or equal to the min area, return the current match
  if(xArea <= xMinArea)
    return 0.33333333333334 * 
      (fMatch->Evaluate(A) + fMatch->Evaluate(B) + fMatch->Evaluate(C));

  // Otherwise, subdivide the triangle into three points, and recurse for each
  // of the triangles
  double m[4], a[4];
  SMLVec3d AB = 0.5 * (A + B);
  SMLVec3d AC = 0.5 * (A + C);
  SMLVec3d BC = 0.5 * (B + C);

  m[0] = RecursiveGetTriangleMatch(A, AB, AC, xMinArea, fMatch, a[0]);
  m[1] = RecursiveGetTriangleMatch(AB, B, BC, xMinArea, fMatch, a[1]);
  m[2] = RecursiveGetTriangleMatch(AC, BC, C, xMinArea, fMatch, a[2]);
  m[3] = RecursiveGetTriangleMatch(AC, AB, BC, xMinArea, fMatch, a[3]);

  xArea = a[0] + a[1] + a[2] + a[3];
  return (a[0]*m[0] + a[1]*m[1] + a[2]*m[2] + a[3]*m[3]) / xArea;
}


/ ** 
 * Integrate a function over the boundary, making sure that irregular quads are 
 * super-sampled. This results in better quality integration 
 * /
double AdaptivelyIntegrateFunctionOverBoundary(
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, 
  double xMinQuadArea, EuclideanFunction *fMatch)
{
  // Match
  double xMatch = 0.0, xArea = 0.0;

  // Integrate the match over each quad
  MedialBoundaryQuadIterator *itQuad = xGrid->NewBoundaryQuadIterator();
  for(; !itQuad->IsAtEnd(); ++(*itQuad))
    {
    // Get the vector at each vertex
    SMLVec3d x00 = GetBoundaryPoint(itQuad, xAtoms, 0, 0).X;
    SMLVec3d x01 = GetBoundaryPoint(itQuad, xAtoms, 0, 1).X;
    SMLVec3d x11 = GetBoundaryPoint(itQuad, xAtoms, 1, 1).X;
    SMLVec3d x10 = GetBoundaryPoint(itQuad, xAtoms, 1, 0).X;

    // Call the recursive procedure for each sub-quad
    double A1, A2;
    double M1 = RecursiveGetTriangleMatch( x00, x01, x11, xMinQuadArea, fMatch, A1);
    double M2 = RecursiveGetTriangleMatch( x00, x11, x10, xMinQuadArea, fMatch, A2);

    // Add the weighted area to the match
    xMatch += A1 * M1 + A2 * M2;
    xArea += A1 + A2;
    }

  delete itQuad;

  // Scale the match by the area
  return xMatch;
}
*/

/** Compute gradient of a function over the boundary */
double ComputeFunctionJetOverBoundary(
  MedialIterationContext *xGrid, MedialAtom *xAtoms, 
  double *xWeights, EuclideanFunction *fMatch, SMLVec3d *xOutGradient) 
{
  // Compute the image match
  double xMatch = IntegrateFunctionOverBoundary(xGrid, xAtoms, xWeights, fMatch);
  
  // Create a boundary point iterator and iterate over all points
  for(MedialBoundaryPointIterator it(xGrid); !it.IsAtEnd(); ++it)
    {
    // Evaluate the image gradient at this point
    size_t iBnd = it.GetIndex();
    SMLVec3d X = GetBoundaryPoint(it, xAtoms).X;
    fMatch->ComputeGradient(X, xOutGradient[iBnd]);
    xOutGradient[iBnd] *= xWeights[iBnd];
    }

  // Return match scaled by total weight
  return xMatch;
}

void TestTriangleAreaPartialDerivative()
{
  SMLVec3d X[3], Xp[3];
  for(size_t i=0;i<3;i++) for(size_t j=0;j<3;j++)
    {
    X[i][j] = rand() * 1.0 / RAND_MAX;
    Xp[i][j] = rand() * 1.0 / RAND_MAX;
    }

  double eps = 0.001;
  double A0 = TriangleArea(X[0], X[1], X[2]);
  double A1 = 
    TriangleArea(X[0] + eps * Xp[0], X[1] + eps * Xp[1], X[2] + eps * Xp[2]);
  double A2 = 
    TriangleArea(X[0] - eps * Xp[0], X[1] - eps * Xp[1], X[2] - eps * Xp[2]);

  cout << "Numerical derivative : " << 0.5 * (A1 - A2) / eps << endl;
  cout << "Analytical derivative : " << TriangleAreaPartialDerivative(
    X[0], X[1], X[2], Xp[0], Xp[1], Xp[2], A0) << endl;
}
  
void TestTetrahedronVolumePartialDerivative()
{
  double eps = 0.001;
  SMLVec3d X[4], Xp[4], X1[4], X2[4];
  for(size_t i=0;i<4;i++) for(size_t j=0;j<3;j++)
    {
    X[i][j] = rand() * 1.0 / RAND_MAX;
    Xp[i][j] = rand() * 1.0 / RAND_MAX;
    X1[i][j] = X[i][j] + eps * Xp[i][j];
    X2[i][j] = X[i][j] - eps * Xp[i][j];
    }

  double V0 = TetrahedronVolume( X[0], X[1], X[2], X[3] );
  double V1 = TetrahedronVolume( X1[0], X1[1], X1[2], X1[3] );
  double V2 = TetrahedronVolume( X2[0], X2[1], X2[2], X2[3] );

  cout << "Testing tetrahedron partial derivative computation" << endl;
  cout << "Numerical derivative : " << 0.5 * (V1 - V2) / eps << endl;
  cout << "Analytical derivative : " << 
    TetrahedronVolumePartialDerivative(
    X[0], X[1], X[2], X[3], 
    Xp[0], Xp[1], Xp[2], Xp[3]) << endl;
}
  

