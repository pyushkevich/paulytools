#ifndef __MedialAtomGrid_h_
#define __MedialAtomGrid_h_

#include "MedialAtom.h"
#include "MedialAtomIterators.h"

using namespace std;

/** 
 * Helper function to access a boundary site in an atom array using a
 * boundary point iterator 
 */
inline BoundaryAtom &
GetBoundaryPoint(MedialBoundaryPointIterator &itBoundary, MedialAtom *xAtoms)
{
  return xAtoms[itBoundary.GetAtomIndex()].xBnd[itBoundary.GetBoundarySide()];
}

/**
 * A helper function that accesses the boundary site for a triangle iterator
 */
inline BoundaryAtom &
GetBoundaryPoint(MedialBoundaryTriangleIterator &itTriangle, MedialAtom *xAtoms, size_t i)
{
  return xAtoms[itTriangle.GetAtomIndex(i)].xBnd[itTriangle.GetBoundarySide()];
}

/** Helper function to compute the coordinate of an internal medial point
 * that is referenced by an iterator */
SMLVec3d GetInternalPoint(MedialInternalPointIterator &it, MedialAtom *xAtoms);

// This method computes the weights for integration over the domain of the medial
// surface. Because the domain may be non-uniform, we must scale all integrals in
// du dv by these weights
double ComputeMedialDomainAreaWeights( MedialIterationContext *xGrid, 
  MedialAtom *xAtoms, double *xWeights);

/**
 * This function computes the area associated with each triangle on the
 * boundary of the object and assigns a third of this area to each of the 
 * vertices of the triangle. Thus for each vertex on the boundary a weight
 * is generated. The total of these weights, equal to the area of the boundary
 * surface is the return value 
 */
double ComputeMedialBoundaryAreaWeights( 
  MedialIterationContext *xGrid, MedialAtom *xAtoms, double *xWeights);

double ComputeMedialInternalVolumePartialDerivativeWeights(
  MedialIterationContext *xGrid, SMLVec3d *xPoints, SMLVec3d *dPoints, size_t nCuts, 
  double *dWeights, double *dInternalProfileWeights);

/**
 * Compute the derivative of the medial boundary weights */
double ComputeMedialBoundaryAreaPartialDerivative(
  MedialIterationContext *xGrid, MedialAtom *xAtoms, MedialAtom *dAtoms, 
  double *xWeights, double *dWeights);

/** 
 * Interpolate a list of internal medial points from the m-rep.
 */
void ComputeMedialInternalPoints(
  MedialIterationContext *xGrid, MedialAtom *xAtoms, size_t nCuts, SMLVec3d *xPoints);

// This is a measure that can be computed over a volume (just a R3 function)
class EuclideanFunction {
public:
  virtual double Evaluate(const SMLVec3d &x) = 0;
  virtual void ComputeGradient(const SMLVec3d &x, SMLVec3d &G)
    { G.fill(0.0); }
};

/**
 * Compute an area of a triangle
 */
inline double TriangleArea(const SMLVec3d &A, const SMLVec3d &B, const SMLVec3d &C)
{
  return 0.5 * vnl_cross_3d(B - A, C - A).magnitude();
}

/** 
 * Compute the partial derivative of the triangle area with respect to 
 * some parameter p.
 */
inline double TriangleAreaPartialDerivative(
  const SMLVec3d &A, const SMLVec3d &B, const SMLVec3d &C, 
  const SMLVec3d &Ap, const SMLVec3d &Bp, const SMLVec3d &Cp,
  double xArea)
{
  SMLVec3d N = 0.5 * vnl_cross_3d(B - A, C - A);
  return - 0.5 * 
    ( dot_product(Ap, vnl_cross_3d(N, (B - C))) + 
      dot_product(Bp, vnl_cross_3d(N, (C - A))) +
      dot_product(Cp, vnl_cross_3d(N, (A - B))) ) / xArea;
}


void TestTriangleAreaPartialDerivative();


/**
 * Compute the volume of a prism formed by two triangles
 */
inline double PrismVolume( SMLVec3d A1, SMLVec3d B1, SMLVec3d C1, 
  SMLVec3d A2, SMLVec3d B2, SMLVec3d C2) 
{
  // TODO: Get/derive the correct formula!!!
  double xArea1 = TriangleArea(A1, B1, C1);
  double xArea2 = TriangleArea(A2, B2, C2);
  double d = ((A1 + B1 + C1) - (A2 + B2 + C2)).magnitude() / 3.0;
  return d * 0.5 * (xArea1 + xArea2);
}

/** Compute the volume of a convex cell. If the cell is not convex, the result
 * will most likely be nonsense */
double CellVolume(
  const SMLVec3d &X000, const SMLVec3d &X001, 
  const SMLVec3d &X010, const SMLVec3d &X011, 
  const SMLVec3d &X100, const SMLVec3d &X101, 
  const SMLVec3d &X110, const SMLVec3d &X111);

/**
 * This function computes the volume associated with each cell in the medial
 * interior and assigns a portion of this cell to each internal point that is
 * the vertex of the cell. The sum of these resulting internal point weights
 * is equal to the volume of the m-rep
 */
double ComputeMedialInternalVolumeWeights(
  MedialIterationContext *xGrid, SMLVec3d *xPoints, size_t nCuts, 
  double *xWeights, double *xProfileIntervalWeights = NULL);

/** This method integrates any given function over the interior of mrep */
double IntegrateFunctionOverInterior (
  MedialIterationContext *xGrid, SMLVec3d *xPoints, 
  double *xWeights, size_t nCuts, EuclideanFunction *fMatch);

/** Integrate any function over the interior of a medial atom */
double IntegrateFunctionOverBoundary (
  MedialIterationContext *xGrid, MedialAtom *xAtoms, 
  double *xWeights, EuclideanFunction *fMatch);

/** Compute gradient of a function over the boundary */
double ComputeFunctionJetOverBoundary(
  MedialIterationContext *xGrid, MedialAtom *xAtoms, 
  double *xWeights, EuclideanFunction *fMatch, SMLVec3d *xOutGradient);


/** 
 * Integrate a function over the boundary, making sure that irregular quads are 
 * super-sampled. This results in better quality integration, but is more costly
 * to compute.
 */
double AdaptivelyIntegrateFunctionOverBoundary(
  MedialIterationContext *xGrid, MedialAtom *xAtoms, 
  double xMinQuadArea, EuclideanFunction *fMatch);


#endif
