#ifndef __MedialAtomGrid_h_
#define __MedialAtomGrid_h_

#include <vector>
#include "MedialAtom.h"
using namespace std;

/**
 * An abstract iterator for parsing points in the medial atom grid
 */
class MedialAtomIterator
{
public:
  virtual size_t GetIndex() = 0;
  virtual bool IsEdgeAtom() = 0;

  virtual MedialAtomIterator &operator++() = 0;
  virtual bool IsAtEnd() = 0;
  virtual void GoToBegin() = 0;
};

/**
 * An iterator for parsing quad cells in a medial atom grid
 */
class MedialQuadIterator
{
public:
  // Get index for one of the four points associated with this quad
  // (i and j are between 0 and 1 both
  virtual size_t GetAtomIndex(size_t i, size_t j) = 0;

  virtual MedialQuadIterator &operator++() = 0;
  virtual bool IsAtEnd() = 0;
  virtual void GoToBegin() = 0;
};

/**
 * An iterator for parsing boundary points associated with a medial grid
 */
class MedialBoundaryPointIterator
{
public:
  virtual size_t GetIndex() = 0;
  virtual size_t GetOppositeIndex() = 0;
  virtual size_t GetAtomIndex() = 0;
  virtual bool IsEdgeAtom() = 0;
  virtual size_t GetBoundarySide() = 0;

  virtual MedialBoundaryPointIterator &operator++() = 0;
  virtual bool IsAtEnd() = 0;
  virtual void GoToBegin() = 0;
};

/**
 * An iterator for parsing quad cells in a medial atom grid
 */
class MedialBoundaryQuadIterator
{
public:
  // Get index for one of the four points associated with this quad
  // (i and j are between 0 and 1 both
  virtual size_t GetBoundaryIndex(size_t i, size_t j) = 0;
  virtual size_t GetAtomIndex(size_t i, size_t j) = 0;
  virtual size_t GetBoundarySide() = 0;

  virtual MedialBoundaryQuadIterator &operator++() = 0;
  virtual bool IsAtEnd() = 0;
  virtual void GoToBegin() = 0;
};

/** An iterator that goes through the internal points in an m-rep */
class MedialInternalPointIterator
{
public:
  virtual size_t GetIndex() = 0;
  virtual size_t GetAtomIndex() = 0;
  virtual size_t GetMaxDepth() = 0;
  virtual size_t GetDepth() = 0;

  // Are we on a crest (edge of medial atom)?
  virtual bool IsEdgeAtom() = 0;

  // Get a non-negative value t that is zero on medial axis, 1 on boundary
  virtual double GetRelativeDistanceToMedialAxis() = 0;

  // This is only defined for points that are not on the medial axis and 
  // not on the crest (edge of the medial surface)
  virtual size_t GetBoundarySide() = 0;

  virtual MedialBoundaryQuadIterator &operator++() = 0;
  virtual bool IsAtEnd() = 0;
  virtual void GoToBegin() = 0;
};

/** An iterator that goes through the internal cells in an m-rep */
class MedialInternalCellIterator
{
public:
  virtual size_t GetAtomIndex(size_t i, size_t j) = 0;
  virtual size_t GetInternalPointIndex(size_t i, size_t j, size_t k) = 0;
  virtual size_t GetProfileIntervalIndex(size_t i, size_t j) = 0;

  /** Get the depth of the side k (0 inner, 1 outer) */
  virtual size_t GetDepth(size_t k) = 0;
  
  virtual MedialInternalCellIterator &operator++() = 0;
  virtual bool IsAtEnd() = 0;
  virtual void GoToBegin() = 0;
};

/** An iterator that goes through point pairs on the profiles that connect the
 * medial axis to the boundary. There are four such intervals for every
 * internal cell */
class MedialProfileIntervalIterator
{
public:
  virtual size_t GetInnerPointIndex() = 0;
  virtual size_t GetOuterPointIndex() = 0;
  virtual size_t GetIndex() = 0;

  virtual MedialProfileIntervalIterator &operator++() = 0;
  virtual bool IsAtEnd() = 0;
  virtual void GoToBegin() = 0;
};

/**
 * A generic grid of medial atoms with iterators. This class should be
 * able to represent cartesian and polar grids transparently to the
 * users of the class. The class does not store any medial atoms or associated
 * data; it's just used to index the data properly
 */
class MedialAtomGrid
{
public:
  // Create a point iterator 
  virtual MedialAtomIterator *NewAtomIterator() = 0;
  virtual MedialQuadIterator *NewQuadIterator() = 0;
  virtual MedialBoundaryPointIterator *NewBoundaryPointIterator() = 0;
  virtual MedialBoundaryQuadIterator *NewBoundaryQuadIterator() = 0;
  virtual MedialInternalCellIterator *NewInternalCellIterator(size_t nCuts) = 0;
  virtual MedialInternalPointIterator *NewInternalPointIterator(size_t nCuts) = 0;
  virtual MedialProfileIntervalIterator *NewProfileIntervalIterator(size_t nCuts) = 0;

  // Get the number of medial atoms in this grid
  virtual size_t GetNumberOfAtoms() = 0;
  virtual size_t GetNumberOfQuads() = 0;
  virtual size_t GetNumberOfBoundaryPoints() = 0;
  virtual size_t GetNumberOfBoundaryQuads() = 0;
  virtual size_t GetNumberOfCells(size_t nCuts) = 0;
  virtual size_t GetNumberOfInternalPoints(size_t nCuts) = 0;
  virtual size_t GetNumberOfProfileIntervals(size_t nCuts) = 0;
};

/** Helper function to access a boundary site in an atom array using a
 * boundary point iterator */
inline BoundaryAtom &
GetBoundaryPoint(MedialBoundaryPointIterator *itBoundary, MedialAtom *xAtoms)
{
  return 
    xAtoms[itBoundary->GetAtomIndex()].xBnd[itBoundary->GetBoundarySide()];
}

/** Helper function to access a boundary site in an atom array using a
 * boundary point iterator */
inline BoundaryAtom &
GetBoundaryPoint(MedialBoundaryQuadIterator *itBQuad, MedialAtom *xAtoms, 
  unsigned int i, unsigned int j)
{
  return 
    xAtoms[itBQuad->GetAtomIndex(i,j)].xBnd[itBQuad->GetBoundarySide()];
}

/** Helper function to compute the coordinate of an internal medial point
 * that is referenced by an iterator */
SMLVec3d GetInternalPoint(MedialInternalPointIterator *it, MedialAtom *xAtoms);


// This method computes the weights for integration over the domain of the medial
// surface. Because the domain may be non-uniform, we must scale all integrals in
// du dv by these weights
double ComputeMedialDomainAreaWeights( MedialAtomGrid *xGrid, 
  MedialAtom *xAtoms, double *xWeights);

/**
 * This function computes the area associated with each triangle on the
 * boundary of the object and assigns a third of this area to each of the 
 * vertices of the triangle. Thus for each vertex on the boundary a weight
 * is generated. The total of these weights, equal to the area of the boundary
 * surface is the return value 
 */
double ComputeMedialBoundaryAreaWeights( 
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, double *xWeights);

double ComputeMedialInternalVolumePartialDerivativeWeights(
  MedialAtomGrid *xGrid, SMLVec3d *xPoints, SMLVec3d *dPoints, size_t nCuts, 
  double *dWeights, double *dInternalProfileWeights);

/**
 * Compute the derivative of the medial boundary weights */
double ComputeMedialBoundaryAreaPartialDerivative(
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, MedialAtom *dAtoms, 
  double *xWeights, double *dWeights);

/** 
 * Interpolate a list of internal medial points from the m-rep.
 */
void ComputeMedialInternalPoints(
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, size_t nCuts, SMLVec3d *xPoints);

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
  MedialAtomGrid *xGrid, SMLVec3d *xPoints, size_t nCuts, 
  double *xWeights, double *xProfileIntervalWeights = NULL);

/** This method integrates any given function over the interior of mrep */
double IntegrateFunctionOverInterior (
  MedialAtomGrid *xGrid, SMLVec3d *xPoints, 
  double *xWeights, size_t nCuts, EuclideanFunction *fMatch);

/** Integrate any function over the interior of a medial atom */
double IntegrateFunctionOverBoundary (
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, 
  double *xWeights, EuclideanFunction *fMatch);

/** Compute gradient of a function over the boundary */
double ComputeFunctionJetOverBoundary(
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, 
  double *xWeights, EuclideanFunction *fMatch, SMLVec3d *xOutGradient);


/** 
 * Integrate a function over the boundary, making sure that irregular quads are 
 * super-sampled. This results in better quality integration, but is more costly
 * to compute.
 */
double AdaptivelyIntegrateFunctionOverBoundary(
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, 
  double xMinQuadArea, EuclideanFunction *fMatch);


#endif
