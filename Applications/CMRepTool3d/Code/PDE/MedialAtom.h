#ifndef __MedialAtom_h_
#define __MedialAtom_h_

#include <smlmath.h>

#include "MedialAtomGrid.h"
#include "GeometryDescriptor.h"

/** 
 * A Boundary atom
 */
struct BoundaryAtom 
{
  SMLVec3d X, N;
};

/**
 * The new medial atom representation
 */
struct MedialAtom
{
  // The coordinates of the atom in the domain
  double u, v;
  
  // The position on the medial surface and corresponding partial derivatives
  SMLVec3d X, Xu, Xv, Xuu, Xuv, Xvv;

  // The differential geometry descriptor of the medial surface
  GeometryDescriptor G;

  // The radius function and its partial derivatives
  double R, Ru, Rv, Ruu, Ruv, Rvv;

  // The normal vector and the Riemannian gradient of R on the surface
  SMLVec3d N, xGradR;

  // The Riemannian laplacian of R
  double xLapR;

  // The badness of the atom
  double xBadness;

  // Whether this is a 'crest' atom, and whether it's valid at all
  bool flagCrest, flagValid;

  // The two associated boundary 'atoms'
  BoundaryAtom xBnd[2];

  /** Compute the differential geometric quantities from X and its 1st and 2nd
   * derivatives. This does not compute the normal vector */
  void ComputeDifferentialGeometry();

  /** Compute the normal vector */
  void ComputeNormalVector();

  /** Given the differential geometry and the normal vector are computed,
   * compute GradR and the boundary sites. */
  bool ComputeBoundaryAtoms();
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

/**
 * This function computes the area associated with each triangle on the
 * boundary of the object and assigns a third of this area to each of the 
 * vertices of the triangle. Thus for each vertex on the boundary a weight
 * is generated. The total of these weights, equal to the area of the boundary
 * surface is the return value 
 */
double ComputeMedialBoundaryAreaWeights( 
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, double *xWeights);

/** 
 * Interpolate a list of internal medial points from the m-rep.
 */
void ComputeMedialInternalPoints(
  MedialAtomGrid *xGrid, MedialAtom *xAtoms, size_t nCuts, SMLVec3d *xPoints);

// This is a measure that can be computed over a volume (just a R3 function)
class EuclideanFunction {
public:
  virtual double Evaluate(const SMLVec3d &x) = 0;
  virtual SMLVec3d ComputeGradient(const SMLVec3d &x)
    { return SMLVec3d(0.0); }
};

/**
 * Compute an area of a triangle
 */
inline double TriangleArea(const SMLVec3d &A, const SMLVec3d &B, const SMLVec3d &C)
{
  return 0.5 * cross_3d(B - A, C - A).magnitude();
}

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
