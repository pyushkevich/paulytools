#ifndef __MedialAtom_h_
#define __MedialAtom_h_

#include <smlmath.h>
#include <vnl/vnl_cross.h>

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

  // The index of the atom in a rectangular grid
  size_t uIndex, vIndex;
  
  // The position on the medial surface and corresponding partial derivatives
  SMLVec3d X, Xu, Xv, Xuu, Xuv, Xvv;

  // The differential geometry descriptor of the medial surface
  GeometryDescriptor G;

  // The phi function and its first partial derivatives
  double F, Fu, Fv;

  // The radius function and its partial derivatives
  double R; //, Ru, Rv, Ruu, Ruv, Rvv;

  // Area element = sqrt(g)
  double aelt;

  // The normal vector and the Riemannian gradient of R on the surface
  SMLVec3d N, xGradR, xGradPhi;

  // The Riemannian laplacian of R
  double xLapR;

  // The badness of the atom; the term sqrt(1-xGradRMagSqr)
  double xGradRMagSqr, xNormalFactor;

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
  bool ComputeBoundaryAtoms(bool flagEdgeAtom);

  /** An alternative method for subdivision surface atoms, where the gradient
   * of phi is passed in rather than calculated from Cu and Cv */
  bool ComputeBoundaryAtoms(const SMLVec3d &gradPhi, bool flagEdgeAtom);
};

/**
 * This function performs an operation C = A + p*B on medial atoms
 */
void AddScaleMedialAtoms(
  const MedialAtom &A, const MedialAtom &B, double p, MedialAtom &C);

void MedialAtomCentralDifference(
  const MedialAtom &A, const MedialAtom &B, double eps, MedialAtom &C);

#endif
