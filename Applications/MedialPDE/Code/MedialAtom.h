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
  double F, Fu, Fv, Fuu, Fuv, Fvv;

  // The radius function and its partial derivatives
  double R, Ru, Rv; 

  // Area element = sqrt(g)
  double aelt;

  // Gaussian and mean curvature
  double xMeanCurv, xGaussCurv;

  // The normal vector and the Riemannian gradient of R on the surface
  SMLVec3d N, xGradR; //, xGradPhi;

  // The Riemannian laplacian of R
  double xLapR;

  // The magnitude of gradR and the term sqrt(1-xGradRMagSqr)
  double xGradRMagSqr, xNormalFactor;

  // Whether this is a 'crest' atom, and whether it's valid at all
  bool flagCrest, flagValid;

  // A 'user' int that can be used to flag atoms
  int user_data;

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

  /**
   * These terms are used to compute Gateaux derivatives of the atom's
   * boundary nodes with respect to different variations
   */
  struct DerivativeTerms {
    double x1_2R, Ru, Rv, x1_2F, Ru_R, Rv_R, Ru_2F, Rv_2F;
    double g1iRi, g2iRi;
    SMLVec3d N_2g, N_2nt, Xu_aelt, Xv_aelt;
    double inv_mag_gradR;
  };

  /** 
   * Compute directional derivative of the contravariant tensor given the
   * directional derivatives of X, Xu and Xv contained in dAtom
   */
  void ComputeMetricTensorDerivatives(MedialAtom &dAtom) const;

  /** Compute the derivatives of the Christoffel symbols */
  void ComputeChristoffelDerivatives(MedialAtom &dAtom) const;

  /**
   * Compute terms that are common in derivative computations (i.e. terms that
   * are the same when computing the Gateaux derivative of the boundary atoms
   * with respect to different variations in X and phi
   */
  void ComputeCommonDerivativeTerms(DerivativeTerms &dt) const;

  /** 
   * Compute the derivatives of the boundary atoms with respect to a
   * variation. The method will not affect the current atom but will modify
   * the terms of the 'derivative' Atom passed in as the parameter
   */
  void ComputeBoundaryAtomDerivatives(MedialAtom &dAtom, const DerivativeTerms &dt) const;

  /**
   * Set all the derivative terms in the atom to equal zero. This is useful when the
   * atom is unaffected by some variation
   */
  void SetAllDerivativeTermsToZero();

  MedialAtom();
};

/**
 * This function performs an operation C = A + p*B on medial atoms
 */
void AddScaleMedialAtoms(
  const MedialAtom &A, const MedialAtom &B, double p, MedialAtom &C);

void MedialAtomCentralDifference(
  const MedialAtom &A, const MedialAtom &B, double eps, MedialAtom &C);

#endif
