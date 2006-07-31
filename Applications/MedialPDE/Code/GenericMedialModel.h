#ifndef __GenericMedialModel_h_
#define __GenericMedialModel_h_

#include "MedialAtom.h"
#include "MedialAtomIterators.h"
#include "MedialException.h"
#include "CoarseToFineMappingDescriptor.h"
#include "AffineTransformDescriptor.h"
#include "Registry.h"
#include <iostream>
#include <vector>

class AffineTransformDescriptor;

/**
 * This class represents a generic medial model, whether Cartesian-based or
 * subdivision surface based. The user can inquire basic information about the
 * model, such as the number of atoms, etc.
 */
class GenericMedialModel
{
public:
  // Vector typedef
  typedef vnl_vector<double> Vec;
  
  size_t GetNumberOfAtoms() const
    { return xIterationContext->GetNumberOfAtoms(); }

  size_t GetNumberOfBoundaryPoints() const
    { return xIterationContext->GetNumberOfBoundaryPoints(); }

  size_t GetNumberOfInternalPoints(size_t ncuts) const
    { return xIterationContext->GetNumberOfInternalPoints(ncuts); }

  size_t GetNumberOfTriangles() const
    { return xIterationContext->GetNumberOfTriangles(); }

  size_t GetNumberOfBoundaryTriangles() const
    { return xIterationContext->GetNumberOfBoundaryTriangles(); }

  size_t GetNumberOfCells(size_t ncuts) const
    { return xIterationContext->GetNumberOfCells(ncuts); }

  size_t GetNumberOfProfileIntervals(size_t ncuts) const
    { return xIterationContext->GetNumberOfProfileIntervals(ncuts); }
  
  MedialAtomIterator GetAtomIterator() const
    { return MedialAtomIterator(xIterationContext); }

  MedialTriangleIterator GetMedialTriangleIterator() const
    { return MedialTriangleIterator(xIterationContext); }

  MedialBoundaryPointIterator GetBoundaryPointIterator() const
    { return MedialBoundaryPointIterator(xIterationContext); }

  MedialBoundaryTriangleIterator GetBoundaryTriangleIterator() const
    { return MedialBoundaryTriangleIterator(xIterationContext); }

  MedialInternalPointIterator GetInternalPointIterator(size_t ncuts) const
    { return MedialInternalPointIterator(xIterationContext, ncuts); }

  MedialInternalCellIterator GetInternalCellIterator(size_t ncuts) const
    { return MedialInternalCellIterator(xIterationContext, ncuts); }

  MedialProfileIntervalIterator GetProfileIntervalIterator(size_t ncuts) const
    { return MedialProfileIntervalIterator(xIterationContext, ncuts); }

  /** Get the array of medial atoms - child implements */
  MedialAtom *GetAtomArray() const 
    { return xAtoms; }

  /** Get the iteration context */
  MedialIterationContext *GetIterationContext() const
    { return xIterationContext; }

  /** Get the boundary point corresponding to the current position of a
   * boundary point iterator */
  BoundaryAtom &GetBoundaryPoint(const MedialBoundaryPointIterator &it)
    { return xAtoms[it.GetAtomIndex()].xBnd[it.GetBoundarySide()]; }

  /** Get the number of optimizable coefficients that define this model */
  virtual size_t GetNumberOfCoefficients() const = 0;

  /** Get the array of (optimizable) coefficients that define this model */
  virtual const Vec GetCoefficientArray() const = 0;

  /** Get one of the coefficients defining this model */
  virtual double GetCoefficient(size_t i) const = 0;

  /** Set one of the coefficients defining this model */
  virtual void SetCoefficient(size_t i, double x) = 0;

  /** Set the array of coefficients defining this model */
  virtual void SetCoefficientArray(const Vec &xData) = 0;

  /** 
   * Compute the atoms from given a set of coefficients. In PDE-based models this step
   * involves solving the PDE. The optional parameter is an array of atoms to be used as
   * the initialization (previous solution, nearby solution).
   */
  virtual void ComputeAtoms(MedialAtom *xInitialSolution = NULL) = 0;

  /**
   * This method is called before multiple calls to the ComputeJet routine. Given a  
   * of variation (direction in which directional derivatives will be computed), and 
   * an array of atoms corresponding to the variation, this method will precompute the
   * parts of the atoms that do not depend on the position in coefficient space where 
   * the gradient is taken. This method must be called for each variation in the Jet.
   */
  virtual void PrepareAtomsForVariationalDerivative(
    const Vec &xVariation, MedialAtom *dAtoms) const = 0;

  /**
   * This method computes the 'gradient' (actually, any set of directional derivatives)
   * for the current values of the coefficients. As in ComputeAtoms, a set of initial 
   * values can be passed in. By default, the previous solution is used to initialize.
   */
  virtual void ComputeAtomGradient(std::vector<MedialAtom *> &dAtoms) = 0;

  /**
   * Get a pointer to the affine transform descriptor corresponding to this class. 
   * The descriptor is a lightweight object, and its allocation should be
   * managed by the child of GenericMedialModel.
   */
  virtual const AffineTransformDescriptor *GetAffineTransformDescriptor() const = 0;

  /**
   * Get a pointer to the coarse-to-fine masking descriptor corresponding to
   * this class. This descriptor is managed internally by the child of
   * GenericMedialModel.
   */
  virtual const CoarseToFineMappingDescriptor *
    GetCoarseToFineMappingDescriptor() const = 0;

  /** Get the center of rotation for the model */
  virtual SMLVec3d GetCenterOfRotation() 
    { 
    return this->GetAffineTransformDescriptor()->GetCenterOfRotation(
      this->GetCoefficientArray());
    }


protected:

  // Hidden constructor
  GenericMedialModel()
    { xIterationContext = NULL; xAtoms = NULL; }
  
  // The context (data structure) used to facilitate iteration over atoms
  // This data structure must be initialized before iterators can be created
  MedialIterationContext *xIterationContext;

  // Medial atoms array (managed by children)
  MedialAtom *xAtoms;
};

#endif
