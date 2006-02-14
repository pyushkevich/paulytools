#ifndef __GenericMedialModel_h_
#define __GenericMedialModel_h_

#include "MedialAtom.h"
#include "MedialAtomIterators.h"

/**
 * This class represents a generic medial model, whether Cartesian-based or
 * subdivision surface based. The user can inquire basic information about the
 * model, such as the number of atoms, etc.
 */
class GenericMedialModel
{
public:
  
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

  /** Get the boundary point corresponding to the current position of a
   * boundary point iterator */
  BoundaryAtom &GetBoundaryPoint(const MedialBoundaryPointIterator &it)
    { return xAtoms[it.GetAtomIndex()].xBnd[it.GetBoundarySide()]; }

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
