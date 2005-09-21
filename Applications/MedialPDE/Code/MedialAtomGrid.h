#ifndef __MedialAtomGrid_h_
#define __MedialAtomGrid_h_

#include <vector>
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

#endif
