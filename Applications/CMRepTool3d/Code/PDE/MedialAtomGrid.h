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

  virtual int GetMaxDepth() = 0;
  virtual int GetDepth() = 0;

  // Are we on the boundary
  virtual bool IsBoundaryPoint() = 0;
  
  // Are we on the medial axis
  virtual bool IsMedialPoint() = 0;
  
  // Are we on a crest (edge of medial atom)?
  virtual bool IsEdgeAtom() = 0;

  // Get a non-negative value t that is zero on medial axis, 1 on boundary
  virtual double GetRelativeDistanceToMedialAxis() = 0;

  // This is only defined for points that are not on the medial axis and 
  // not on the crest (edge of the medial surface)
  virtual int GetBoundarySide() = 0;

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

  /** Get the depth of the side k (0 inner, 1 outer) */
  virtual int GetDepth(size_t k);
  
  virtual MedialBoundaryQuadIterator &operator++() = 0;
  virtual bool IsAtEnd() = 0;
  virtual void GoToBegin() = 0;
};

SMLVec3d GetInternalPoint(MedialInternalPointIterator *it, MedialAtom *xAtoms)
{
  // If a medial atom, return the atom
  size_t iAtom = it->GetAtomIndex();
  if(it->IsMedialPoint())
    return xAtoms[iAtom].X;

  // If a boundary point return it
  size_t iSide = it->GetBoundarySide();
  if(it->IsBoundaryPoint())
    return xAtoms[iAtom].xBnd[iSide].X;

  // Interpolate between the medial atom and the boundary one
  SMLVec3d X = xAtoms[iAtom].X, Y = xAtoms[iAtom].xBnd[iSide].X;
  return X + (Y - X) * it->GetRelativeDistanceToMedialAxis();
}

SMLVec3d GetInternalPoint(MedialInternalCellIterator *it, MedialAtom *xAtoms,
  size_t i, size_t j, size_t k)
{
  //

}

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
  // virtual MedialCellIterator *NewCellIterator(size_t depth) = 0;
  // virtual MedialInternalPointIterator *NewInternalPointIterator(size_t depth) = 0;

  // Get the number of medial atoms in this grid
  virtual size_t GetNumberOfAtoms() = 0;
  virtual size_t GetNumberOfQuads() = 0;
  virtual size_t GetNumberOfBoundaryPoints() = 0;
  virtual size_t GetNumberOfBoundaryQuads() = 0;
  virtual size_t GetNumberOfCells(size_t depth) = 0;
  virtual size_t GetNumberOfInternalPoints(size_t depth) = 0;
};

/**
 * An atom grid that lives on a cartesian unit square
 */
class CartesianMedialAtomGrid : public MedialAtomGrid
{
public:
  // Initialize, specifying the number of atoms in x and y directions
  CartesianMedialAtomGrid(size_t m, size_t n);

  size_t GetAtomIndex(size_t i, size_t j)
    { return xIndex[j] + i; }

  size_t GetBoundaryPointIndex(size_t k, size_t iSide)
    { return xBndIndex[(k << 1) + iSide]; }
  
  size_t GetGridSize(size_t d) 
    { return d==0 ? m : n; }

  size_t GetInternalPointIndex(size_t iAtom, int iDepth, int nCuts)
    {
    // Get the internal point associated with this atom. Internal points are
    // folded in the following way: first we have all the medial surface points
    // themselves; then for each boundary point, we have a list of nCuts+1 * 2
    if(iDepth == 0)
      return iAtom;
    else if(iDepth > 0)
      return numAtoms + iDepth * numBndPts + GetBoundaryPointIndex(iAtom, 1);
    else
      return numAtoms - iDepth * numBndPts + GetBoundaryPointIndex(iAtom, 0);
    }

  // Create various iterators
  MedialAtomIterator *NewAtomIterator();
  MedialQuadIterator *NewQuadIterator();
  MedialBoundaryPointIterator *NewBoundaryPointIterator();
  MedialBoundaryQuadIterator *NewBoundaryQuadIterator();
  MedialCellIterator *NewInternalCellIterator(size_t nCuts);
  MedialInternalPointIterator *NewInternalPointIterator(size_t nCuts);

  // Get the number of iterated objects 
  size_t GetNumberOfAtoms()
    { return numAtoms; }

  size_t GetNumberOfQuads()
    { return numQuads; }

  size_t GetNumberOfBoundaryPoints()
    { return numBndPts; }

  size_t GetNumberOfBoundaryQuads()
    { return numBndQuads; }

  size_t GetNumberOfCells(size_t nCuts)
    { return (nCuts + 1) * 2 * numQuads; }

  size_t GetNumberOfInternalPoints(size_t nCuts)
    { return (nCuts + 1) * GetNumberOfBoundaryPoints() + GetNumberOfAtoms(); }

private:
  // A way to map a pair of indices to a linear index
  vector<size_t> xIndex, xBndIndex;

  // Dimensions of the grid and various sizes
  size_t m, n, nEdge, nInner;
  size_t numQuads, numBndQuads, numBndPts;
};


void TestCartesianGrid();

#endif
