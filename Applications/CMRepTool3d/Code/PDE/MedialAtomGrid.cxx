#include "MedialAtomGrid.h"
#include <iostream>

using namespace std;

/**
 * Cartesian Atom Iterator
 */
class CartesianMedialAtomIterator : public MedialAtomIterator
{
public:
  CartesianMedialAtomIterator(CartesianMedialAtomGrid *grid)
    {
    this->grid = grid;
    m = grid->GetGridSize(0);
    n = grid->GetGridSize(1);
    GoToBegin();
    }

  size_t GetIndex()
    { return k; }

  MedialAtomIterator &operator++()
    { 
    k++; i++; 
    if(i >= m) 
      { j++; i = 0; } 
    return *this;
    } 
  
  bool IsAtEnd()
    { return j >= n; }
  
  void GoToBegin()
    { i = 0; j = 0; k = 0; }

  bool IsEdgeAtom()
    { return i == 0 || j == 0 || i == m - 1 || j == n - 1; }

private:
  // Pointer to the grid
  CartesianMedialAtomGrid *grid;

  // Indices and sizes
  size_t i, j, m, n, k;
};

/**
 * Cartesian Quad Iterator
 */
class CartesianMedialQuadIterator : public MedialQuadIterator
{
public:
  CartesianMedialQuadIterator(CartesianMedialAtomGrid *grid)
    {
    this->grid = grid;
    i = j = 0;
    m = grid->GetGridSize(0);
    n = grid->GetGridSize(1);
    }
    
  size_t GetAtomIndex(size_t di, size_t dj)
    { return grid->GetAtomIndex(i + di, j + dj); }

  MedialQuadIterator &operator++()
    { 
    if(i == m-2) 
      { j++; i = 0; } 
    else i++; 
    return *this;
    }
    
  bool IsAtEnd()
    { return j >= n-1; }
  
  void GoToBegin()
    { i = j = 0; }

private:
  // The grid
  CartesianMedialAtomGrid *grid;

  // Indices and sizes
  size_t i, j, m, n;
};

/**
 * An iterator for parsing boundary points associated with a medial grid
 */
class CartesianMedialBoundaryPointIterator 
: public MedialBoundaryPointIterator
{
public:
  CartesianMedialBoundaryPointIterator(CartesianMedialAtomGrid *grid)
    : itAtom(grid) 
    { this->GoToBegin(); }
  
  size_t GetIndex()
    { return iIndex; }

  size_t GetOppositeIndex()
    {
    if(itAtom.IsEdgeAtom())
      return iIndex;
    else if(iSide == 0)
      return iIndex + 1;
    else return iIndex - 1;
    }

  size_t GetAtomIndex()
    { return itAtom.GetIndex(); }

  bool IsEdgeAtom()
    { return itAtom.IsEdgeAtom(); }

  size_t GetBoundarySide()
    { return iSide; }

  MedialBoundaryPointIterator &operator++() 
    {
    if(itAtom.IsEdgeAtom() || iSide == 1)
      { ++itAtom; iSide = 0; }
    else
      iSide = 1;
    iIndex++;
    }
  
  bool IsAtEnd()
    { return itAtom.IsAtEnd(); }
  
  void GoToBegin() 
    { itAtom.GoToBegin(); iIndex = 0; iSide = 0; }

private:
  CartesianMedialAtomIterator itAtom;
  size_t iIndex, iSide;
};

/**
 * An iterator for parsing quad cells in a medial atom grid
 */
class CartesianMedialBoundaryQuadIterator : public MedialBoundaryQuadIterator
{
public:
  CartesianMedialBoundaryQuadIterator(CartesianMedialAtomGrid *grid)
    : itQuad(grid)
    { this->grid = grid; this->GoToBegin(); }
  
  // Get index for one of the four points associated with this quad
  // (i and j are between 0 and 1 both). To make sure the quads are equally
  // winded, we swap i and j if iSide is 1
  size_t GetBoundaryIndex(size_t i, size_t j)
    { return grid->GetBoundaryPointIndex(GetAtomIndex(i, j), iSide); }

  size_t GetAtomIndex(size_t i, size_t j)
    { 
    return (iSide == 0) 
      ? itQuad.GetAtomIndex(i, j) : itQuad.GetAtomIndex(j, i); 
    }

  size_t GetBoundarySide()
    { return iSide; }

  MedialBoundaryQuadIterator &operator++()
    {
    if(iSide == 0) 
      { iSide = 1; }
    else 
      { iSide = 0; ++itQuad; }
    iIndex++;
    return *this;
    }
  
  bool IsAtEnd()
    { return itQuad.IsAtEnd(); }
  
  void GoToBegin()
    { itQuad.GoToBegin(); iSide = 0; iIndex = 0; }
  
private:
  CartesianMedialQuadIterator itQuad;
  CartesianMedialAtomGrid *grid;
  size_t iSide, iIndex;
};

/** An iterator that goes through the internal points in an m-rep */
class CartesianMedialInternalPointIterator : public MedialInternalPointIterator
{
public:
  CartesianMedialInternalPointIterator(CartesianMedialAtomGrid *grid, size_t nCuts)
    : itAtom(grid), itBnd(grid);
    {
    // Max depth is nCuts + 1, i.e. point index counting from 0 (medial atom)
    // to -maxdepth or +maxdepth
    iMaxDepth = nCuts + 1;
    GoToBegin();
    }
  
  size_t GetIndex()
    { return iIndex; }

  size_t GetAtomIndex()
    { 
    if(depth == 0) return itAtom.GetIndex(); 
    else return itBnd.GetAtomIndex();
    }
  
  int GetMaxDepth() 
    { return iMaxDepth; }
  
  /** Get the signed depth of the point; 0 means on the medial axis */
  int GetDepth() 
    { 
    if(iDepth == 0) return 0;
    if(itBnd->GetBoundarySide() == 0) return -iDepth;
    else return iDepth;
    }

  // We first iterate over the medial atoms, then over the atoms closest to
  // the medial axis, and so on, until we get to the boundary
  MedialBoundaryQuadIterator &operator++()
    {
    if(depth == 0)
      if(itAtom.IsAtEnd())
        depth = 1;
      else
        ++itAtom;
    else
      if(itBnd.IsAtEnd())
        { iDepth++; itBnd.GoToBegin(); }
      else
        ++itBnd;
    iIndex++;
    }
  
  bool IsAtEnd() 
    { return iDepth > iMaxDepth; } 
  
  void GoToBegin()
    {
    itAtom.GoToBegin(); itBnd.GoToBegin();
    iDepth = 0; iIndex = 0;
    }
    
private:
  // We need both an atom and a boundary point iterator
  CartesianMedialAtomIterator itAtom;
  CartesianMedialBoundaryPointIterator itBnd;
  
  size_t iIndex;
  int iDepth, iMaxDepth;
};

/** An iterator that goes through the internal cells in an m-rep */
class CartesianMedialInternalCellIterator : public MedialInternalCellIterator
{
public:
  CartesianMedialInternalCellIterator(CartesianMedialAtomGrid *grid, size_t nCuts) 
    : itQuad(grid)
    {
    this->grid = grid;
    iMaxDepth = nCuts;
    GoToBegin();
    }
    
  size_t GetAtomIndex(size_t i, size_t j)
    { return itQuad.GetAtomIndex(i, j); }

  size_t GetInternalPointIndex(size_t i, size_t j, size_t k)
    { return grid->GetInternalPointIndex(GetAtomIndex(i, j), GetDepth(k)); }
    
  /** Get the depth of the side k (0 inner, 1 outer) */
  int GetDepth(size_t k)
    {
    int depth = k + iDepth;
    return (itQuad.GetBoundarySide() == 0) ? -depth : depth;
    }
  
  MedialBoundaryQuadIterator &operator++()
    {
    if(itQuad.IsAtEnd()) 
      { ++iDepth; itQuad.GoToBegin(); }
    else ++itQuad;
    }

  bool IsAtEnd()
    { return iDepth > iMaxDepth; }
  
  void GoToBegin()
    {
    iDepth = 0;
    itQuad->GoToBegin();
    }

private:
  CartesianMedialAtomGrid *grid;
  CartesianMedialBoundaryQuadIterator itQuad;
  size_t iDepth;
};

/**
 * Cartesian Medial Atom Grid
 */
CartesianMedialAtomGrid::CartesianMedialAtomGrid(size_t m, size_t n)
: xIndex(n+1,0)
{
  xIndex.resize(n + 1);
  this->m = m; this->n = n;
  for(size_t j = 0, k = 0; j <= n; j++, k+=m)
    xIndex[j] = k;

  // Compute the numbers of things
  nInner = (m - 2) * (n - 2);
  nEdge = m * n - nInner;
  
  numAtoms = m * n;
  numQuads = (m - 1) * (n - 1);
  numBndPts = 2 * nInner + nEdge;
  numBndQuads = numQuads * 2;

  // Compute the boundary site mapping
  xBndIndex.resize(2 * m * n);
  CartesianMedialBoundaryPointIterator itBnd(this);
  while(!itBnd.IsAtEnd())
    {
    size_t i = 2 * itBnd.GetAtomIndex();
    xBndIndex[i] = itBnd.GetIndex();
    xBndIndex[i+1] = itBnd.GetOppositeIndex();
    if(!itBnd.IsEdgeAtom())
      ++itBnd;
    ++itBnd;
    }
}

MedialAtomIterator *CartesianMedialAtomGrid::NewAtomIterator()
{
  return new CartesianMedialAtomIterator(this);
}

MedialQuadIterator *CartesianMedialAtomGrid::NewQuadIterator()
{
  return new CartesianMedialQuadIterator(this);
}

MedialBoundaryPointIterator *CartesianMedialAtomGrid::NewBoundaryPointIterator()
{
  return new CartesianMedialBoundaryPointIterator(this);
}

MedialBoundaryQuadIterator *CartesianMedialAtomGrid::NewBoundaryQuadIterator()
{
  return new CartesianMedialBoundaryQuadIterator(this);
}

/*
MedialCellIterator *CartesianMedialAtomGrid::NewCellIterator(size_t depth)
{
  return new CartesianMedialCellIterator(this, depth);
}

MedialInternalPointIterator *
CartesianMedialAtomGrid::NewInternalPointIterator(size_t depth)
{
  return new CartesianMedialInternalPointIterator(this, depth);
}
*/

void TestCartesianGrid()
{
  CartesianMedialAtomGrid grid(5, 3);
  
  CartesianMedialAtomIterator itAtom(&grid);
  while(!itAtom.IsAtEnd())
    {
    cout << "Atom " << itAtom.GetIndex() << "; edge: " << itAtom.IsEdgeAtom() << endl;
    ++itAtom;
    }

  CartesianMedialQuadIterator itQuad(&grid);
  while(!itQuad.IsAtEnd())
    {
    cout << "Quad " 
      << itQuad.GetAtomIndex(0, 0) << ", "
      << itQuad.GetAtomIndex(0, 1) << ", "
      << itQuad.GetAtomIndex(1, 1) << ", "
      << itQuad.GetAtomIndex(1, 0) << endl;
    ++itQuad;
    }

  CartesianMedialBoundaryPointIterator itBnd(&grid);
  while(!itBnd.IsAtEnd())
    {
    cout << "Bnd " << itBnd.GetIndex() << "; Opp: " << itBnd.GetOppositeIndex() 
      << "; Atom: " << itBnd.GetAtomIndex() << "; Edge: " << itBnd.IsEdgeAtom()
      << "; Side: " << itBnd.GetBoundarySide() << endl;
    ++itBnd;
    }

  CartesianMedialBoundaryQuadIterator itBQuad(&grid);
  while(!itBQuad.IsAtEnd())
    {
    cout << "BQuad " 
      << itBQuad.GetBoundaryIndex(0, 0) << ", "
      << itBQuad.GetBoundaryIndex(0, 1) << ", "
      << itBQuad.GetBoundaryIndex(1, 1) << ", "
      << itBQuad.GetBoundaryIndex(1, 0) << "; Atoms: " 
      << itBQuad.GetAtomIndex(0, 0) << ", "
      << itBQuad.GetAtomIndex(0, 1) << ", "
      << itBQuad.GetAtomIndex(1, 1) << ", "
      << itBQuad.GetAtomIndex(1, 0) << ", Side: " 
      << itBQuad.GetBoundarySide() << endl;
    ++itBQuad;
    }
}




