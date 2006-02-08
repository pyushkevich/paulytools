/**
 * This class represents a sparse matrix whose zero elements can not
 * be altered once they have been set. This matrix can only be 
 * initialized using a mutable sparse matrix, such as those provided 
 * by VNL
 */
#ifndef __SparseMatrix_h_
#define __SparseMatrix_h_

#include <vnl/vnl_sparse_matrix.h>
#include <iostream>

template<class TVal>
class ImmutableSparseMatrix
{
public:
  // Typedefs
  typedef ImmutableSparseMatrix<TVal> Self;
  typedef vnl_sparse_matrix<TVal> SourceType;

  // Default constructor
  ImmutableSparseMatrix();

  // Destructor
  ~ImmutableSparseMatrix();

  // Assignment operator for VNL
  Self& operator= (SourceType &src);

  // Row iterator goes through nonzero elements of rows
  class RowIterator {
  public:
    RowIterator(Self *p, size_t row)
      { this->p = p; iPos = p->xRowIndex[row]; iEnd = p->xRowIndex[row+1]; }
    
    bool IsAtEnd()
      { return iPos == iEnd; }

    TVal &Value()
      { return p->xSparseValues[iPos]; }

    size_t Column()
      { return p->xColIndex[iPos]; }

    RowIterator &operator++()
      { ++iPos; return *this; }
    
  private:
    Self *p;
    size_t iPos, iEnd;
  };
    
  // Get the row iterator
  RowIterator Row(size_t iRow)
    { return RowIterator(this, iRow); } 

  // Compute the matrix product C = A * B
  static void Multiply(Self &C, const Self &A, const Self &B);

  // Compare two matrices
  bool operator == (const Self &B);

  // Print to the standard stream
  void PrintSelf(std::ostream &out) const;

private:
  
  /** Representation for the sparse matrix */
  TVal *xSparseValues;
  size_t *xRowIndex, *xColIndex, nRows, nColumns, nSparseEntries;

  // Method to reinitialize storage if needed
  void Reset();

  // Iterator is our friend
  friend class RowIterator;
};

// Print the matrix to an output stream
template<class TVal>
std::ostream& operator << (std::ostream &out, const ImmutableSparseMatrix<TVal> &A)
  { A.PrintSelf(out); }

#endif
