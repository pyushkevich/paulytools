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
#include <list>
#include <vector>

template<class TVal>
class ImmutableSparseArray
{
public:
  // Typedefs
  typedef ImmutableSparseArray<TVal> Self;

  // Typedef for import from VNL
  typedef vnl_sparse_matrix<TVal> VNLSourceType;

  // Typedefs for import from STL structures
  typedef std::pair<size_t, TVal> STLEntryType;
  typedef std::list<STLEntryType> STLRowType;
  typedef std::vector<STLRowType> STLSourceType;

  // Default constructor
  ImmutableSparseArray();

  // Destructor
  ~ImmutableSparseArray();

  // Assignment operator for VNL
  void SetFromVNL(VNLSourceType &src);

  // Assignment operator that takes an array of lists
  void SetFromSTL(STLSourceType &src, size_t nColumns);

  // Set all the arrays in the matrix to external pointers. The caller must
  // relinquish the control of the pointers to the sparse matrix, which will
  // delete the data at some point
  void SetArrays(size_t rows, size_t cols, size_t *xRowIndex, size_t *xColIndex, TVal *data);

  // Pointers to the data stored inside the matrix
  size_t *GetRowIndex() { return xRowIndex; }
  size_t *GetColIndex() { return xColIndex; }
  TVal *GetSparseData() { return xSparseValues; }

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

  // Set i-th non-zero value in a row
  TVal &GetValueBySparseIndex(size_t iRow, size_t iNZInRow)
    { return xSparseValues[xRowIndex[iRow] + iNZInRow]; }

  // Get the number of sparse values
  size_t GetNumberOfSparseValues() { return nSparseEntries; }

  // Get the number of rows and columns
  size_t GetNumberOfRows() { return nRows; }
  size_t GetNumberOfColumns() { return nColumns; }  

protected:
  
  /** Representation for the sparse matrix */
  TVal *xSparseValues;
  size_t *xRowIndex, *xColIndex, nRows, nColumns, nSparseEntries;

  // Method to reinitialize storage if needed
  void Reset();

  // Iterator is our friend
  friend class RowIterator;
};


template<class TVal>
class ImmutableSparseMatrix : public ImmutableSparseArray<TVal>
{
public:
  // Typedefs
  typedef ImmutableSparseMatrix<TVal> Self;

  // Compare two matrices
  bool operator == (const Self &B);

  // Compute the matrix product C = A * B
  static void Multiply(Self &C, const Self &A, const Self &B);

  // Print to the standard stream
  void PrintSelf(std::ostream &out) const;
};

// Print the matrix to an output stream
template<class TVal>
std::ostream& operator << (std::ostream &out, const ImmutableSparseMatrix<TVal> &A)
  { A.PrintSelf(out); }


  
#endif
