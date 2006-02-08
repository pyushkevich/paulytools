#include "SparseMatrix.h"
#include <algorithm>
#include <vnl/vnl_sparse_matrix.txx>

template<class TVal>
ImmutableSparseMatrix<TVal>
::ImmutableSparseMatrix()
{
  xSparseValues = NULL;
  xRowIndex = NULL;
  xColIndex = NULL;
  nRows = nColumns = nSparseEntries = 0;
}

template<class TVal>
ImmutableSparseMatrix<TVal>
::~ImmutableSparseMatrix()
{
  Reset();
}

template<class TVal>
void
ImmutableSparseMatrix<TVal>
::Reset()
{
  nRows = nColumns = nSparseEntries = 0;
  if(xSparseValues)
    { 
    delete xSparseValues; 
    delete xRowIndex; 
    delete xColIndex; 
    xSparseValues = NULL;
    xRowIndex = xColIndex = NULL;
    }
}

template<class TVal>
ImmutableSparseMatrix<TVal> &
ImmutableSparseMatrix<TVal>
::operator = (SourceType &src)
{
  size_t i, j;
  
  // Delete the sparse matrix storage if it exists
  Reset();

  // Set the number of rows and columns
  nRows = src.rows();
  nColumns = src.columns();

  // Allocate the row index (number of rows + 1)
  xRowIndex = new size_t[src.rows() + 1];

  // Fill the row index with indices into the sparse data
  xRowIndex[0] = 0;
  for(i = 0; i < src.rows(); i++)
    xRowIndex[i+1] = xRowIndex[i] + src.get_row(i).size();

  // Set the number of non-zero elements
  nSparseEntries = xRowIndex[src.rows()];

  // Initialize the data and column index arrays
  xColIndex = new size_t[nSparseEntries];
  xSparseValues = new TVal[nSparseEntries];

  // Fill the arrays
  size_t k = 0;
  for(i = 0; i < src.rows(); i++)
    {
    typename SourceType::row &r = src.get_row(i);
    for(j = 0; j < r.size(); j++, k++)
      {
      xColIndex[k] = r[j].first;
      xSparseValues[k] = r[j].second;
      }
    }
}

template<class TVal>
bool
ImmutableSparseMatrix<TVal>
::operator == (const Self &B)
{
  // Metastructure must match
  if(nColumns != B.nColumns || nRows != B.nRows 
    || nSparseEntries != B.nSparseEntries)
    return false;

  // Compare row indices, etc
  for(size_t i = 0; i < nRows; i++)
    {
    // Row size must match
    if(xRowIndex[i+1] != B.xRowIndex[i+1])
      return false;
    
    // Column entries and values must match
    for(size_t j = xRowIndex[i]; j < xRowIndex[i+1]; j++)
      if(xColIndex[j] != B.xColIndex[j] ||
        xSparseValues[j] != B.xSparseValues[j])
        return false;
    }

  return true;
}


template<class TVal>
void 
ImmutableSparseMatrix<TVal>
::Multiply(Self &C, const Self &A, const Self &B)
{
  size_t i, j, k, l, q;

  // Of course, check compatibility
  assert(A.nColumns == B.nRows);

  // This is a horrible cheat, but we will create a vnl_sparse_matrix
  // into which we will stick in intermediate products
  vnl_sparse_matrix<TVal> T(A.nRows, B.nColumns);

  for(i = 0; i < A.nRows; i++) 
    for(j = A.xRowIndex[i]; j < A.xRowIndex[i+1]; j++)
      {
      // Here we are looking at the element A(ik). Its product
      // with B(kq) contributes to C(iq). So that's what we do
      k = A.xColIndex[j];

      // Loop over data in B(k*)
      for(l = B.xRowIndex[k]; l < B.xRowIndex[k+1]; l++)
        {
        // Here we found a non-zero element B(kq).
        q = B.xColIndex[l];

        // Add the product to C(iq)
        T(i, q) += A.xSparseValues[j] * B.xSparseValues[l];
        }
      }
 
  // Now, just use the assignment operator to compact the sparse matrix
  C = T;
}

template<class TVal>
void
ImmutableSparseMatrix<TVal>
::PrintSelf(std::ostream &out) const 
{
  size_t i, j;
  out << "ImmutableSparseMatrix: [ ";
  for(i = 0; i < nRows; i++) 
    for(j = xRowIndex[i]; j < xRowIndex[i+1]; j++)
      out << "(" << i << "," << xColIndex[j] << "," << xSparseValues[j] << ") ";
  out << "]";
}

// Need an instantiation of sparse matrix
template vnl_sparse_matrix<int>;
template ImmutableSparseMatrix<double>;
template ImmutableSparseMatrix<int>;

