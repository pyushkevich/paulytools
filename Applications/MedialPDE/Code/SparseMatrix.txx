#include <algorithm>
#include <cassert>

template<class TVal>
ImmutableSparseArray<TVal>
::ImmutableSparseArray()
{
  xSparseValues = NULL;
  xRowIndex = NULL;
  xColIndex = NULL;
  nRows = nColumns = nSparseEntries = 0;
}

template<class TVal>
ImmutableSparseArray<TVal>
::~ImmutableSparseArray()
{
  Reset();
}

template<class TVal>
void
ImmutableSparseArray<TVal>
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
void
ImmutableSparseArray<TVal>
::SetFromSTL(STLSourceType &src, size_t nColumns)
{
  size_t i, j;
  
  // Delete the sparse matrix storage if it exists
  Reset();

  // Set the number of rows and columns
  this->nRows = src.size();
  this->nColumns = nColumns;

  // Allocate the row index (number of rows + 1)
  xRowIndex = new size_t[nRows + 1];

  // Fill the row index with indices into the sparse data
  xRowIndex[0] = 0;
  for(i = 0; i < nRows; i++)
    xRowIndex[i+1] = xRowIndex[i] + src[i].size();

  // Set the number of non-zero elements
  nSparseEntries = xRowIndex[nRows];

  // Initialize the data and column index arrays
  xColIndex = new size_t[nSparseEntries];
  xSparseValues = new TVal[nSparseEntries];

  // Fill the arrays
  size_t k = 0;
  for(i = 0; i < nRows; i++)
    {
    typename STLRowType::iterator it = src[i].begin();
    for(; it != src[i].end(); ++it, k++)
      {
      xColIndex[k] = it->first;
      xSparseValues[k] = it->second;
      }
    }

}

template<class TVal>
void
ImmutableSparseArray<TVal>
::SetFromVNL(VNLSourceType &src)
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
    typename VNLSourceType::row &r = src.get_row(i);
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
  if(this->nColumns != B.nColumns || this->nRows != B.nRows 
    || this->nSparseEntries != B.nSparseEntries)
    return false;

  // Compare row indices, etc
  for(size_t i = 0; i < this->nRows; i++)
    {
    // Row size must match
    if(this->xRowIndex[i+1] != B.xRowIndex[i+1])
      return false;
    
    // Column entries and values must match
    for(size_t j = this->xRowIndex[i]; j < this->xRowIndex[i+1]; j++)
      if(this->xColIndex[j] != B.xColIndex[j] ||
        this->xSparseValues[j] != B.xSparseValues[j])
        return false;
    }

  return true;
}

template<class TVal>
ImmutableSparseArray<TVal>::ImmutableSparseArray(const ImmutableSparseArray<TVal> &src)
{
  // Make copies of all non-arrays
  nRows = src.nRows;
  nColumns = src.nColumns;
  nSparseEntries = src.nSparseEntries;

  // If the source object is NULL, there is nothing to do
  if(src.xSparseValues == NULL) 
  {
    xSparseValues = NULL;
    xRowIndex = NULL;
    xColIndex = NULL;
    return;
  }

  // Allocate the arrays
  xRowIndex = new size_t[nRows + 1];
  xColIndex = new size_t[nSparseEntries];
  xSparseValues = new TVal[nSparseEntries];

  // Copy the array contennts
  std::copy(src.xRowIndex, src.xRowIndex + nRows + 1, xRowIndex);
  std::copy(src.xColIndex, src.xColIndex + nSparseEntries, xColIndex);
  std::copy(src.xSparseValues, src.xSparseValues + nSparseEntries, xSparseValues);
}

template<class TVal>
ImmutableSparseArray<TVal> &
ImmutableSparseArray<TVal>::operator= (const ImmutableSparseArray<TVal> &src)
{
  // Check if this is the same object (or both are reset)
  if(xSparseValues == src.xSparseValues)
    return *this;

  // Clear all the data
  Reset();

  // If the source object is NULL, there is nothing to do
  if(src.xSparseValues == NULL) return *this;

  // Make copies of all non-arrays
  nRows = src.nRows;
  nColumns = src.nColumns;
  nSparseEntries = src.nSparseEntries;

  // Allocate the arrays
  xRowIndex = new size_t[nRows + 1];
  xColIndex = new size_t[nSparseEntries];
  xSparseValues = new TVal[nSparseEntries];

  // Copy the array contennts
  std::copy(src.xRowIndex, src.xRowIndex + nRows + 1, xRowIndex);
  std::copy(src.xColIndex, src.xColIndex + nSparseEntries, xColIndex);
  std::copy(src.xSparseValues, src.xSparseValues + nSparseEntries, xSparseValues);
}


template<class TVal>
typename ImmutableSparseMatrix<TVal>::Vec 
ImmutableSparseMatrix<TVal>
::MultiplyTransposeByVector(const Vec &b)
{
  // Make sure the dimensions match
  assert(b.size() == this->nRows);

  // Initialize the vector
  Vec c(this->nColumns, 0);

  // Iterate over rows and columns of the matrix
  for(size_t i = 0; i < this->nRows; i++)
    for(size_t j = this->xRowIndex[i]; j < this->xRowIndex[i+1]; j++)
      c[this->xColIndex[j]] += this->xSparseValues[j] * b[i];

  return c;
}

template<class TVal>
typename ImmutableSparseMatrix<TVal>::Vec 
ImmutableSparseMatrix<TVal>
::MultiplyByVector(const Vec &b)
{
  // Make sure the dimensions match
  assert(b.size() == this->nColumns);

  // Initialize the vector
  Vec c(this->nRows, 0);

  // Iterate over rows and columns of the matrix
  for(size_t i = 0; i < this->nRows; i++)
    for(size_t j = this->xRowIndex[i]; j < this->xRowIndex[i+1]; j++)
      c[i] += this->xSparseValues[j] * b[this->xColIndex[j]];

  return c;
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
  C.SetFromVNL(T);
}

template<class TVal>
void
ImmutableSparseMatrix<TVal>
::PrintSelf(std::ostream &out) const 
{
  size_t i, j;
  out << "ImmutableSparseArray: [ ";
  for(i = 0; i < this->nRows; i++) 
    for(j = this->xRowIndex[i]; j < this->xRowIndex[i+1]; j++)
      out << "(" << i << "," << this->xColIndex[j] << "," << this->xSparseValues[j] << ") ";
  out << "]";
}

template<class TVal>
void
ImmutableSparseArray<TVal>
::SetArrays(size_t rows, size_t cols, size_t *xRowIndex, size_t *xColIndex, TVal *data)
{
  Reset();
  this->nRows = rows; 
  this->nColumns = cols;
  this->nSparseEntries = xRowIndex[rows];
  this->xRowIndex = xRowIndex; 
  this->xColIndex = xColIndex;
  this->xSparseValues = data;
}
