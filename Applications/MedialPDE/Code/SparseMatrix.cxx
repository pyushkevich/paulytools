#include "SparseMatrix.h"
#include "SparseMatrix.txx"
#include <vnl/vnl_sparse_matrix.txx>

// Need an instantiation of sparse matrix
template vnl_sparse_matrix<int>;
template ImmutableSparseArray<double>;
template ImmutableSparseArray<int>;
template ImmutableSparseMatrix<double>;
template ImmutableSparseMatrix<int>;
