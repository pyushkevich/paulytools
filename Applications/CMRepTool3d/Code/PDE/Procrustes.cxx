#include "Procrustes.h"
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_trace.h>

typedef vnl_matrix<double> Mat;
typedef vnl_vector<double> Vec;

/**
 * Compute procrustes match from A to B 
 */
void ExtendedProcrustesAnalysis(const Mat &A, const Mat &B, 
  Mat &Ahat, Mat &R, Vec &t, double &s)
{
  // Get the data dimensions
  size_t k = A.cols(), p = A.rows();
  
  // Initialize the one-matrix
  Mat Z(p, p);
  Z.fill(-1.0 / p);
  Z.fill_diagonal(1.0 - 1.0 / p);

  // Compute the matrix that should be SVD'd
  Mat S = A.transpose() * Z * B;
  
  // Compute the SVD of S
  vnl_svd<double> svd(S);

  // Compute the rotation matrix
  R = svd.U().transpose() * svd.V();

  // Compute the scale factor
  s = vnl_trace(R.transpose() * S) / vnl_trace(A.transpose() * Z * A);

  // Compute the translation
  t = (B - s * A * R).transpose() * Vec(p, 1.0 / p);

  // Compute the transformed shape
  Ahat = s * A * R + outer_product(Vec(p, 1.0), t);
}

/**
 * Compute the generalized procrustes alignment of a set of shapes
 */
void GeneralizedProcrustesAnalysis(size_t m, Mat *A, Mat *R, Vec *t, double *s)
{
  // Get the data dimensions
  size_t k = A[0].cols(), p = A[0].rows(), i, j;
  
  // Start with the initial Procrustes mean
  Mat C(A[0]);

  bool done = false;
  while(!done)
    {
    // Map every matrix onto the centroid
    for(i = 0; i < m; i++)
      ExtendedProcrustesAnalysis(A[i], C, A[i], R[i], t[i], s[i]);

    // Compute the new Procrustes mean
    Mat C1(p, k, 0.0);
    for(i = 0; i < m; i++)
      C1 += A[i];
    C1 /= p;

    // Compare the procrustes means
    if((C1 - C).array_inf_norm() < 0.00001)
      done = true;
    C = C1;
    }
}
