/******************************************************************
 * MATRIX Library                                                 *
 ******************************************************************
 * Author:                      Paul Yushkevich
 *
 * Date:                            Apr 1, 1999
 *
 * Description                  Basic and numerical matrix operations   
 *                                  See http://www.cs.unc.edu/~pauly/matrix
 *                                  
 *  Sources:                        Uses my compilation of CLAPACK
 *                                  
 * Dependencies:                CLAPACK
 ******************************************************************
 * matrix.cpp
 *  -----------
 * Definitions for matrix and vector classes
 ******************************************************************/
#include <matrix.h>
#include <memory.h>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <stdarg.h>

using namespace std;

// Begin namespace
NAMESPACE_PAULY_START

/***********************************************************************
 * CLAPACK IMPORTS
 * These functions are not present in Intel Extensions
 ***********************************************************************/

extern "C" {
  //int dgesv_(long *n, long *nrhs, double *a, long *lda, long *ipiv, double *b, long *ldb, long *info);
  //int dposv_(char *uplo, long *n, long *nrhs, double *a, long *lda, double *b, long *ldb, long *info);
  //int dsyev_(char *jobz, char *uplo, long *n, double *a,long *lda, double *w, double *work, long *lwork, long *info);
  //int dgesvd_(char *jobu, char *jobvt, long *m, long *n, double *a, long *lda, double *s, double *u, long * ldu, double *vt, long *ldvt, double *work, long *lwork, long *info);

}


/*

 **********************************************************************
 * If you decide to use SIMD extensions, then use the following #define
 * and link with MLK libraries.  Otherwise, link with clapack.lib
 **********************************************************************

*/
#ifdef _INTEL_SIMD_ 

// Include mkl
  #include <mkl.h>

// Declare lapack_int as the int 
typedef int lapack_int;

#else

typedef long lapack_int;
/*
extern "C" {
    int dgemm_(char *transa, char *transb, long *m, long *n, long *k, 
                  double *alpha, double *a, long *lda, double *b, long *ldb, 
                  double *beta, double *c, long *ldc); 

    int dgetrf_(long *m, long *n, double *a, long *lda, long *ipiv, long *info);
   int dgetri_(long *n, double *a, long *lda, long *ipiv, double *work, long *lwork, long *info);
   int dgetrs_(char *,long *, long *, double *a, long *lda, long *ipiv, double *work, long *lwork, long *info);
    int dpotrf_(char *,long *, double *a, long *, long *);
}
*/

int dgemm_(char *transa, char *transb, long *m, long *n, long *k, 
           double *alpha, double *a, long *lda, double *b, long *ldb, 
           double *beta, double *c, long *ldc)
{
  assert(0);
  return 0;
}

int dgetrf_(long *m, long *n, double *a, long *lda, long *ipiv, long *info)
{
  assert(0);
  return 0;
}

int dgetri_(long *n, double *a, long *lda, long *ipiv, double *work, long *lwork, long *info)
{
  assert(0);
  return 0;
}

int dgetrs_(char *,long *, long *, double *a, long *lda, long *ipiv, double *work, long *lwork, long *info)
{
  assert(0);
  return 0;
}

int dpotrf_(char *,long *, double *a, long *, long *)
{
  assert(0);
  return 0;
}

int dsyev_(char *jobz, char *uplo, long *n, double *a,long *lda, double *w, double *work, long *lwork, long *info)
{
  assert(0);
  return 0;
}

int dgesvd_(char *jobu, char *jobvt, long *m, long *n, double *a, long *lda, double *s, double *u, long * ldu, double *vt, long *ldvt, double *work, long *lwork, long *info)
{
  assert(0);
  return 0;
}

  #define DGETRI dgetri_ 
  #define DGETRF dgetrf_
  #define DGETRS dgetrs_
  #define DGEMM dgemm_ 
  #define DPOTRF dpotrf_

#endif

/*
#ifdef WIN32

#else
extern "C" {
   int dposv_(char *uplo, lapack_int *n, lapack_int *nrhs, doublereal *a, lapack_int *lda, doublereal *b, lapack_int *ldb, lapack_int *info);
   int dgesv_(lapack_int *n, lapack_int *nrhs, doublereal *a, lapack_int *lda, lapack_int *ipiv, doublereal *b, lapack_int *ldb, lapack_int *info);
   int dgesvd_(char *jobu, char *jobvt, lapack_int *m, lapack_int *n, doublereal *a, lapack_int *lda, doublereal *s, doublereal *u, lapack_int *
        ldu, doublereal *vt, lapack_int *ldvt, doublereal *work, lapack_int *lwork, lapack_int *info);
   int dsyev_(char *jobz, char *uplo, lapack_int *n, doublereal *a,lapack_int *lda, doublereal *w, doublereal *work, lapack_int *lwork, lapack_int *info);
   int dgetri_(lapack_int *n, doublereal *a, lapack_int *lda, lapack_int *ipiv, doublereal *work, lapack_int *lwork, lapack_int *info);
    
    int dgemm_(char *transa, char *transb, lapack_int *m, lapack_int *n, lapack_int *k, 
        doublereal *alpha, doublereal *a, lapack_int *lda, doublereal *b, lapack_int *ldb, 
        doublereal *beta, doublereal *c, lapack_int *ldc); 
}
#endif
*/

void Matrix::initData(double *inData,bool columnMajorOrder) {
  if (inData)
    {
    if (columnMajorOrder)
      {
      memcpy(data,inData,sizeof(double)*nCells);
      }
    else
      {
      for (int iRow=0;iRow<nRows;iRow++)
        {
        for (int iCol=0;iCol<nCols;iCol++)
          {
          cols[iCol][iRow] = *inData;
          inData++;
          }
        }
      }
    }
  else
    {
    memset(data,0,sizeof(double)*nCells);
    }
}

// Initializes a matrix to new size
void Matrix::initStorage(int rows,int columns) {
  nRows = rows;
  nCols = columns;
  nCells = nRows * nCols;
  if (nCells > 0)
    {
    data = new double[nCells];
    cols = new double*[nCols];

    cols[0] = data;
    for (int i=1;i<nCols;i++)
      {
      cols[i] = cols[i-1] + nRows;
      }
    }
  else
    {
    cols = NULL;
    data = NULL;
    }
}

Matrix::Matrix(int rows,int columns,double firstValue...) {
  initStorage(rows,columns);

  // Read in values.  Values are in column major order
  va_list ap;
  va_start(ap, firstValue);

  for (int iCol=0;iCol<nCols;iCol++)
    {
    for (int iRow=0;iRow<nRows;iRow++)
      {
      if (iCol==0 && iRow==0)
        cols[iCol][iRow] = firstValue;
      else
        cols[iCol][iRow] = va_arg(ap,double);
      }
    }

  va_end(ap);
}


/*
Matrix Matrix::operator*(const Matrix &A) const {
   dassert(nCols == A.nRows);

   register int i,j,k;

   // Initialize destination matrix with all zeros
   Matrix C(nRows,A.nCols);

   // Compute product
   for(i=0;i<C.nRows;i++) {
      for(j=0;j<C.nCols;j++) {
         for(k=0;k<nCols;k++) {
            C.cols[j][i]  += cols[k][i] * A.cols[j][k];
         }
      }
   }

   return C;
}
*/

#ifdef _INTEL_SIMD_ 

void Matrix::multiply(const Matrix &A, const Matrix &B, Matrix &C) {
  static char code = 'c';
  static double alpha = 1.0;
  static double beta = 0.0;
  lapack_int m = A.nRows;
  lapack_int n = B.nCols;
  lapack_int k = A.nCols;
  cblas_dgemm(CblasColMajor,
              CblasNoTrans,CblasNoTrans,
              A.nRows,B.nCols,A.nCols,
              1.0,
              A.data,A.nRows,
              B.data,A.nCols,
              0.0,
              C.data,A.nRows);
}

#else 


void Matrix::multiply(const Matrix &A, const Matrix &B, Matrix &C) {
  register int i,j,k;

  // Compute product
  for (i=0;i<C.nRows;i++)
    {
    for (j=0;j<C.nCols;j++)
      {
      C.cols[j][i] = A.cols[0][i] * B.cols[j][0];
      for (k=1;k<A.nCols;k++)
        {
        C.cols[j][i]  += A.cols[k][i] * B.cols[j][k];
        }
      }
    }

}

#endif
/*
void Matrix::multiply(const Matrix &A, const Vector &X, Vector &R) {
    register int i,k;

   for(i=0;i<A.nRows;i++) {
        R.data[i] = A.cols[0][i] * X.data[0];
      for(k=1;k<A.nCols;k++) {
            R.data[i]  += A.cols[k][i] * X.data[k];
      }
   }
}*/

/*
void Matrix::multiply(const Matrix &A, const Matrix &B, Matrix &C) {
    static char code = 'c';
    static double alpha = 1.0;
   static double beta = 0.0;
    lapack_int m = A.nRows;
    lapack_int n = B.nCols;
    lapack_int k = A.nCols;
   dgemm_(&code,&code,&m,&n,&k,&alpha,
          A.data,&m,B.data,&k,
          &beta,C.data,&m);
}*/

Matrix Matrix::operator*(const Matrix &A) const {
  // Initialize target matrix
  Matrix C(nRows,A.nCols);

  multiply(*this,A,C);

  return C;
}


/**
 * Get transpose of a matrix
 */
Matrix Matrix::t() const {
  Matrix T(nCols,nRows);

  for (int iRow=0;iRow<T.nRows;iRow++)
    {
    for (int iCol=0;iCol<T.nCols;iCol++)
      {
      T.cols[iCol][iRow] = cols[iRow][iCol];
      }
    }

  return T;
}

int Matrix::ipFactorGenericLU(lapack_int pivot[]) {
  lapack_int n = nRows;
  lapack_int info = 0;

  DGETRF(&n,&n,data,&n,pivot,&info);
  return info;
}

int Matrix::ipFactorSPDLU() {
  lapack_int n = nRows;
  lapack_int info = 0;

  DPOTRF("u",&n,data,&n,&info);
  return info;
}


/*
int Matrix::factorGenericLU(Matrix &LU,lapack_int pivot[]) {
   // Copy our contents to the LU matrix
   LU = *this;

   // Data used in call
   lapack_int   n       = nRows;
   lapack_int   nrhs    = 0;
   double *a    = LU.data;
   lapack_int   lda = n;   
   lapack_int   info    = 0;

   // Perform LU call
    DGETRF(&n,&nrhs,a,&lda,pivot,&info);
   //dgesv_(&n,&nrhs,a,&lda,pivot,b,&ldb,&info);

   return (int) info;
}

int Matrix::factorSPDLU(Matrix &LU) {
   // Copy our contents to the LU matrix
   LU = *this;

   // Data used in call
   lapack_int n = nRows;
   lapack_int nrhs = 0;
   double *a = LU.data;
   lapack_int lda = n;   
   double *b = NULL;
   lapack_int ldb = n;
   lapack_int info = 0;

   // Perform LU call
   char uplo = 'U';

   
    DPOTRF(&uplo,&n,a,&lda,&info);

   return (int) info;
}

*/

/**
 * LU Factorization of a matrix - 'public' call
 */
int Matrix::factorLU(Matrix &L, Matrix &U, Matrix &P,int type) {
  int info;

  P.setSize(nRows,nRows);
  P = 1;

  Matrix LU = *this;

  if (type==SPD)
    {
    info = LU.ipFactorSPDLU();      
    }

  else
    {
    lapack_int *pivot = new lapack_int[nRows];
    info = LU.ipFactorGenericLU(pivot);

    // Create a permutation matrix
    for (int iRow=0;iRow<nRows;iRow++)
      {
      if (pivot[iRow] && pivot[iRow]-1!=iRow)
        {
        P.swapColumns(iRow,pivot[iRow]-1);
        }
      }

    delete pivot;
    }

  // Construct L, U
  L = LU;
  U = LU;

  for (int iRow=0;iRow<nRows;iRow++)
    {
    for (int iCol=0;iCol<iRow;iCol++)
      {
      U.cols[iCol][iRow] = 0;
      L.cols[iRow][iCol] = 0;
      }
    L.cols[iRow][iRow] = 1;
    }

  return(int) info;
}

/**
 * This method is destructive.  It solves the linear system AX = B, places the
 * result in B and the LU decomposition of A in A
 */
int Matrix::ipSolveLinearSystem(Matrix &B) {
  dassert(nRows == nCols);
  dassert(B.nRows == nRows);

  lapack_int n = nRows;
  lapack_int nrhs = B.nCols;
  lapack_int *pivot = new lapack_int[n];
  lapack_int info = 0;

  ipFactorGenericLU(pivot);
  DGETRS("N",&n,&nrhs,data,&n,pivot,B.data,&n,&info);

  delete pivot;
  return info;
}

/**
 * Perform Gaussian elimination, get back LU decomp.
 *
int Matrix::solveGE(const Matrix &B,Matrix &X,int type) {
   dassert(nRows == nCols);
   dassert(B.nRows == nRows);

   // Create a copy of ourselves
   Matrix LU;
   LU = *this;

   // Create a copy of B
   X = B;

   // Data used in call
   long n = nRows;
   long nrhs = B.nCols;
   double *a = data;
   long lda = n;   
   double *b = X.data;
   long ldb = n;
   int info = 0;
   
   // Call to gaussian eliminator
   //if(type==SPD) {
   //   char uplo = 'U';
   //   dposv_(&uplo,&n,&nrhs,a,&lda,b,&ldb,&info);
   //}
   //else {
      int *ipiv = new int[n];
        factorGenericLU(LU,ipiv);
      // dgesv_(&n,&nrhs,a,&lda,ipiv,b,&ldb,&info);
        DGETRS("N",&nRows,&X.nCols,LU.data,&nRows,ipiv,X.data,&nRows,&info);
      delete ipiv;
  // }

   return (int)info;
}
*/

// Reverse a range of rows
void Matrix::reverseRows(int startRow,int endRow) {
  while (endRow > startRow)
    swapRows(startRow++,endRow--);
}

// Reverse a range of columns
void Matrix::reverseColumns(int startColumn,int endColumn) {
  while (endColumn > startColumn)
    swapColumns(startColumn++,endColumn--);
}

// Compute the mean of this feature matrix
Vector Matrix::computeMean() const {
  // Compute the mean
  Vector mean(nCols);
  for (int c=0;c<nCols;c++)
    {
    mean(c) = 0;
    for (int r=0;r<nRows;r++)
      {
      mean(c) += cols[c][r];
      }
    mean(c) /= nRows;
    }
  return mean;
}

// Compute the mean of the matrix and a the matrix where each row has the mean
// removed
Vector Matrix::computeMean(Matrix &minusMean) const {
  Vector mean = computeMean();
  minusMean.setSize(nRows,nCols);
  for (int i=0;i<nRows;i++)
    {
    for (int j=0;j<nCols;j++)
      {
      minusMean(i,j) = cols[j][i] - mean(j);
      }
    }
  return mean;
}

// Compute the mean and covairance matrix of a feature matrix
void Matrix::computeMeanCov(Vector &mean,Matrix &cov) const {
  Matrix noMean;
  mean = computeMean(noMean);

  cov.setSize(nCols,nCols);
  cov.insertMatrix(0,0,(noMean.t() * noMean) / (nRows-1));
}

// Compute the mean and covairance matrix of a feature matrix
void Matrix::computeMeanCorr(Vector &mean,Vector &sdev,Matrix &corr) const {
  int c,r;

  // Compute the mean
  mean = computeMean();

  // Compute the standard deviation of each column
  sdev.setSize(nCols);
  for (c=0;c<nCols;c++)
    {
    double s = 0;
    for (r=0;r<nRows;r++)
      {
      double dx = cols[c][r] - mean(c);
      s += dx*dx;
      }
    sdev(c) = sqrt(s / (nRows-1));
    }

  // A Z-matrix
  Matrix Z(nRows,nCols);
  for (c=0;c<nCols;c++)
    {
    for (r=0;r<nRows;r++)
      {
      Z(r,c) = (cols[c][r] - mean(c)) / sdev(c);
      }
    }

  corr = (Z.t() * Z) / (nRows-1);
}

/**
 * This method performs PCA on a matrix whose rows are samples and columns
 * are different features.  We expect this to be a tall matrix...  
 */
int Matrix::computePCA(Vector &mean,Vector &eValues,Matrix &eVectors) const {
  // Compute the covariance matrix
  Matrix cov;
  computeMeanCov(mean,cov);

  // Compute the PCA
  int rc = cov.factorEV(eValues,eVectors);

  // Rerrange rows and columns
  eVectors.reverseAllColumns();
  eValues.reverseAllRows();

  return rc;
}

/**
 * This is a method that computes PCA, but before doing so converts the feature matrix
 * to a whitened matrix (uses correlation analysis instead of covariance)
 */
int Matrix::computeCorrelationPCA(Vector &mean,Vector &sdev,Vector &eValues,Matrix &eVectors) const {
  // Compute the covariance matrix
  Matrix corr;
  computeMeanCorr(mean,sdev,corr);

  // Compute the PCA
  int rc = corr.factorEV(eValues,eVectors);

  // Rerrange rows and columns
  eVectors.reverseAllColumns();
  eValues.reverseAllRows();

  return rc;
}

/**
 * Compute matrix determinant - Ill posed problem
 */
double Matrix::det(int type) {
  dassert(nRows == nCols);
  double det=1;

  Matrix LU = *this;

  if (type==SPD)
    {
    LU.ipFactorSPDLU();
    }
  else
    {
    lapack_int *pivot = new lapack_int[nRows];
    LU.ipFactorGenericLU(pivot);

    // Compute sign of determinant 
    for (int iRow=0;iRow<nRows;iRow++)
      {
      if (pivot[iRow] != iRow+1)
        det = 0 - det;
      }
    delete pivot;
    }

  for (int iRow=0;iRow<nRows;iRow++)
    {
    det*=LU(iRow,iRow);
    }

  return det;
}

/** 
 * Compute the logarithm of the determinant.  The sign of the determinant is returned in the 
 * second parameter
 */
double Matrix::logDeterminant(double &sign,int type) {
  dassert(nRows == nCols);
  sign=1;
  double ldet = 0;

  Matrix LU = *this;

  if (type==SPD)
    {
    LU.ipFactorSPDLU();
    }
  else
    {
    lapack_int *pivot = new lapack_int[nRows];
    LU.ipFactorGenericLU(pivot);

    // Compute sign of determinant 
    for (int iRow=0;iRow<nRows;iRow++)
      {
      if (pivot[iRow] != iRow+1)
        sign = 0 - sign;
      }
    delete pivot;
    }

  for (int iRow=0;iRow<nRows;iRow++)
    {
    double lu = LU(iRow,iRow);
    if (lu < 0)
      {
      sign = 1 - sign;
      lu = -lu;
      }
    ldet += log(lu);
    }

  return ldet;
}

// Prints a matrix
void Matrix::print()
{
  for (int iRow=0;iRow<nRows;iRow++)
    {
    for (int iCol=0;iCol<nCols;iCol++)
      {
      cout << cols[iCol][iRow] << "\t";           
      }
    cout << "\n";
    }
  cout << "\n";
}

void Matrix::printMatlab(FILE *f) const {
  fprintf(f,"[");
  for (int iRow=0;iRow<nRows;iRow++)
    {
    if (iRow)    fprintf(f,";\n");
    for (int iCol=0;iCol<nCols;iCol++)
      {
      if (iCol) fprintf(f," ");
      fprintf(f,"%.8lg",cols[iCol][iRow]);            
      }       
    }
  fprintf(f,"]\n");
}

// Prints a matrix
void Matrix::printMatlab(const char *fname) const
{
  FILE *f = fopen(fname,"wb");
  printMatlab(f);
  fclose(f);
}


/**
 * Insert matrix into this matrix with top left of insreted matrix at row,col
 */
void Matrix::insertMatrix(int row,int col,const Matrix &A) {
  dassert(row + A.nRows <= nRows);
  dassert(col + A.nCols <= nCols);

  for (int iRow=row;iRow<row + A.nRows;iRow++)
    {
    for (int iCol=col;iCol<col + A.nCols;iCol++)
      {
      cols[iCol][iRow] = A.cols[iCol-col][iRow-row];
      }
    }
}

// Extracts a sub-matrix
void Matrix::extractMatrix(int row,int col,Matrix &A)  const{
  dassert(row + A.nRows <= nRows);
  dassert(col + A.nCols <= nCols);

  for (int iRow=row;iRow<row + A.nRows;iRow++)
    {
    for (int iCol=col;iCol<col + A.nCols;iCol++)
      {
      A.cols[iCol-col][iRow-row] = cols[iCol][iRow];
      }
    }
}

/**
 * Get column vector
 */
Vector Matrix::getColumn(int col)  const{
  dassert(col < nCols);

  Vector v(nRows);

  for (int iRow=0;iRow<nRows;iRow++)
    {
    v.data[iRow] = cols[col][iRow];
    }

  return v;
}

/**
 * Get row vector
 */
Vector Matrix::getRow(int row)  const{
  dassert(row < nRows);

  Vector v(nCols);

  for (int iCol=0;iCol<nCols;iCol++)
    {
    v.data[iCol] = cols[iCol][row];
    }

  return v;
}

/**
 * Swap two rows
 */
void Matrix::swapRows(int r1,int r2) {
  dassert(r1 < nRows);
  dassert(r2 < nRows);

  for (int iCol=0;iCol<nCols;iCol++)
    {
    double tmp = cols[iCol][r1];
    cols[iCol][r1] = cols[iCol][r2];
    cols[iCol][r2] = tmp;
    }
}

/**
 * Swap two columns
 */
void Matrix::swapColumns(int c1,int c2) {
  dassert(c1 < nCols);
  dassert(c2 < nCols);

  for (int iRow=0;iRow<nRows;iRow++)
    {
    double tmp = cols[c1][iRow];
    cols[c1][iRow] = cols[c2][iRow];
    cols[c2][iRow] = tmp;
    }
}

/**
 * SVD Factorization

    The SVD is written

         A = U * SIGMA * transpose(V)

    where SIGMA is an M-by-N matrix which is zero except for its
    min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
    V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
    are the singular values of A; they are double and non-negative, and
    are returned in descending order.  The first min(m,n) columns of
    U and V are the left and right singular vectors of A.

    Note that the routine returns transpose(V), not V.
 */
int Matrix::factorSVD(Matrix &U,Matrix &Vt,Vector &sv,int type) {

  // Have to create a copy of A so it does not get destroyed
  Matrix A(*this);

  // Compute the smaller of two dimensions
  int nsv = (nRows < nCols) ? nRows : nCols;

  // Set size of U,Vt
  U.setSize(nRows,nRows);
  Vt.setSize(nCols,nCols);
  sv.setSize(nsv);

  // Input to svd method (get all rows/cols of u,vt)
  char jobu = 'A';
  char jobvt = 'A';
  long m = nRows;
  long n = nCols;
  double *a = A.data;
  long lda = n;
  double *s = sv.data;
  double *u = U.data;
  long ldu = m;
  double *vt = Vt.data;
  long ldvt = n;
  long lwork = 6*(m+n);
  double *work = new double[lwork];
  long info;

  // Call svd method
  dgesvd_(&jobu,&jobvt,&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,work,&lwork,&info);

  delete work;

  return(int) info;
}

/**
 * Last biggie - eigenvectors/eigenvalues
 * A = V * LAMBDA * V.t
 * LAMBDA contains eigenvalues of A on it diagonals
 * A must be symmetrical.
 */
int Matrix::factorEV(Vector &l,Matrix &V,int type) {
  dassert(nCols == nRows);

  V = *this;
  l.setSize(nCols);

  // Parameters
  char jobz = 'V';
  char uplo = 'U';
  long n = nCols;
  double *a = V.data;
  long lda = n;
  double *w = l.data;
  long lwork = 6*n;
  double *work = new double[lwork];
  long info = 0;

  dsyev_(&jobz, &uplo, &n, a,
         &lda, w, work, &lwork, 
         &info);

  delete work;
  return(int)info;
}

int Matrix::ipInverse() {
  dassert(nCols == nRows);

  lapack_int n = nCols;
  lapack_int *ipiv = new lapack_int[n];
  lapack_int lwork = 64 * n;
  double *work = new double[lwork];
  lapack_int info;

  // Compute inverse
  ipFactorGenericLU(ipiv);
  DGETRI(&n,data,&n,ipiv,work,&lwork,&info);

  // Delete trash
  delete[] ipiv;
  delete[] work;

  return(int)info;
}

/**
 * Oops - one more - inverse calc.
 


int Matrix::inverse(Matrix &Ainv,int type) {
   dassert(nCols == nRows);
   
   // Do LU factorization
   
   lapack_int n = nCols;
   lapack_int lda = n;
   lapack_int *ipiv = new lapack_int[n];
   lapack_int lwork = 4 * n;
   double *work = new double[lwork];
   lapack_int info;

   factorGenericLU(Ainv,ipiv);
   double *a = Ainv.data;

   // Compute inverse
   DGETRI(&n,a,&lda,ipiv,work,&lwork,&info);

   // Delete trash
   delete[] ipiv;
   delete[] work;

   return (int)info;
}
*/

// Norms of the matrix
double Matrix::oneNorm() const {
  double rtn = 0;
  for (int i=0;i<nCells;i++)
    rtn += fabs(data[i]);
  return rtn;
}

double Matrix::infinityNorm() const {
  double rtn = 0;
  for (int i=0;i<nCells;i++)
    {
    double val = fabs(data[i]);
    rtn = (val > rtn) ? val : rtn;
    }
  return rtn;
}

double Matrix::twoNorm() const {
  double rtn = 0;
  for (int i=0;i<nCells;i++)
    rtn += data[i] * data[i];
  return sqrt(rtn);
}

double Matrix::pNorm(double p) const {
  double rtn = 0;
  for (int i=0;i<nCells;i++)
    rtn += pow(fabs(data[i]),p);
  return pow(rtn,1/p);
}


Vector::Vector(int rows,double firstValue,...) : Matrix(rows,1) {
  // Read in values.  Values are in column major order
  va_list ap;
  va_start(ap, firstValue);

  data[0] = firstValue;
  for (int iRow=1;iRow<nRows;iRow++)
    {
    data[iRow] = va_arg(ap,double);
    }

  va_end(ap);
}


/* (C) Copr. 1986-92 Numerical Recipes Software 42,. */

/************************************************************************/
/* unif_rand_dbl.C -- contains routines to generate uniform random      */
/*      floats or reals in range [0.0, 1.0], and uniform random       */
/*      lapack_ints in the range [0, n], for n specified.                  */
/*                                                                      */
/*      Contents:  rand3() -- uniform random floats in range [0.0, 1.0] */
/*                 unif_rand_dbl() -- as above, but returns reals     */
/*                 randint() -- returns uniform random ints,            */
/*                      in range [0, n] for n specified.                */
/*                                                                      */
/*      Author:  rand3() is from Numerical Recipes in C, 2nd ed.        */
/*              mods and other code by A. Thall                         */
/*      Date:  2. Feb. 1997                                             */
/************************************************************************/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

// formerly, ran3(), from numeric recipes.

float rand3(long *idum)
{
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;

  if (*idum < 0 || iff == 0)
    {
    iff=1;
    mj=MSEED-(*idum < 0 ? -*idum : *idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++)
      {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
      }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++)
        {
        ma[i] -= ma[1+(i+30) % 55];
        if (ma[i] < MZ) ma[i] += MBIG;
        }
    inext=0;
    inextp=31;
    *idum=1;
    }
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return(float) mj*FAC;
}

/************************************************************************/
/* unif_rand_dbl() returns a uniform random variate between [0.0, 1.0]. */
/*      (based on NumRec routine rand3(), based itself on this-or-that  */
/*      from Knuth.  Not a linear congruence generator.                 */
/* To initialize/reinitialize, pass it a negative long int; it has a    */
/*      memory, so passing it the same initializer multiple times       */
/*      during a run of the program will produce different values.      */
/************************************************************************/
double unif_rand_dbl(long *idum)
{
  double highorder = (double) rand3(idum);
  double loworder = (double) rand3(idum);

  return highorder + loworder*FAC;
}

/************************************************************************/
/* randint() -- returns a uniformly distributed random lapack_int in the   */
/*      range [0, thismax], based on a scaled unif_rand_dbl() value.    */
/************************************************************************/
int randint(int thismax, long *seed)
{
  double scaleval = (thismax + 1)*unif_rand_dbl(seed);
  return(int) scaleval;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


Matrix *allocMatrixArray(int nMatrices,int rows,int columns) {
  Matrix *array = new Matrix[nMatrices];
  for (int i=0;i<nMatrices;i++)
    {
    array[i].setSize(rows,columns);
    }
  return array;
}


Vector *allocVectorArray(int nVectors,int size) {
  Vector *array = new Vector[nVectors];
  for (int i=0;i<nVectors;i++)
    {
    array[i].setSize(size);
    }
  return array;
}


// Begin namespace
NAMESPACE_PAULY_END


/*
void main(void) {
   long seed = 100;
   int r;
   
   Matrix A(8,8);
   Matrix B(8,8);

   for(int i=0;i<8;i++) {
      for(int j=0;j<8;j++)  {
         A.cell(i,j) = unif_rand_dbl(&seed);    
         B.cell(i,j) = unif_rand_dbl(&seed);
      }
   }
         
   Matrix L,U,P;

   A.factorLU(L,U,P);

   Matrix LU = L*U;
   LU.print();
   
   Matrix PLU = P*LU;
   PLU.print();
   
   Matrix Q = A - PLU;
   Q.print();

   double det = A.det();
   printf("%ld",det);

   Matrix Ai;
   A.inverse(Ai);
   Ai.print();

   Matrix ONE = A*Ai;
   ONE.print();

   Matrix X;
   A.solveGE(B,X);
   (A*X-B).print();

   Matrix V;
   Vector LL;
   Matrix S = A.t() * A;

   S.factorEV(LL,V);
   LL.print();

   for(r=0;r<8;r++) {
      Matrix lI(8,8);
      lI = LL(r);

      double d = (S-lI).det();
      printf("%lg\n\n",d);

      Vector x = V.getColumn(r);
      Vector z = S*x-LL(r)*x;
      z.print();
   }

   Matrix UU,VVt;
   Vector SS;
   A.factorSVD(UU,VVt,SS);

   for(r=0;r<8;r++) {
      SS.cell(r) *= SS(r);
   }

   (SS-LL).print();

   ((UU.t() * UU) + (VVt.t() * VVt)).print();
}
*/

/*
__declspec(naked) 
void simdDotProduct4(float *x,float *y,double *result) {
    __asm {
        movaps  xmm0, dword ptr[esp+4]      ;
        mulps       xmm0, dword ptr[esp+8]      ;
        movaps  xmm1, xmm0                      ;
        shufps  xmm1, xmm1, 4Eh             ;
        addps       xmm0, xmm1                      ;
        movaps  xmm1, xmm0                      ;
        shufps  xmm1, xmm1, 11h             ;
        addps       xmm0, xmm1
        movss       dword ptr[esp+0Ch], xmm0    ;
    }
}
*/

