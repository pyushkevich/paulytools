/******************************************************************
 * MATRIX Library                                                 *
 ******************************************************************
 * Author:                 Paul Yushkevich
 *
 * Date:                   Apr 1, 1999
 *
 * Description             Basic and numerical matrix operations   
 *	                        See http://www.cs.unc.edu/~pauly/matrix
 *									
 *	Sources:                Uses my compilation of CLAPACK
 *									
 * Dependencies:           CLAPACK, INTEL MKL (optional)
 *
 * Flags							_INTEL_SIMD_				use Intel extensions
 *									_MATRIX_BOUNDS_CHECK		use parameter bounds checking
 ******************************************************************
 * matrix.h
 *	---------
 * Declarations for matrix and vector classes
 ******************************************************************/
#ifndef MATRIX_H
#define MATRIX_H

#include <ostream>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include <mylibs.h>

// Begin namespace
NAMESPACE_PAULY_START

class Vector;
class Matrix;

/**
 * Debug macros
 */
#ifdef _MATRIX_BOUNDS_CHECK_BALKAKA
#define dassert(a) assert(a);
#else
#define dassert(a) ;
#endif

// Types of possible matrices
#define GENERAL      0
#define TRIANGULAR   -1
#define SPD          -2

/*

 **********************************************************************
 * If you decide to use SIMD extensions, then use the following #define
 * and link with MLK libraries.  Otherwise, link with clapack.lib
 **********************************************************************

*/
#ifdef _INTEL_SIMD_	
   typedef int lapack_int;
#else
   typedef long lapack_int;
#endif


/**
 * Matrix of double
 */
class Matrix {
protected:
   // Number of rows and columns in the matrix
   int nRows,nCols,nCells;
   
   // The data array
   double *data;
   
   // The columns index into the data - for fast access
   double **cols;

   // A couple init methods
   void initStorage(int nRows,int nCols);
   void initData(double *data = NULL,bool columnMajorOrder = true);

   int ipFactorGenericLU(lapack_int pivot[]);
   int ipFactorSPDLU();

public:
	// This constructor initializes the matrix with no data
   // Do not attempt matrix operations with such matrices - they wll
   // crash (assertions will fail).
   Matrix() : nRows(0),nCols(0),nCells(0),data(NULL),cols(NULL) {
   }
   
   // This is the right way to initialize the matrix.  Data contains values
   // in column major order, unless the flag is changed
   // This is also a default constructor
   Matrix(int rows,int columns,double *data=NULL,bool columnMajorOrder=true) {
      initStorage(rows,columns);
      initData(data,columnMajorOrder);
   }

   // Creates a matrix and reads in values in column major order
   // Make sure that constants are specified correctly, e.g. MatrixTmp<double> (1,2,3,3) will 
   // cause problems but MatrixTmp<double> (1,2,3.0,3.0) will not
   Matrix(int rows,int columns,double firstValue,...);
   

	// Copy constructor
   Matrix(const Matrix &m) {
      initStorage(m.nRows,m.nCols);
      initData(m.data);
   }

   // Destructor
   virtual ~Matrix() {
      if(data) {
         delete[] data;
         delete[] cols;
      }
   }

	// This method returns the data pointer (in column major order).  This is useful for GL and 
	// for interfacing with other packages.
	double *getDataArray() const{
		return data;
	}

   // Read only member access - inline for extra speed
	const double operator() (int row, int col) const {
      dassert(data && cols);
		dassert(row >= 0 && row < nRows);
      dassert(col >= 0 && col < nCols);
		return cols[col][row];
	}

   // Read-Write member access - inline for extra speed
   double& operator() (int row,int col) {
      dassert(data && cols);
		dassert(row >= 0 && row < nRows);
      dassert(col >= 0 && col < nCols);
      return cols[col][row];
   }

	// Set all elements of a matrix to a number; different from M=I
	void setAll(double k) {
     for(int i=0;i < nCells;i++)
         data[i] = k;
   }

	// Resizes matrix
   void setSize(int rows,int columns) {
      if(rows != nRows || columns != nCols) {
         if(data) {
            delete[] data;
            delete[] cols;
         }
         initStorage(rows,columns);
      }
      initData(NULL);      
   }

   // Matrix addition
   virtual void operator += (const Matrix &A) {
      dassert(nRows == A.nRows && nCols == A.nCols);
      for(int i=0;i < nCells;i++)
         data[i] += A.data[i];
   }

   virtual Matrix operator +(const Matrix &A) const { 
      dassert(nRows == A.nRows && nCols == A.nCols);
      Matrix C;
      C.initStorage(nRows,nCols);

      for(int i=0;i < nCells;i++)
         C.data[i] = data[i] + A.data[i];

      return C;
   }

   // Matrix substraction
	virtual void operator -= (const Matrix &A) {
      dassert(nRows == A.nRows && nCols == A.nCols);
      for(int i=0;i < nCells;i++)
         data[i] -= A.data[i];
   }   

   virtual Matrix operator -(const Matrix &A) const { 
      dassert(nRows == A.nRows && nCols == A.nCols);
      Matrix C;
      C.initStorage(nRows,nCols);

      for(int i=0;i < nCells;i++)
         C.data[i] = data[i] - A.data[i];
      
      return C;
   }

   // Matrix multiplication
   virtual void operator *= (const Matrix &A) {
      *this = *this * A;
   }

   virtual Matrix operator * (const Matrix &A) const;

   // Multiplication by constant
   virtual void operator *= (const double k) {
      for(int i=0;i < nCells;i++)
         data[i] *= k;
   }

   virtual Matrix operator *(const double k) const { 
      Matrix C;
      C.initStorage(nRows,nCols);

      for(int i=0;i < nCells;i++)
         C.data[i] = data[i] * k;
      
      return C;
   }

   virtual void operator /= (const double k) {
      for(int i=0;i < nCells;i++)
         data[i] /= k;
   }

   virtual Matrix operator /(const double k) const { 
      Matrix C;
      C.initStorage(nRows,nCols);

      for(int i=0;i < nCells;i++)
         C.data[i] = data[i] / k;
      
      return C;
   }

   // Negation
   virtual void negate() {
      for(int i=0;i < nCells;i++)
         data[i] = -data[i];
   }

   virtual Matrix operator -() const { 
      Matrix C;
      C.initStorage(nRows,nCols);

      for(int i=0;i < nCells;i++)
         C.data[i] = -data[i];
      
      return C;
   }

	// This assignment operator assign the matrix k*I to a square matrix
   Matrix& operator= (const double k) {
      dassert(nRows==nCols);
      initData(NULL);
      for(int i=0;i<nRows;i++)
         cols[i][i] = k;
      return *this;
   }

   // Copy operator
   Matrix& operator= (const Matrix &A) {
      if(nRows != A.nRows || nCols != A.nCols) {
         if(data) {
            delete[] data;
            delete[] cols;
         }
         initStorage(A.nRows,A.nCols);
      }
      initData(A.data);
      return *this;
   }

	// Transpose operator - returns the transpose
	virtual Matrix t() const;

   // Compute the determinant
	double det(int type=GENERAL);
	double logDeterminant(double &outSign,int type=GENERAL);

   // Solve A*x = B system
   // int solveGE(const Matrix &A,Matrix &C,int type=GENERAL);

	/**
	 * Solve System of Linear Equations (In Place)
	 *
	 * It solves the linear system AX = B, places the result in B and the 
	 * LU decomposition of A in A.  THIS METHOD IS DESTRUCTUVE
    */
	int ipSolveLinearSystem(Matrix &B);

   // LU Factorization
   int factorLU(Matrix &L, Matrix &U, Matrix &P,int type=GENERAL);

   // SVD Factorization
   int factorSVD(Matrix &U,Matrix &Vt,Vector &sv,int type=GENERAL);

   // Eigenvalue analysis
   int factorEV(Vector &l,Matrix &V,int type=GENERAL);

   // Get inverse
   int ipInverse();

	// Compute the mean and covariance of a feature matrix
	Vector computeMean() const;
	Vector computeMean(Matrix &woMean) const;
	void computeMeanCov(Vector &mean,Matrix &cov) const ;
	void computeMeanCorr(Vector &mean,Vector &sigma,Matrix &corr) const ;

	// Compute the PCA of this matrix, where each row is a data sample.
	int computePCA(Vector &mean,Vector &eVal,Matrix &eVec) const ;
	int computeCorrelationPCA(Vector &mean,Vector &sdev,Vector &eVal,Matrix &eVec) const ;

   // get/extract submatrices
	virtual void insertMatrix(int row,int col,const Matrix &m);
	virtual void extractMatrix(int row,int col,Matrix &m) const;
   virtual Vector getColumn(int col) const;
   virtual Vector getRow(int row) const;

	// Apply a double to double routine to every element in the matrix
	// e.g. m.applyToElements(sqrt)
	void applyToElements(double (*func)(double)) {
		for(int i=0;i<nCols;i++)
			for(int j=0;j<nRows;j++) 
				cols[i][j] = func(cols[i][j]);
	}

	int rows() const {
      return nRows;
   }

	int columns() const {
      return nCols;
   }

   // Norms of the matrix
   double oneNorm() const ;
   double infinityNorm() const;
   double pNorm(double p) const;
	double twoNorm() const;

	// Swap two rows/columns
   virtual void swapRows(int r1,int r2);
	virtual void swapColumns(int c1,int c2);
	
	// Reverse rows in a range (-1 here means all following rows
	void reverseRows(int startRow,int endRow);
	void reverseColumns(int startColumn,int endColumn);

	// Flip the order of all rows in a matrix
	void reverseAllRows() {
		reverseRows(0,nRows-1);
	}

	// Flip the order of all columns in a matrix
	void reverseAllColumns() {
		reverseColumns(0,nCols-1);
	}

   // Output operators
   void print();

	// This prints the matrix out in a matlab fashion
	void printMatlab(FILE *f = stdout) const;
	void printMatlab(const char *fname) const;


	static void multiply(const Matrix &A, const Matrix &B, Matrix &C);
	static void multiply(const Matrix &A, const Vector &X, Vector &R);
};

// Backwards multiplication operator
inline Matrix operator * (const double k,const Matrix &A) {
   return A*k;
}


/**
 * Generic vector class - a Nx1 matrix
 */
class Vector : public Matrix
{
protected:
public:
   Vector() : Matrix() {};

   // Creates a vector, feeds it some data
   Vector(int dim,double *data=NULL) : Matrix(dim,1,data) {};

   // Creates a vector, may be initialized with doubles (!!!)
   Vector(int dim,double firstValue,...);

	Vector(const Matrix &m) : Matrix(m) {
		dassert(m.columns() == 1);
	}

	// Read only access operator
	double operator() (int row) const {
		return data[row];
	}

   // Read-write access operator
   double &operator() (int row) {
      return data[row];
   }

   // Set size
   void setSize(int rows) {
      Matrix::setSize(rows,1);
   }

	// Cross product
	virtual Vector cross (const Vector& v) {
		dassert(nRows == v.nRows && nRows >= 3);
		Vector result(nRows);

		result.data[0] = data[1]*v.data[2] - data[2]*v.data[1];
		result.data[1] = data[2]*v.data[0] - data[0]*v.data[2];
		result.data[2] = data[0]*v.data[1] - data[1]*v.data[0];

		return result;
	}

	// Dot product
	double dotProduct(const Vector &v) const {
      dassert(nRows == v.nRows);

      double dp = 0.0;
		for(int i=0;i<nRows;i++)
			dp += data[i] * v.data[i];
		return dp;
	}

   // Normalizes the vector
   virtual void normalize() {
		double len = twoNorm();
		if(len != 0)
			(*this) /= len;
	}

	// Return number of rows
   virtual int size() const {
	  return rows();
	}
};

Matrix *allocMatrixArray(int nMatrices,int rows,int columns);
Vector *allocVectorArray(int nVectors,int size);

double unif_rand_dbl(long *idum);
int randint(int thismax, long *seed);

// End namespace
NAMESPACE_PAULY_END

#endif
