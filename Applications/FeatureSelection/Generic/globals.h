#ifndef _FSELECT_GLOBALS_H_
#define _FSELECT_GLOBALS_H_

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
                                                                                                                                                                        
// Shorthand for matrix and vector classes that doesn't clash with other code
typedef vnl_matrix<double> Mat;
typedef vnl_vector<double> Vec;

// Insert the element elt diagonally into matrix A starting at row r, col c
void insertDiagonal(Mat &A,int r, int c, int n,double elt);

// Insert the element elt vertically into matrix A at row r and column c
void insertVertical(Mat &A,int r, int c, int n,double elt);

// Insert a range of elements into a vector
void fillRange(Vec &A, int r, int n, double elt);

// Generate a vector of random elements in range 0..mlt
Vec randVector(int n,double mlt=1.0);

#endif //_GLOBALS_H_
