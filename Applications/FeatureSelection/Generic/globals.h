#ifndef _GLOBALS_H_
#define _GLOBALS_H_

#include "matrix.h"

// Insert the element elt diagonally into matrix A starting at row r, col c
void insertDiagonal(pauly::Matrix &A,int r, int c, int n,double elt);

// Insert the element elt vertically into matrix A at row r and column c
void insertVertical(pauly::Matrix &A,int r, int c, int n,double elt);

// Generate a vector of random elements in range 0..mlt
pauly::Vector randVector(int n,double mlt=1.0);

#endif //_GLOBALS_H_
