#include "globals.h"

typedef pauly::Matrix Mat;
typedef pauly::Vector Vec;

void insertDiagonal(Mat &A,int r, int c, int n,double elt) {
	for(int i=0;i<n;i++)
		A(r+i,c+i) = elt;
}

void insertVertical(Mat &A,int r, int c, int n,double elt) {
	for(int i=r;i<r+n;i++)
		A(i,c) = elt;
}

Vec randVector(int n,double mlt) {
	Vec v(n);
	for(int i=0;i<n;i++)
		v(i) = (mlt * rand()) / RAND_MAX;
	return v;
}
