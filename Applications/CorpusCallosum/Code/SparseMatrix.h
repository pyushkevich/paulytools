
#ifndef _SparseMatrix_H
#define _SparseMatrix_H

#include <iostream>

using namespace std;

class SparseMatrix {
	public:
	const int dim;
	// the number of none zero elements before this row
	int *rowIndex;
	// the column index of none zero elements in each row
	int *columnIndex;
	double *dataArray;
		
	public:
	SparseMatrix (const int dim);
	~SparseMatrix ();
	
	// row: matrix's row index
	// column: the index for the three non-zero element, 0 through 2
	void SparseMatrix::updateElement (const int row, const int column, const double dataValue);	
	
	friend ostream& operator<< (ostream& out, const SparseMatrix& sm);
};

#endif

