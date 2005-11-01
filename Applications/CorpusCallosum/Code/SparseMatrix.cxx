#include "SparseMatrix.h"

SparseMatrix::SparseMatrix (const int _dim) : dim(_dim) {
	rowIndex = new int[(dim + 3) + 1];
	columnIndex = new int[3*(dim + 3)] ;
	dataArray = new double[3*(dim + 3)];
	
	// add 1
	rowIndex[0] = 1;
	for (int i = 0; i < dim + 3; ++i) {
		rowIndex[i+1] = 3*i + 4;
	}
	
	// the first and last row
	// fortran indexing convention applies
	columnIndex[0] = 1;
	columnIndex[1] = 2;
	columnIndex[2] = 3;
	columnIndex[3*(dim + 3) - 3] = dim + 1;
	columnIndex[3*(dim + 3) - 2] = dim + 2;
	columnIndex[3*(dim + 3) - 1] = dim + 3;
	
	// the rows in between
	for (int i = 1; i < dim + 2; ++i){
		columnIndex[3*i] = i ;
		columnIndex[3*i + 1] = i + 1;
		columnIndex[3*i + 2] = i + 2;
	}
}

SparseMatrix::~SparseMatrix () {
	delete[] rowIndex;
	delete[] columnIndex;
	delete[] dataArray;
}

void SparseMatrix::updateElement (const int row, const int column, const double dataValue) {
	dataArray[3*row + column] = dataValue;
}

ostream& operator<< (ostream& out, const SparseMatrix& sm) {
	out << "rowIndex = " << endl;
	for (int i = 0; i < sm.dim + 4; ++i) {
		out << sm.rowIndex[i] << endl;
	}
	
	out << "columnIndex = " << endl;
	for (int i = 0; i < sm.dim + 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			out << sm.columnIndex[3*i + j] << '\t';
		}
		out << endl;
	}
	
	out << "dataArray = " << endl;
	for (int i = 0; i < sm.dim + 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			out << sm.dataArray[3*i + j] << '\t';
		}
		out << endl;
	}
	out << flush;
}

