#include <iostream>
#include <cstdlib>
#include "SparseMatrix.h"

using namespace std;

int main (int argc, char *argv[]) {
	if (argc < 5) {
		cerr << "Usage: " << argv[0] << " dim row column value" << endl;
		exit(1);
	}
	
	int dim = atoi(argv[1]);
	int row = atoi(argv[2]);
	int col = atoi(argv[3]);
	double value = atof(argv[4]);
	SparseMatrix A(dim);
	A.updateElement(row, col, value);
	cout << A << endl;
	
	return 0;
}

