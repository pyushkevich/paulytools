#include "matrix.h"
#include "FeatureSelection.h"
#include "SoPlexWrapper.h"
#include "CrossValidation.h"

#include "mathlink.h"

#include <ctime>
#include <iostream>
#include <vector>
#include <map>
#include <string>
using namespace std;

// ---- Global variables

struct FSInit {
	Mat A,B,O;
	double alpha;
};

// Feature selection drivers
map <int, WindowSelection *> wsDriver;
map <int, FSInit *> wsInit;

// Shared provider
SoPlexWrapper *sharedProvider = NULL, *fsProvider = NULL;


// Handle index
int handle = 0;

void readMMAMatrix(Mat &M) {
	long nRows,nCols;
	double **m;

	// Load the data into matrices
	MLCheckFunction(stdlink, "List", &nRows);
	m = new double *[nRows]; 
	for(int i=0;i<nRows;i++) {
		MLGetRealList(stdlink, &(m[i]), &nCols);
	}

	// Set the matrix values
	M.setSize(nRows,nCols);
	for(int r=0;r<nRows;r++)
		for(int c=0;c<nCols;c++)
			M(r,c) = m[r][c];


	// Clean up the data
	for(i=0;i<nRows;i++) {
		MLDisownRealList(stdlink, m[i], nCols);
	}
	delete m;
}

void readMMAVector(Vec &v) {
	long n;

	// Read the vector
	MLCheckFunction(stdlink,"List",&n);
	v.setSize(n);

	// Set the matrix values
	for(int r=0;r<n;r++) {
		double x;
		MLGetReal(stdlink,&x);
		v(r) = x;
	}
}

/**
 * This is a class used to hold the results of a feature search
 */
struct WinSearchResult {
  double obj,gamma;
  Vec v,w;
};


int lpDeleteFeatureSelector(int) {
	try {

		if(wsDriver[handle]) {
			delete wsDriver[handle];
			wsDriver[handle] = NULL;

			delete wsInit[handle];
			wsInit[handle] = NULL;
		}
	
		return handle;

	} catch(...) {
		return -1;
	}
}

int lpNewFeatureSelector() {
	try {
		FSInit *initData = new FSInit;

		readMMAMatrix(initData->A);
		readMMAMatrix(initData->B);
		readMMAMatrix(initData->O);
		MLGetReal(stdlink,&initData->alpha);

		// SoPlexWrapper *wrapper = new SoPlexWrapper();
		if(!fsProvider)
			fsProvider = new SoPlexWrapper();

		WindowSelection *ws = new WindowSelection();
		ws->setLPProvider(fsProvider);
		ws->initialize(initData->A,initData->B,initData->O);
		ws->setAlpha(initData->alpha);

		handle++;
		wsDriver[handle] = ws;
		wsInit[handle] = initData;

		return handle;	

	} catch(...) {
		cout << "Exception in NewFeatureSelector" <<  endl;
		return -1;
	}	
}

void lpRunFeatureSelection(int handle, double lambda, double eta, int nRuns) {
	int i;

	if(!wsDriver[handle]) {
		// Return error
		MLPutFunction(stdlink,"List",0);
		return;
	}

	WindowSelection *ws = wsDriver[handle];

	// Keep the best attempt
	WinSearchResult wsr;
	bool readyToCompare = false;

	// Perform nRuns attempts
	for(i=0;i<nRuns;i++) {	
		
		// Try to run the next experiment
		try {			
			ws->setLinearModulation(lambda,eta);
			if(ws->run()) {					
				if(!readyToCompare || ws->getObjectiveValue() < wsr.obj) {
					wsr.obj = ws->getObjectiveValue();
					wsr.gamma = ws->getSeparatingPlane(wsr.w);
					ws->getWindowVector(wsr.v);
					readyToCompare = true;
				}
			} else {
				cerr << "Window selection failed for lambda=" << lambda << ", eta = " << eta << endl;
			}
		}
		catch(...) {
			cerr << "Exception caught for lambda = " << lambda << ", eta = " << eta << endl;			

			// Reinitialize problem
			try{
				delete wsDriver[handle];
				delete fsProvider;
			} catch(...) {}

			fsProvider = new SoPlexWrapper();
			wsDriver[handle] = new WindowSelection();
			wsDriver[handle]->setLPProvider(fsProvider);
			wsDriver[handle]->initialize(wsInit[handle]->A,wsInit[handle]->B,wsInit[handle]->O);
			wsDriver[handle]->setAlpha(wsInit[handle]->alpha);
		}
	}

	// Compute the mask vector
	Vec z = wsInit[handle]->O * wsr.v;
	int *idx = new int[z.size()];

	// Count the number of windows and the number of features
	int nw = 0, nf = 0;
	for(i=0;i<wsr.v.size();i++) {
		if(wsr.v(i) > 0) nw++;
	}
	for(i=0;i<z.size();i++) {
		if(z(i) > 0) {
			nf++;
			idx[i] = 1;
		}
		else 
			idx[i] = 0;
	}

	// Print out the results of the optimization
	cout << " ************************************************************ " << endl;
	cout << " Parameters          : " << lambda << ", " << eta << endl;
	cout << " Best objective      : " << wsr.obj << endl;
	cout << " Selected windows    : " << nw << endl;
	cout << " Selected features   : " << nf << endl;
	cout << " ************************************************************ " << endl;

	// Return the results
	MLPutFunction(stdlink,"List",4);
	MLPutReal(stdlink,wsr.obj);
	MLPutReal(stdlink,wsr.gamma);
	MLPutRealList(stdlink,wsr.v.getDataArray(),wsr.v.size());
	MLPutIntegerList(stdlink,idx,z.size());

	delete idx;
}	



void lpSeparation(void) {
	Mat A,B;

	// Read in the data matrices
	readMMAMatrix(A);
	readMMAMatrix(B);

	try {
		// Separation problem
		if(sharedProvider == NULL)
			sharedProvider = new SoPlexWrapper;

		LinearSeparation sp;
		sp.setLPProvider(sharedProvider);
		sp.initialize(A,B);

		if(sp.run()) {
			Vec w;
			double g = sp.getSeparatingPlane(w);

			MLPutFunction(stdlink,"List",3);
			MLPutReal(stdlink,sp.getObjectiveValue());
			MLPutReal(stdlink,g);
			MLPutRealList(stdlink,w.getDataArray(),w.rows());
			return;
		}
	} catch(...) {
		cerr << "Exception caught in lpSeparation" << endl;
		try{
			delete sharedProvider;			
		} catch(...) {
		}
		sharedProvider = NULL;
	}

	MLPutFunction(stdlink,"List",0);
}

int main(int argc, char *argv[])
{
	return MLMain(argc, argv);
}