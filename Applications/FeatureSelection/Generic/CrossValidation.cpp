#include "CrossValidation.h"
#include <vector>
using namespace std;

FeatureSetCrossValidation::FeatureSetCrossValidation() {
	lpWrapper = NULL;
}

FeatureSetCrossValidation::~FeatureSetCrossValidation() {

}

void FeatureSetCrossValidation::setLPProvider(LPWrapper *provider) {
	this->lpWrapper = provider;
}

double FeatureSetCrossValidation::leaveSomeOut(const Mat &A, const Mat &B, const Vec &w, double pLeftOut, int N) 
{
	// Ret the numbers of rows and stuff
	int m,n,k; 
	m = A.rows();
	k = B.rows();
	n = A.columns();

	// Select the elements rows of A and B
	vector<int> fsel;
	unsigned int j;
	for(j=0;j<w.rows();j++) {
		if(w(j)!=0)
			fsel.push_back(j);
	}

	// Create the matrices X, Y
	Mat X(A.rows(),fsel.size()),Y(B.rows(),fsel.size());
	for(j=0;j<fsel.size();j++) {
		X.insertMatrix(0,j,A.getColumn(fsel[j]));   
		Y.insertMatrix(0,j,B.getColumn(fsel[j]));
	}

	int cl,el;

	// Random threshold
	int randThresh = (int)floor(RAND_MAX * pLeftOut);

	Mat Dtrain[2],Dtest[2];

	// Keep track of attempts and errors
	int att = 0;
	int err = 0;

	// Loop N times
	for(int i=0;i<N;i++) {

		// Select rows from X and Y
		for(cl=0;cl<2;cl++) {

			const Mat &D = cl ? Y : X;

			// Train set and test set
			vector<int> train,test;

			int el;
			for(el=0;el<D.rows();el++) {

				// Generate a random number
				int p = rand();
				if(p > randThresh) {
					train.push_back(el);
				} else {
					test.push_back(el);
				}
			}

			// Build matrices
			Dtrain[cl].setSize(train.size(),fsel.size());
			Dtest[cl].setSize(test.size(),fsel.size());

			for(el=0;el<Dtrain[cl].rows();el++)
				for(int col=0;col<D.columns();col++) 
					Dtrain[cl](el,col) = D(train[el],col);
			for(el=0;el<Dtest[cl].rows();el++) 
				for(int col=0;col<D.columns();col++) 
					Dtest[cl](el,col) = D(test[el],col);
		}

		// Perform classification
		LinearSeparation ls;
		ls.setLPProvider(lpWrapper);
		ls.initialize(Dtrain[0],Dtrain[1]);

		// Run LS
		ls.run();
		
		// Get the results
		Vec margin;
		double gamma = ls.getSeparatingPlane(margin);
		// margin.t().print();

		// Cross-validate
		for(cl=0;cl<2;cl++) {
			for(el=0;el<Dtest[cl].rows();el++) {
				Vec x = Dtest[cl].getRow(el);
				double score = x.dotProduct(margin)-gamma;
				int cls = score < 0 ? 1 : 0;
				att++;
				err += cls==cl ? 0 : 1;
			}
		}
	}

	return (1.0 * err) / att;
}

