/**
This code runs Mathematica experiments
*/

#include "matrix.h"
#include "mathematicaio.h"
#include "FeatureSelection.h"
#include "SoPlexWrapper.h"
#include "CrossValidation.h"

#include <ctime>
#include <iostream>
#include <vector>
#include <map>
#include <string>
using namespace std;

// A global LP provider, initialized in main
LPWrapper *lpProvider = NULL;
LPWrapper *lpProviderCVX = NULL;

// A global option for alpha parameter
double alpha = 5.0;


void usage(void) {
	cout << "Bad input. Please read manual.html for usage instructions." << endl;
}

/**
 * This is a class used to hold the results of a feature search
 */
struct WinSearchResult {
  double obj,gamma;
  Vec v,w;
  double xv;
};

/**
 * Make a bindary mask from such a result
 */
string makeMask(const WinSearchResult &R) {

	// Construct a textual mask
	string mask;
	for(int i=0;i<R.v.size();i++) {
		// mask.push_back(R.v(i)==0 ? '0' : '1');
        mask += (R.v(i)==0 ? '0' : '1');
	}
	return mask;
}

void findFeatureWindow(const Mat &A, const Mat &B,const Mat &O, const Vec &lam, 
		               const Vec &eta, vector<WinSearchResult> &out)
{
	// Create a feature selection object
	WindowSelection ws;
	ws.setLPProvider(lpProvider);
	ws.setAlpha(alpha);
	ws.initialize(A,B,O);

	// Run the window selection experiments
	for(int i=0;i<lam.size();i++) {
		ws.setLinearModulation(lam(i),eta(i));
		if(ws.run()) {
			WinSearchResult wsr;
			
			wsr.obj = ws.getObjectiveValue();
			wsr.gamma = ws.getSeparatingPlane(wsr.w);
			ws.getWindowVector(wsr.v);

			out.push_back(wsr);
		} else {
			cerr << "Window selection failed for lambda=" << lam(i) << ", eta = " << eta(i) << endl;
		}
	}	
}
/*
class FSMainDriver {
public:
	void setProblem(const Mat &A, const Mat &B,const Mat &O, const Mat &p);
	void setXVal(int n, double plo);
	void setOutput(const char *file);
	void run();

private:
	WindowSelection *sel;
	Mat A,B,O,p,pLeft;

};
*/
void findFeaturesFSV(const Mat &A, const Mat &B, const Vec &L,vector<WinSearchResult> &out) 
{
	// Create a feature selection object
	FeatureSelection fs;
	fs.setLPProvider(lpProvider);
	fs.setAlpha(alpha);
	fs.initialize(A,B);
	
	// Run the feature selection experiments
	for(int i=0;i<L.size();i++) {
		fs.setLambda(L(i));
		if(fs.run()) {
			WinSearchResult wsr;
			
			wsr.obj = fs.getObjectiveValue();
			wsr.gamma = fs.getSeparatingPlane(wsr.w);
			fs.getWindowVector(wsr.v);

			out.push_back(wsr);
		} else {
			cerr << "Feature selection failed for lambda=" << L(i) << endl;
		}
	}
}

void padVector(const Vec &v,Vec &p,int reps) {
  p.set_size(v.size()*reps);
  for(int r=0;r<v.size();r++)
    for(int i=0;i<reps;i++)
      p(r*reps+i)=v(r);
}

void unpadResults(const vector<WinSearchResult> &res,int reps,vector<WinSearchResult> &out) {
  out.clear();
  for(int i=0;i<res.size()/reps;i++) {
    double spanMin = 1e300;
    int spanIdx = -1;
    for(int j=0;j<reps;j++) {
      if(spanMin > res[i*reps+j].obj) {
	spanMin = res[i*reps+j].obj;
	spanIdx = j;
      }
    }
    out.push_back(res[i*reps+spanIdx]);
  }
}

void processParameterSet(Mat &A, Mat &B, Mat &O, Mat&p, double alpha, int nRuns, FILE *fout, 
						 bool isXVOn = false, double pXVLeaveOut = 0.5, int nXVRuns = 100)
{
	int q,i;

	// This map stores the crossvalidation result for each window tested
	map<string,double> xvMap;

	// This flag indicates that a reinitialization is necessary
	bool reinit = false;

	// Create a new window selection
	WindowSelection *ws = new WindowSelection;
	ws->setLPProvider(lpProvider);
	ws->setAlpha(alpha);
	ws->initialize(A,B,O);	

	// Dump out the results and crossvalidate on the run
	fprintf(fout,"{");

	// Number of failures
	int failCount = 100;

	// Go through the iterations (values of parameter p)
	for(q=0;q<p.rows();q++) {
		
		// Do we need to reinitialize?
		if(reinit) {

			// Create a feature selection object
			lpProvider = new SoPlexWrapper();
			lpProviderCVX = new SoPlexWrapper();
			ws = new WindowSelection();
			ws->setLPProvider(lpProvider);
			ws->setAlpha(alpha);
			ws->initialize(A,B,O);

			reinit = false;
		}

		// Try to process the next parameter pair
		try {			
			// Keep the best attempt
			WinSearchResult wsr;

			// Perform nRuns attempts
			for(i=0;i<nRuns;i++) {				
				ws->setLinearModulation(p(q,0),p(q,1));
				if(ws->run()) {					
					if(i==0 || ws->getObjectiveValue() < wsr.obj) {
						wsr.obj = ws->getObjectiveValue();
						wsr.gamma = ws->getSeparatingPlane(wsr.w);
						ws->getWindowVector(wsr.v);
					}

					WinSearchResult wsr2;
					wsr2.obj = ws->getObjectiveValue();
					wsr2.gamma = ws->getSeparatingPlane(wsr2.w);
					ws->getWindowVector(wsr2.v);

					// Get the resulting feature mask
					Vec z = O * wsr2.v;

					// Count the number of windows and the number of features
					int nw = 0, nf = 0, j;
					for(j=0;j<wsr2.v.size();j++) 
						if(wsr2.v(j) > 0) nw++;
					for(j=0;j<z.size();j++) 
						if(z(j) > 0) nf++;

					// Print out the results of the optimization
					cout << " ............................................................ " << endl;
					cout << " Parameters          : " << p(q,0) << ", " << p(q,1) << endl;
					cout << " Objective           : " << wsr2.obj << endl;
					cout << " Selected windows    : " << nw << endl;
					cout << " Selected features   : " << nf << endl;
					cout << " ............................................................ " << endl;

				} else {
					cerr << "Window selection failed for lambda=" << p(q,0) << ", eta = " << p(q,1) << endl;
				}
			}

			// We now have the best attempt in wsr
			string mask = makeMask(wsr);

			// Get the resulting feature mask
			Vec z = O * wsr.v;

			// If nothing is associated with the map, do xvalidation
			wsr.xv = 0.0;
			if(isXVOn) {
				if(xvMap[mask]==xvMap["dummy"] ) {
					FeatureSetCrossValidation cv;
					cv.setLPProvider(lpProviderCVX);
					xvMap[mask] = cv.leaveSomeOut(A,B,z,pXVLeaveOut,nXVRuns);
				}
				wsr.xv = xvMap[mask];
			}

			// Count the number of windows and the number of features
			int nw = 0, nf = 0;
			for(i=0;i<wsr.v.size();i++) 
				if(wsr.v(i) > 0) nw++;
			for(i=0;i<z.size();i++) 
				if(z(i) > 0) nf++;

			// Print out the results of the optimization
			cout << " ************************************************************ " << endl;
			cout << " Parameters          : " << p(q,0) << ", " << p(q,1) << endl;
			cout << " Best objective      : " << wsr.obj << endl;
			cout << " X-validation        : " << wsr.xv << endl;
			cout << " Selected windows    : " << nw << endl;
			cout << " Selected features   : " << nf << endl;
			cout << " ************************************************************ " << endl;

			fprintf(fout,"{%lg,%lg,%lg,\"%s\",%lg}%c\n", p(q,0),p(q,1),wsr.obj,mask.c_str(),wsr.xv,q==p.rows()-1 ? ' ' : ',');
			fflush(fout);

		} catch (...) {
			cerr << "Exception caught !!!" << endl;
			reinit = true;

			if(failCount-- > 0)
				q--;
		}
	}	

	try {
		delete ws;
	} catch(...) {
		cerr << "Deletion exception caught !!!" << endl;
	}

	fprintf(fout,"}");
}

void leaveRowOut(const Mat &A, Mat &B, int row) {
	B.set_size(A.rows()-1,A.columns());

	int ia = 0;
	for(int j=0;j<A.rows();j++) {
		if(j==row) 
			continue;
        for(int k=0;k<A.columns();k++) {
			B(ia,k) = A(j,k);
		}
		ia++;
	}
}


int main(int argc, char *argv[]) {

	// Read in a text file that contains Mathematica code
	if(argc < 3) {
		usage();
		return -1;
	}

	// Some of the parameters
	char *inFile = NULL, *outFile = NULL;
	bool isXVOn = false;
	bool isHardXVOn = false;
	double pXVLeaveOut = 0.5;
	int nRuns = 10, nXVRuns = 100;

	// Read the options and check for errors
	try {
		for(int p=1;p<argc;p++) {
			if(0 == strcmp(argv[p],"-x"))
				isXVOn = true;
			if(0 == strcmp(argv[p],"-loo"))
				isHardXVOn = true;
			else if(0 == strcmp(argv[p],"-xp"))
				pXVLeaveOut = atof(argv[++p]);
			else if(0 == strcmp(argv[p],"-xn"))
				nXVRuns = atoi(argv[++p]);
			else if(0 == strcmp(argv[p],"-a"))
				alpha = atof(argv[++p]);
			else if(0 == strcmp(argv[p],"-r")) 
				nRuns = atoi(argv[++p]);
			else if(inFile == NULL)
				inFile = argv[p];
			else if(outFile == NULL)
				outFile = argv[p];
			else
				throw "Too many arguments";
		}
	} catch(...) {
		usage();
		throw;
	}

	// Print out options
	cout << "Using " << nRuns << " repetitions per parameter pair" << endl;
	if(isXVOn) {
		cout << "Using cross validation" << endl;
		cout << "Leave out probability: " << pXVLeaveOut << endl;
		cout << "Leave out experiments: " << nXVRuns << endl;
	}
	else {
		cout << "Not using cross validation" << endl;
	}


	// Open files for reading and writing
	FILE *f = fopen(inFile,"rt");
	FILE *fout = fopen(outFile,"wt");
	if(f == NULL || fout == NULL) {
		usage();
		return -1;
	}

	// Feature matrices A, B
	Mat A,B;

	// Window matrix Omega;
	Mat O;

	// Window selection parameter list
	Mat p;

	// Read the matrices and check validity
	try {
		MathematicaIO::readMathematicaMatrix(f,A);
		cout << "A is " << A.rows() << " by " << A.columns() << endl;  

		MathematicaIO::readMathematicaMatrix(f,B);
		cout << "B is " << B.rows() << " by " << B.columns() << endl;

		MathematicaIO::readMathematicaMatrix(f,O);
		cout << "O is " << O.rows() << " by " << O.columns() << endl;  

		MathematicaIO::readMathematicaMatrix(f,p);  
		cout << "p is " << p.rows() << " by " << p.columns() << endl;  

		if(A.columns() != B.columns() || A.columns() != O.rows())
			throw "The dimensions of the input matrices A, B, O are incompatible";
		if(p.columns() != 2)
			throw "The parameter matrix p must have 2 columns";
	} catch(...) {
		usage();
		throw;
	}

	// Randomize
	srand(clock());

	// Rewrite the vectors to do multiple repetitions
	// Vec lambdaPad,etaPad;
	// padVector(p.getColumn(0),lambdaPad,nRuns);  
	// padVector(p.getColumn(1),etaPad,nRuns);

	// Initialize the LP provider (should use factory here)
	
	// TEST
	// LinearSeparation ls;
	// ls.setLPProvider(lpProvider);
	// ls.initialize(A,B);
	// ls.run();

	// Perform the tests
	// vector<WinSearchResult> results,bestResults;
	// findFeatureWindow(A,B,O,lambdaPad,etaPad,results);
	// unpadResults(results,nRuns,bestResults);

	// Create providers
	lpProvider = new SoPlexWrapper();
	lpProviderCVX = new SoPlexWrapper();

	if(!isHardXVOn) {
		processParameterSet(A,B,O,p,alpha,nRuns,fout,isXVOn,pXVLeaveOut,nXVRuns);
	}
	else {
		Mat L;
		int i;

		fprintf(fout,"{\n");	
		for(i=0;i<A.rows();i++) {
			fprintf(fout,"{1,%d,\n",i);			
			leaveRowOut(A,L,i);
			processParameterSet(L,B,O,p,alpha,nRuns,fout,isXVOn,pXVLeaveOut,nXVRuns);
			fprintf(fout,"},\n");			
		}
		for(i=0;i<B.rows();i++) {
			fprintf(fout,"{2,%d,\n",i);			
			leaveRowOut(B,L,i);
			processParameterSet(A,L,O,p,alpha,nRuns,fout,isXVOn,pXVLeaveOut,nXVRuns);
			fprintf(fout,"}%s\n",(i==B.rows()-1 ? "}": ","));			
		}
	}

	fclose(fout);
}

