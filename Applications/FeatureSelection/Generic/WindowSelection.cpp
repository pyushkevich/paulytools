#include "WindowSelection.h"
#include "globals.h"

#include <iostream>
#include <list>
using std::cout;
using std::endl;
using std::list;

WindowSelection::WindowSelection() 
: AbstractLPClassification()
{
	alpha = 5.0;
	eps = 1e-15;
	m_provider = NULL;
}

WindowSelection::~WindowSelection() 
{
}

void WindowSelection::initialize(const Mat &A, const Mat &B, const Mat &Omega) 
{
	assert(m_provider);

	cout << endl << "-------------- Starting Multi-Window Search ------------" << endl;

	// Inverse of the O matrix
	Mat Ot = Omega.transpose();

	// Record the dimensions of various matrices
	m = A.rows();
	k = B.rows();
	n = A.columns();
	nw = Omega.columns();

	// Construct the linear programming matrix
	Mat M(m+k+n,m+k+n+nw+1);

	// Insert the data matrices and window matrix
	M.update(A,0,0);
	M.update(-A*Omega,0,n+m+k);
	M.update(-B,m,0);
	M.update(B*Omega,m,n+m+k);
	M.update(Omega * 2,m+k,n+m+k);

	// Insert all the identities
	insertDiagonal(M,0     ,n ,m+k ,1);
	insertDiagonal(M,m+k   ,0 ,n   ,-1);

	// Insert the vertical unit vectors
	insertVertical(M,0,n+m+k+nw,m,-1);
	insertVertical(M,m,n+m+k+nw,k,1);

	// Now, specify the vector B
	Vec b(m+k+n);
	fillRange(b,0,m+k,1);

	// Construct the basic vector C
	c.set_size(n+m+k+nw+1);
	fillRange(c,n,m,1.0/m);
	fillRange(c,n+m,k,1.0/k);

	// Construct the upper and lower limits
	Vec upper(n+m+k+nw+1), lower(n+m+k+nw+1);
	upper.fill(1e100);
	lower(n+m+k+nw) = -1e100;

	// Create a soplex problem
	m_provider->setProblem(c,M,b,lower,upper);

	// Compute window lengths
	wl.set_size(nw);
	for(int wc=0;wc<nw;wc++) {
		wl(wc) = 0;
		for(int wr=0;wr<n;wr++) {
			wl(wc) += Omega(wr,wc);
		}
	}

	// Set the weight vector
	wgt.set_size(nw);
	wgt.fill(0);

	state = INIT;
}

void WindowSelection::setLinearModulation(double lambda, double eta) {
	assert(state != UNINIT);

	cout << "-------------- New Parameter Value" << endl;
	cout << "lam = " << lambda << endl << "eta = " << eta << endl;

	wgt.fill(eta);
	wgt += lambda*wl;

	state = INIT;
}

void WindowSelection::setArbitraryModulation(const Vec &weight) {
	assert(weight.size()==nw);
	assert(state != UNINIT);

	wgt = weight;
	state = INIT;
}

bool WindowSelection::run() {
	assert(state != UNINIT && m_provider);

	// This list keeps the sequence of results r = [w y z v g]
	typedef list<Vec> VecListType;
	VecListType r;
	r.push_back(randVector(n+m+k+nw+1));

	// This is the status
	bool status;
	double fLast;
	do {

		// Update the objective vector c
		int i,j;
		for(i=0;i<nw;i++) {
			j=n+m+k+i;
			c(j) = wgt(i) * alpha * exp(-alpha * r.back()(j));
		}

		// Compute the pre-optimization objective value
		fLast = dot_product(c,r.back());

		// Perform the linear programming task!
		m_provider->updateObjective(c);
		r.push_back(Vec(n+m+k+nw+1));

		// Get status
		status = m_provider->solve(r.back());

		// Now, continue if the objective has changed
	} while(status && fabs(fLast - dot_product(c,r.back())) > eps);

	// Store the results
	Vec rf = r.back();
	Vec ww = rf.extract(n,0);
	Vec yy = rf.extract(m,n);
	Vec zz = rf.extract(k,n+m);
	Vec vv = rf.extract(nw,n+m+k);

	// Compute the objective function
	objResult = yy.one_norm() / m + zz.one_norm() / k;
	for(int i=0;i<nw;i++) {
		objResult += wgt(i) * (1 - exp(-alpha * vv(i)));
	}
	vResult = vv;
	wResult = ww;
	gResult = rf(m+n+k+nw);

	// Record the status
	state = status ? SUCCESS : FAILURE;
	return status;
}
