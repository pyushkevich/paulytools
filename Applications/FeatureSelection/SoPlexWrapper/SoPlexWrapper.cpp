// Include my wrapper
#include "SoPlexWrapper.h"

// Include the SoPlex library
#include "timer.h"
#include "spxpricer.h"
#include "spxdefaultpr.h"
#include "spxparmultpr.h"
#include "spxdevexpr.h"
#include "spxhybridpr.h"
#include "spxsteeppr.h"
#include "spxweightpr.h"
#include "spxratiotester.h"
#include "spxharrisrt.h"
#include "spxdefaultrt.h"
#include "spxfastrt.h"
#include "spxsimplifier.h"
#include "spxaggregatesm.h"
#include "spxredundantsm.h"
#include "spxrem1sm.h"
#include "spxgeneralsm.h"
#include "spxscaler.h"
#include "spxequilisc.h"
#include "spxsumst.h"
#include "spxweightst.h"
#include "spxvectorst.h"
#include "slufactor.h"

// Use the SoPlex Namespace
using namespace soplex;

bool SoPlexWrapper::solve(Vec &r) {
    
    Timer timer;
    SoPlex::Status status;

    // Solve the problem
    timer.start();
    status = plex->solve();
    timer.stop();

    // Print status
    switch (status)
    {
	case SoPlex::OPTIMAL: 
	    {
		// Return the optimal solution solution
		DVector rtn(r.size());
		plex->getPrimalUnscaled(rtn);
		for(int i=0;i<r.size();i++)
		    r(i) = rtn.get_ptr()[i];

		return true;
	    }

	case SoPlex::UNBOUNDED:
	    std::cout << "LP is unbounded";
	    break;
	case SoPlex::INFEASIBLE:
	    std::cout << "LP is infeasible";
	    break;
	case SoPlex::ABORT_TIME:
	    std::cout << "aborted due to time limit";
	    break;
	case SoPlex::ABORT_ITER:
	    std::cout << "aborted due to iteration limit";
	    break;
	case SoPlex::ABORT_VALUE:
	    std::cout << "aborted due to objective value limit";
	    break;
	default:
	    std::cout << "An error occurred during the solution process";

		// If the status is not optimal, the SoPlex should be 
		// recreated and re-initialized.
		setProblem(c,M,b,l,u);

	    break;
    }
    return false;
}

inline soplex::Vector vecCast(const Vec &v) {
	Real *r = new Real[v.size()];
	for(int i=0;i<v.size();i++)
		r[i] = v(i);
	soplex::Vector w(v.size(),r);
	// delete r;
	return w;
}

void SoPlexWrapper::updateObjective(const Vec &c) {
    this->c = c;

    Vec obj = -c;
    plex->changeObj(vecCast(obj));
		// soplex::Vector(obj.rows(),obj.getDataArray()));
}

void SoPlexWrapper::updateBounds(const Vec &u,const Vec &l) {
    this->u = u;
    this->l = l;
    
    Vec upper = u;
    Vec lower = l;
    plex->changeBounds(vecCast(lower),vecCast(upper));
}

void SoPlexWrapper::updateRow(int i,const Vec &row,double bVal) {
    DSVector dsv;

    for(int j=0;j<row.size();j++) {
	if(row(j)!=0)
	    dsv.add(j,row(j));
	M(i,j) = row(j);
    }
    b(i) = bVal;

    LPRow lpr(bVal,dsv,infinity);
    plex->changeRow(i,lpr);
}

SoPlexWrapper::SoPlexWrapper() {
    plex = NULL;
    pricer = NULL;
    tester = NULL;

	// Param::setVerbose(0);
}

SoPlexWrapper::~SoPlexWrapper() {
    if(plex) {
	delete plex;
	delete pricer;
	delete tester;
    }
}

void SoPlexWrapper::setProblem(const Vec &c,const Mat &M,const Vec &b,
	const Vec &lower,const Vec &upper) 
{
    SPxLP lp;
    DSVector dsv;		
    int i,j;

    // Delete the old soplex if necessary
    if(plex) {
	delete plex;
	delete pricer;
	delete tester;
    }
    
    // Store the passed in values
    this->M = M;
    this->c = c;
    this->b = b;
    this->u = upper;
    this->l = lower;
    
    // Set up the rows
	LPRowSet rset;
	for(i=0;i<M.rows();i++) {
		dsv.clear();
		for(j=0;j<M.columns();j++) 
			if(M(i,j)!=0.0)
				dsv.add(j,M(i,j));

		// Add sparse vector to row (in greater than fashion)
		rset.add(b(i),dsv,infinity);
	}

    lp.addRows(rset);

    // Set up the columns		
    Vec obj = -c;
    lp.changeObj(vecCast(obj));
    lp.changeBounds(vecCast(lower),vecCast(upper));

    // Initialize the SoPlex problem
    // plex = new SPxSolver(SoPlex::LEAVE,SoPlex::COLUMN);
    plex = new SPxSolver();

    // Set up the components of the solver
    pricer = new SPxDefaultPR();
    tester = new SPxDefaultRT();
    plex->setPricer(pricer);
    plex->setTester(tester);

    plex->loadLP(lp);

}
