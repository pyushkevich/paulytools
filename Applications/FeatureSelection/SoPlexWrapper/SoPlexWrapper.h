#ifndef _SOPLEX_WRAPPER_H_
#define _SOPLEX_WRAPPER_H_

// Inlcude the parent class header
#include "LPWrapper.h"

// Include the necessary SoPlex headers
#include "spxdefines.h"
#include "spxsolver.h"
#include "spxdefaultpr.h"
#include "spxdefaultrt.h"

/**
 * This is an LP Solver based on the SoPlex implementation 
 * by Roland Wunderling, located at 
 * <a href="http://www.zib.de/Optimization/Software/Soplex/soplex.php">
 * 	http://www.zib.de/Optimization/Software/Soplex/soplex.php
 * </a>
 */
class SoPlexWrapper : public LPWrapper {
    public:
	SoPlexWrapper();
	~SoPlexWrapper();

	virtual void setProblem(const Vec &c,const Mat &M,
		const Vec &b,const Vec &l, const Vec &u);

	virtual void updateObjective(const Vec &c);
	virtual void updateBounds(const Vec &u,const Vec &l);
	virtual void updateRow(int i,const Vec &row,double bVal);

	virtual bool solve(Vec &res);

    private:
	soplex::SPxSolver *plex;
	soplex::SPxDefaultPR *pricer;
	soplex::SPxDefaultRT *tester;

	// This method also stores all the parameters because
	// it needs to reinitialize the SoPlex at times
	Mat M;
	Vec c,b,u,l;
};
				
#endif // _SOPLEX_WRAPPER_H_
