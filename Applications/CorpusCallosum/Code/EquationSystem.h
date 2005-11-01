
#ifndef _EquationSystem_H
#define _EquationSystem_H

#include <iostream>

using namespace std;

class WaveletRep;

class EquationSystem {
	public:
	const int dim;
	double *coeffConst;
	double **coeff1stOrder;
	double coeff2ndOrder[2];
		
	public:
	EquationSystem (const int dim);
	~EquationSystem ();
	
	void buildEquationSystem (const WaveletRep& fx, const WaveletRep& fy, const WaveletRep& frho);
	
	friend ostream& operator<< (ostream& out, const EquationSystem& es);
};

#endif

