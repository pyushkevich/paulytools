#include "EquationSystem.h"
#include "TemplatedWaveletRep.h"

EquationSystem::EquationSystem (const int _dim) : dim(_dim) {
	coeffConst = new double[dim + 3];
	coeff1stOrder = new double*[dim + 3];
	for (int i = 0; i < dim + 3; ++i) {
		coeff1stOrder[i] = new double[3];
	}
	
	// fixed boundary coefficient values
	coeffConst[0] = coeffConst[dim + 2] = 0.0;
	coeff1stOrder[0][0] = 0.0;
	coeff1stOrder[0][1] = -4;
	coeff1stOrder[0][2] = 0.0;
	coeff1stOrder[dim + 2][0] = 0.0;
	coeff1stOrder[dim + 2][1] = -4;
	coeff1stOrder[dim + 2][2] = 0.0;
}

EquationSystem::~EquationSystem () {
	delete[] coeffConst;
	for (int i = 0; i < dim + 3; ++i) {
		delete[] coeff1stOrder[i];
	}
	delete[] coeff1stOrder;
}

void EquationSystem::buildEquationSystem (const WaveletRep& fx, const WaveletRep& fy, const WaveletRep& frho) {
	
	// the step size of parameter t in x(t), y(t) and rho(t)
	const double dt = 1.0/dim;
	
	double t;
	// the discrete x y and their first and second derivative vectors
	double x, y, rho, xt, yt, xtt, ytt;
	// some intermedial coeffs
	double coeff1, coeff2;
	
	// iterate through all the points other than the two boundary points
	for (int j = 1; j < dim + 2; ++j) {
		t = (j - 1)*dt;
		
		// get function values at t
		x = fx.get(t);
		y = fy.get(t);
		rho = frho.get(t);
		xt = fx.get1stDeriv(t);
		yt = fy.get1stDeriv(t);
		xtt = fx.get2ndDeriv(t);
		ytt = fy.get2ndDeriv(t);
		
		coeff1 = 1/(xt*xt+yt*yt);
		coeff2 = -(xt*xtt+yt*ytt)*coeff1*coeff1;
		
		coeffConst[j] = -rho;
		
		double tmp1 = coeff1/(dt*dt);
		double tmp2 = coeff2/(2*dt);
		coeff1stOrder[j][0] = tmp1 - tmp2;
		coeff1stOrder[j][1] = -2 * tmp1;
		coeff1stOrder[j][2] = tmp1 + tmp2;

	}
	
	for (int j = 0; j < 2; ++j) {
		t = j;
		
		// get function values at t
		x = fx.get(t);
		y = fy.get(t);
		rho = frho.get(t);
		xt = fx.get1stDeriv(t);
		yt = fy.get1stDeriv(t);
		
		coeff1 = 1/(xt*xt+yt*yt);
		
		double tmp = coeff1/(4.0*dt*dt);
		coeff2ndOrder[j] = tmp;
		
	}
}

ostream& operator<< (ostream& out, const EquationSystem& es) {
	out << "coeffConst = " << endl;
	for (int i = 0; i < es.dim + 3; ++i) {
		out << es.coeffConst[i] << endl;
	}
	
	out << "coeff1stOrder = " << endl;
	for (int i = 0; i < es.dim + 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			out << es.coeff1stOrder[i][j] << '\t';
		}
		out << endl;
	}
	
	out << "coeff2ndOrder = " << endl;
	for (int i = 0; i < 2; ++i) {
		out << es.coeff2ndOrder[i] << '\t';
	}
	out << flush;
}

