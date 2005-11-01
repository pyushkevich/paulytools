
#ifndef _Wavelet_H
#define _Wavelet_H

class Wavelet {
	private:
	static const double MUL;
	
	public:
	static double psi (const double x);
	static double psi1stDeriv (const double x);
	static double psi2ndDeriv (const double x);
	static double wavelet (const double x, const int j, const int k);
	static double wavelet1stDeriv (const double x, const int j, const int k);
	static double wavelet2ndDeriv (const double x, const int j, const int k);
};

#endif

