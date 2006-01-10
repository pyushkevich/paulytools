#include <iostream>
#include <cmath>

#include "CMRep2D.h"
#include "MedialAtom.h"
#include "WaveletRep.h"
#include "BSplineCubUni.h"
#include "EquationSystem.h"
#include "SparseMatrix.h"
#include "NewtonMethod.h"
#include "matrix.h"
using namespace std;

CMRep2D::CMRep2D (const int _dim) : dim(_dim) {
	medialAtoms = new MedialAtom[dim + 1];
	areaOfCMRep = 0.0;
	arcLenVar = 0.0;
}

CMRep2D::~CMRep2D () {
	delete[] medialAtoms;
}


double  CMRep2D::buildCMRep2D (const FunctionRep& _fx, const FunctionRep& _fy, const FunctionRep& _frho, double *phi) {
	// find phi
	EquationSystem c(dim);
	c.buildEquationSystem(_fx, _fy, _frho);
	NewtonMethod solver(c);
	
	solver.initialize(phi);
	solver.solve(c);
	solver.getSolution(phi);
	
	// get the min phi value
	double minPhi = 1.0e10;
	for (int i = 1; i < dim + 2; ++i) {
	  if (phi[i] < minPhi) {
		  minPhi = phi[i];
		}
	}

	// if there exists negative phi, reset phi to normal values for computation next time
	if (minPhi <= 0) {
	  //  cout << "negative phi encounted!" << endl;
	  for (int i = 0; i < dim + 3; ++i){
	    phi[i] = 5.0;
	  }
	  return minPhi;
	}
       
	// compute the medial atoms and get arcLength of medial axis
	double increment = 1.0 / dim;
	double t = 0.0;
	arcLenVar = 0.0;
	double len = 0.0;
	for (int i = 0; i < dim + 1; ++i) {
	  double fx = _fx.get(t);
	  double fy = _fy.get(t);
	  double fxt = _fx.get1stDeriv(t);
	  double fyt = _fy.get1stDeriv(t);
	  double st = sqrt(fxt*fxt + fyt*fyt);
	  // compute r and rt
	  double r = sqrt(phi[i + 1]);
	  double rt = phi[i + 2] - phi[i];
	  rt /= r;
	  rt /= 4;
	  rt *= dim;
	  medialAtoms[i].set(fx, fy, r, fxt, fyt, rt);
	  len += st;
	  t += increment;
	}
	len /= dim;
	t = 0.0;
	for (int i = 0; i < dim + 1; ++i ) {
	  double fxt = _fx.get1stDeriv(t);
	  double fyt = _fy.get1stDeriv(t);
	  double st = sqrt(fxt*fxt + fyt*fyt);
	  double tmp = st - len;
	  tmp = tmp*tmp;
	  arcLenVar +=tmp;
	  t += increment;
	  }
	/*
	for (int i = 0; i < dim + 1; ++i ) {
	  double fxt = _fx.get1stDeriv(t);
	  double fyt = _fy.get1stDeriv(t);
	  double st = fxt*fxt + fyt*fyt;
	  double fxtt = _fx.get2ndDeriv(t);
	  double fytt = _fy.get2ndDeriv(t);
	  double tmp = fxt*fxtt + fyt*fytt;
	  arcLenVar += tmp*tmp/st;
	  t += increment;
	  }*/
	arcLenVar /= dim;
	

	return minPhi;
}

void CMRep2D::getBoundary(Vector &bx, Vector &by) {
  bx.setSize(2*dim + 2);
  by.setSize(2*dim + 2);
  for ( int i = 0; i < dim + 1; ++i) {
    bx(i) = medialAtoms[i].locus[0] + medialAtoms[i].radius * medialAtoms[i].normal1[0];
    by(i) = medialAtoms[i].locus[1] + medialAtoms[i].radius * medialAtoms[i].normal1[1];
    bx(2*dim+1-i) = medialAtoms[i].locus[0] + medialAtoms[i].radius * medialAtoms[i].normal2[0];
    by(2*dim+1-i) = medialAtoms[i].locus[1] + medialAtoms[i].radius * medialAtoms[i].normal2[1];
  }  
}

double CMRep2D::checkBoundaryFold () const {
    // the coordinates of boundary points and 
	double *b1Current = new double[2];
	double *b2Current = new double[2];
	double *b1Next = new double[2];
	double *b2Next = new double[2];
	double *v0 = new double[2];
	double *v1 = new double[2];
	double *v2 = new double[2];
	
	// counter for the boundary reverse points
	double minJac = 1.0e10;
	double negativeJac = 0.0;

	// first atom
	for (int k = 0; k < 2; ++k) {
		b1Current[k] = medialAtoms[0].locus[k] + medialAtoms[0].radius * medialAtoms[0].normal1[k];
		b2Current[k] = medialAtoms[0].locus[k] + medialAtoms[0].radius * medialAtoms[0].normal2[k];
		//		cout << "FirstboundaryPointUpper[" << k+1 << "] = " << b1Current[k] << endl;
		//		cout << "FirstboundaryPointLower[" << k+1 << "] = " << b2Current[k] << endl;
	}

	// rest atom
	for (int i = 1; i < dim+1; ++i) {
		double check1 = 0.0;
		double check2 = 0.0;
		double check3 = 0.0;

		for (int k = 0; k < 2; ++k) {
			b1Next[k] = medialAtoms[i].locus[k] + medialAtoms[i].radius * medialAtoms[i].normal1[k];
			b2Next[k] = medialAtoms[i].locus[k] + medialAtoms[i].radius * medialAtoms[i].normal2[k];

			// vector pointing from current point to next point
			double tmp0 = medialAtoms[i].locus[k] - medialAtoms[i-1].locus[k];
			double tmp1 = b1Next[k] - b1Current[k];
			double tmp2 = b2Next[k] - b2Current[k];
			//	double tmp3 = b1Current[k] -b2Current[k];
			//double tmp4 = b1Next[k] - b2Next[k];

			// the inner product of two vectors
			check1 += tmp0*tmp1;
			check2 += tmp0*tmp2;
			//	check3 +=tmp3*tmp4;
			// update current value
			b1Current[k] = b1Next[k];
			b2Current[k] = b2Next[k];
		}
		if (check1 < 0) {
		  negativeJac += check1;
		}
		if (check2 < 0) {
		  negativeJac += check2;
		}
		//		if (check3 < 0) {
		// negativeJac += check3;
		//}
		if(check1 < minJac) {
		  minJac = check1;
		}
		if(check2 < minJac) {
		  minJac = check2;
		}
		//	if(check3 < minJac) {
		//  minJac = check3;
		//	}
	}
	if (minJac < 0) {
	  // cout << "boundary folds!" << endl;
	  return negativeJac;
	}

	delete[] b1Current;
	delete[] b2Current;
	delete[] b1Next;
	delete[] b2Next;
	delete[] v0;
	delete[] v1;
	delete[] v2;

	return minJac;
}

double CMRep2D::computeAreaOverlap (const ImageInterpolator *image, const double ratio, const int n, const int ne) {
	// total points to sample
	const int N = n + ne + 1;
	
	// the coordinates along radius to sample image values [0.0, 1.0]
	double *radial = new double[N];
	// (n + 1) points [0.0, ratio]
	double mul = ratio;
	mul /= n;
	for (int i = 0; i <= n; ++i) {
		radial[i] = i * mul;
	}
	// ne points (ratio, 1.0]
	mul = 1 - ratio;
	mul /= ne;
	for (int i = 1; i <= ne; ++i) {
		radial[i + n] = ratio + i * mul;
	}
	
	// the coordinates in cartisean
	double **u1Current = new double *[N];
	double **u2Current = new double *[N];
	double **u1Next = new double *[N];
	double **u2Next = new double *[N];
	for (int i = 0; i < N; ++i) {
		u1Current[i] = new double[2];
		u2Current[i] = new double[2];
		u1Next[i] = new double[2];
		u2Next[i] = new double[2];
	}
	
	// the corresponding image values:
	double *v1Current = new double[N];
	double *v2Current = new double[N];
	double *v1Next = new double[N];
	double *v2Next = new double[N];
	
	// first atom
	for (int j = 0; j < N; ++j) {
		for (int k = 0; k < 2; ++k) {
			u1Current[j][k] = medialAtoms[0].locus[k] + medialAtoms[0].radius * radial[j] * medialAtoms[0].normal1[k];
			u1Current[j][k] += 128.0;
			u2Current[j][k] = medialAtoms[0].locus[k] + medialAtoms[0].radius * radial[j] * medialAtoms[0].normal2[k];
			u2Current[j][k] += 128.0;
		}
		// look up the image values, use 6-250 pixels because there may be errors near boundary, meantime convert binary (inside, outside)=(1,0) to (1,-1)
		if (u1Current[j][0] > 6.0 && u1Current[j][0] < 250.0 && u1Current[j][1] > 6.0 && u1Current[j][1] < 250.0) {
			v1Current[j] = 2.0*image->Evaluate(u1Current[j])-1.0;
		} else {
			v1Current[j] = -1.0;
		}
		if (u2Current[j][0] > 6.0 && u2Current[j][0] < 250.0 && u2Current[j][1] > 6.0 && u2Current[j][1] < 250.0) {
			v2Current[j] = 2.0*image->Evaluate(u2Current[j])-1.0;
		} else {
			v2Current[j] = -1.0;
		}
	}
	
	double areaOverlap = 0.0;
	// rest of the atoms
	for (int i = 1; i < dim + 1; ++i) {
		for (int j = 0; j < N; ++j) {
			for (int k = 0; k < 2; ++k) {
				u1Next[j][k] = medialAtoms[i].locus[k] + medialAtoms[i].radius * radial[j] * medialAtoms[i].normal1[k];
				u1Next[j][k] += 128.0;
				u2Next[j][k] = medialAtoms[i].locus[k] + medialAtoms[i].radius * radial[j] * medialAtoms[i].normal2[k];
				u2Next[j][k] += 128.0;
			}
			// look up the image values, use 6-250 pixels because there may be errors near boundary, meantime convert binary (inside, outside)=(1,0) to (1,-1)
			if (u1Next[j][0] > 6.0 && u1Next[j][0] < 250.0 && u1Next[j][1] > 6.0 && u1Next[j][1] < 250.0) {
				v1Next[j] = 2.0*image->Evaluate(u1Next[j])-1.0;
			} else {
				v1Next[j] = -1.0;
			}
			if (u2Next[j][0] > 6.0 && u2Next[j][0] < 250.0 && u2Next[j][1] > 6.0 && u2Next[j][1] < 250.0) {
				v2Next[j] = 2.0*image->Evaluate(u2Next[j])-1.0;
			} else {
				v2Next[j] = -1.0;
			}
		}
		
		// compute the area overlap
		for (int j = 0; j < N - 1; ++j) {
			// lower triangle
			double a0 = u1Current[j + 1][0] - u1Current[j][0];
			double a1 = u1Current[j + 1][1] - u1Current[j][1];
			double b0 = u1Next[j][0] - u1Current[j][0];
			double b1 = u1Next[j][1] - u1Current[j][1];
			double tmp1 = fabs(a0*b1 - a1*b0);
			tmp1 *= 0.5;
			// upper triangle
			a0 = u1Next[j][0] - u1Next[j + 1][0];
			a1 = u1Next[j][1] - u1Next[j + 1][1];
			b0 = u1Current[j + 1][0] - u1Next[j + 1][0];
			b1 = u1Current[j + 1][1] - u1Next[j + 1][1];
			double tmp2 = fabs(a0*b1 - a1*b0);
			tmp2 *= 0.5;
			double tmp = tmp1 + tmp2;
			tmp *= v1Current[j] + v1Current[j + 1] + v1Next[j] + v1Next[j + 1];
			tmp *= 0.25;
			areaOverlap += tmp;
			
			// lower triangle
			a0 = u2Current[j + 1][0] - u2Current[j][0];
			a1 = u2Current[j + 1][1] - u2Current[j][1];
			b0 = u2Next[j][0] - u2Current[j][0];
			b1 = u2Next[j][1] - u2Current[j][1];
			tmp1 = fabs(a0*b1 - a1*b0);
			tmp1 *= 0.5;
			// upper triangle
			a0 = u2Next[j][0] - u2Next[j + 1][0];
			a1 = u2Next[j][1] - u2Next[j + 1][1];
			b0 = u2Current[j + 1][0] - u2Next[j + 1][0];
			b1 = u2Current[j + 1][1] - u2Next[j + 1][1];
			tmp2 = fabs(a0*b1 - a1*b0);
			tmp2 *= 0.5;
			tmp = tmp1 + tmp2;
			tmp *= v2Current[j] + v2Current[j + 1] + v2Next[j] + v2Next[j + 1];
			tmp *= 0.25;
			areaOverlap += tmp;
		}
		
		// copy next to current
		for (int j = 0; j < N; ++j) {
			for (int k = 0; k < 2; ++k) {
				u1Current[j][k] = u1Next[j][k];
				u2Current[j][k] = u2Next[j][k];
			}
			v1Current[j] = v1Next[j];
			v2Current[j] = v2Next[j];
		}
	}
	
	// clean up
	for (int i = 0; i < N; ++i) {
		delete[] u1Current[i];
		delete[] u2Current[i];
		delete[] u1Next[i];
		delete[] u2Next[i];
	}
	delete[] u1Current;
	delete[] u2Current;
	delete[] u1Next;
	delete[] u2Next;
	delete[] v1Current;
	delete[] v2Current;
	delete[] v1Next;
	delete[] v2Next;
	delete[] radial;
	
	return areaOverlap;
}


double CMRep2D::overlapAndCMRep (const ImageInterpolator *image, const double ratio, const int n, const int ne, const double dx) {
	// total points to sample
	const int N = n + ne + 1;
	
	// the coordinates along radius to sample image values [0.0, 1.0]
	double *radial = new double[N];
	// (n + 1) points [0.0, ratio]
	double mul = ratio;
	mul /= n;
	for (int i = 0; i <= n; ++i) {
		radial[i] = i * mul;
	}
	// ne points (ratio, 1.0]
	mul = 1 - ratio;
	mul /= ne;
	for (int i = 1; i <= ne; ++i) {
		radial[i + n] = ratio + i * mul;
	}
	
	// the coordinates in cartisean
	double **u1Current = new double *[N];
	double **u2Current = new double *[N];
	double **u1Next = new double *[N];
	double **u2Next = new double *[N];
	for (int i = 0; i < N; ++i) {
		u1Current[i] = new double[2];
		u2Current[i] = new double[2];
		u1Next[i] = new double[2];
		u2Next[i] = new double[2];
	}
	
	// the corresponding image values:
	double *v1Current = new double[N];
	double *v2Current = new double[N];
	double *v1Next = new double[N];
	double *v2Next = new double[N];
	
	// first atom
	for (int j = 0; j < N; ++j) {
		for (int k = 0; k < 2; ++k) {
			u1Current[j][k] = medialAtoms[0].locus[k] + medialAtoms[0].radius * radial[j] * medialAtoms[0].normal1[k] + dx*(1.0-k);
			u1Current[j][k] += 128.0;
			u2Current[j][k] = medialAtoms[0].locus[k] + medialAtoms[0].radius * radial[j] * medialAtoms[0].normal2[k] + dx*(1.0-k);
			u2Current[j][k] += 128.0;
		}
		// look up the image values, use 6-250 pixels because there may be errors near boundary, meantime convert binary (inside, outside)=(1,0) to (1,-1)
		if (u1Current[j][0] > 6.0 && u1Current[j][0] < 250.0 && u1Current[j][1] > 6.0 && u1Current[j][1] < 250.0) {
			v1Current[j] = 2.0*image->Evaluate(u1Current[j])-1.0;
		} else {
			v1Current[j] = -1.0;
		}
		if (u2Current[j][0] > 6.0 && u2Current[j][0] < 250.0 && u2Current[j][1] > 6.0 && u2Current[j][1] < 250.0) {
			v2Current[j] = 2.0*image->Evaluate(u2Current[j])-1.0;
		} else {
			v2Current[j] = -1.0;
		}
	}
	
	double areaOverlap = 0.0;
	areaOfCMRep = 0.0;
	// rest of the atoms
	for (int i = 1; i < dim + 1; ++i) {
		for (int j = 0; j < N; ++j) {
			for (int k = 0; k < 2; ++k) {
				u1Next[j][k] = medialAtoms[i].locus[k] + medialAtoms[i].radius * radial[j] * medialAtoms[i].normal1[k] + dx*(1.0-k);
				u1Next[j][k] += 128.0;
				u2Next[j][k] = medialAtoms[i].locus[k] + medialAtoms[i].radius * radial[j] * medialAtoms[i].normal2[k] + dx*(1.0-k);
				u2Next[j][k] += 128.0;
			}
			// look up the image values, use 6-250 pixels because there may be errors near boundary, meantime convert binary (inside, outside)=(1,0) to (1,-1)
			if (u1Next[j][0] > 6.0 && u1Next[j][0] < 250.0 && u1Next[j][1] > 6.0 && u1Next[j][1] < 250.0) {
				v1Next[j] = 2.0*image->Evaluate(u1Next[j])-1.0;
			} else {
				v1Next[j] = -1.0;
			}
			if (u2Next[j][0] > 6.0 && u2Next[j][0] < 250.0 && u2Next[j][1] > 6.0 && u2Next[j][1] < 250.0) {
				v2Next[j] = 2.0*image->Evaluate(u2Next[j])-1.0;
			} else {
				v2Next[j] = -1.0;
			}
		}
		
		// compute the area overlap
		for (int j = 0; j < N - 1; ++j) {
			// lower triangle
			double a0 = u1Current[j + 1][0] - u1Current[j][0];
			double a1 = u1Current[j + 1][1] - u1Current[j][1];
			double b0 = u1Next[j][0] - u1Current[j][0];
			double b1 = u1Next[j][1] - u1Current[j][1];
			double tmp1 = fabs(a0*b1 - a1*b0);
			tmp1 *= 0.5;
			// upper triangle
			a0 = u1Next[j][0] - u1Next[j + 1][0];
			a1 = u1Next[j][1] - u1Next[j + 1][1];
			b0 = u1Current[j + 1][0] - u1Next[j + 1][0];
			b1 = u1Current[j + 1][1] - u1Next[j + 1][1];
			double tmp2 = fabs(a0*b1 - a1*b0);
			tmp2 *= 0.5;
			double tmp = tmp1 + tmp2;
			areaOfCMRep += tmp;
			tmp *= v1Current[j] + v1Current[j + 1] + v1Next[j] + v1Next[j + 1];
			tmp *= 0.25;
			areaOverlap += tmp;
			
			// lower triangle
			a0 = u2Current[j + 1][0] - u2Current[j][0];
			a1 = u2Current[j + 1][1] - u2Current[j][1];
			b0 = u2Next[j][0] - u2Current[j][0];
			b1 = u2Next[j][1] - u2Current[j][1];
			tmp1 = fabs(a0*b1 - a1*b0);
			tmp1 *= 0.5;
			// upper triangle
			a0 = u2Next[j][0] - u2Next[j + 1][0];
			a1 = u2Next[j][1] - u2Next[j + 1][1];
			b0 = u2Current[j + 1][0] - u2Next[j + 1][0];
			b1 = u2Current[j + 1][1] - u2Next[j + 1][1];
			tmp2 = fabs(a0*b1 - a1*b0);
			tmp2 *= 0.5;
			tmp = tmp1 + tmp2;
			areaOfCMRep += tmp;
			tmp *= v2Current[j] + v2Current[j + 1] + v2Next[j] + v2Next[j + 1];
			tmp *= 0.25;
			areaOverlap += tmp;
		}
		
		// copy next to current
		for (int j = 0; j < N; ++j) {
			for (int k = 0; k < 2; ++k) {
				u1Current[j][k] = u1Next[j][k];
				u2Current[j][k] = u2Next[j][k];
			}
			v1Current[j] = v1Next[j];
			v2Current[j] = v2Next[j];
		}
	}
	
	// clean up
	for (int i = 0; i < N; ++i) {
		delete[] u1Current[i];
		delete[] u2Current[i];
		delete[] u1Next[i];
		delete[] u2Next[i];
	}
	delete[] u1Current;
	delete[] u2Current;
	delete[] u1Next;
	delete[] u2Next;
	delete[] v1Current;
	delete[] v2Current;
	delete[] v1Next;
	delete[] v2Next;
	delete[] radial;
	
	return areaOverlap;
}

double CMRep2D::getAreaOfCMRep() {
  return areaOfCMRep;
}

double CMRep2D::getArcLenVar() {
  return arcLenVar;
}

