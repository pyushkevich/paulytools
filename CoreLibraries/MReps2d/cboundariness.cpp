/******************************************************************
 * MMODEL Library                                                 *
 ******************************************************************
 * Author:						Paul Yushkevich
 *
 * Date:							Feb 15, 1999
 *
 * Description					Medial model algorithms and data structures
 *
 *
 *	Sources:                Various papers on medial models at UNC,
 *                         Pythagoreah hodographs, etc.
 *
 * Dependencies:				PY Matrix, Optima, Registry libs, CLAPACK
 ******************************************************************
 * cboundariness.cpp
 *	-----------------
 * Continous boundariness measure
 ******************************************************************/
#include <math.h>
#include "cboundariness.h"

// Begin namespace
NAMESPACE_PAULY_START

const int ktStep = 1000;
const double ktStdDev = 3.0;
const int ktSize = (int)(2*ktStdDev*ktStep)+1;
double *kernelTable;
bool tableBuilt = false;

void buildKernelTable() {
	double step = 1.0 / ktStep;
	double mult = 2 / sqrt(2*M_PI);
	int index = 0;

	kernelTable = new double[ktSize];

	for(double x=-ktStdDev;x<=ktStdDev;x+=step) {
		kernelTable[index] = mult*x*exp(-x*x);
		index++;
	}
}

double getContinuousBoundariness(MedialInterpolant &mi,cimage im,double steps,double rho) {	
	if(!tableBuilt) {
		buildKernelTable();
		tableBuilt=true;
	}

	// Standard deviation computation
	double sigma = sqrt(mi.start.s() * mi.end.s()) * rho;
	double mult = 1.0 / (sigma*sigma);

	// Later we will take steps along the normal that are half a pixel long
	int res = (aux_cimage_xdim(im) > aux_cimage_ydim(im)) ? 
					aux_cimage_xdim(im) : aux_cimage_ydim(im);
	double nStepLength = 1.0 / (2*res);
	double zStepSize = nStepLength / sigma;

	// Find total sum of boundariness
	double boundariness = 0.0;

	// Apply kernel every so ofter along the boundary track
	double stepSize = 1.0 / steps;
	for(double t=0;t<=1.0;t+=stepSize) {
		for(int side=0;side<2;side++) {
			double x,y,nx,ny;
			mi.q[side]->getInterpolation(t,x,y,nx,ny);

			double nLength = sqrt(nx*nx+ny*ny);
			nx /= nLength;ny /= nLength;


			// A local boundariness measure
			double localBoundariness = 0.0;

			// Take steps along the normal vector
			double nxStepLength = nx*zStepSize*sigma*res;
			double nyStepLength = ny*zStepSize*sigma*res;
			x = res * (-ktStdDev*nx*sigma + x);
			y = res * (-ktStdDev*ny*sigma + y);

			// Take steps in the table too
			double indexStep = zStepSize*ktStep;
			double index = 0.0;
			for(double z=-ktStdDev;z<=ktStdDev;z+=zStepSize) {
				double intensity = aux_cimage_get(im,(int)x,(int)-y);
				double kernelValue = kernelTable[(int)index];
				localBoundariness += kernelValue*intensity;

				x += nxStepLength;y += nyStepLength;
				index += indexStep;
			}

			boundariness += fabs(localBoundariness);
		}
	}

	boundariness *= mult;
	return -boundariness;
}


double getContinuousEndness(EndInterpolant &mi,cimage im,double steps,double rho) {
	if(!tableBuilt) {
		buildKernelTable();
		tableBuilt=true;
	}

	// Standard deviation computation
	double sigma = mi.mp.s() * rho;
	double mult = 1.0 / (sigma*sigma);

	// Later we will take steps along the normal that are half a pixel long
	int res = (aux_cimage_xdim(im) < aux_cimage_ydim(im)) ?
					aux_cimage_xdim(im) : aux_cimage_ydim(im);
	double nStepLength = 1.0 / (2*res);
	double zStepSize = nStepLength / sigma;

	// Find total sum of boundariness
	double boundariness = 0.0;

	// Apply kernel every so ofter along the boundary track
	double stepSize = 1.0 / steps;
	for(double t=0;t<=1.0;t+=stepSize) {
		for(int side=0;side<2;side++) {
			double x,y,nx,ny;
			mi.q[side]->getInterpolation(t,x,y,nx,ny);

			double nLength = sqrt(nx*nx+ny*ny);
			nx /= nLength;ny /= nLength;


			// A local boundariness measure
			double localBoundariness = 0.0;

			// Take steps along the normal vector
			double nxStepLength = nx*zStepSize*sigma*res;
			double nyStepLength = ny*zStepSize*sigma*res;
			x = res * (-ktStdDev*nx*sigma + x);
			y = res * (-ktStdDev*ny*sigma + y);

			// Take steps in the table too
			double indexStep = zStepSize*ktStep;
			double index = 0.0;
			for(double z=-ktStdDev;z<=ktStdDev;z+=zStepSize) {
            // This should be treated in a better fashion.  For now, the
            // images are assumed to have a dark border.
            if(x < 0 || x >= res || y > 0 || -y >= res)
               continue;

				double intensity = aux_cimage_get(im,(int)x,(int)-y);
				double kernelValue = kernelTable[(int)index];
				localBoundariness += kernelValue*intensity;

				x += nxStepLength;y += nyStepLength;
				index += indexStep;
			}
			
			boundariness += fabs(localBoundariness);
		}
	}

	boundariness *= mult;
	return -boundariness;
}





// End namespace
NAMESPACE_PAULY_END

