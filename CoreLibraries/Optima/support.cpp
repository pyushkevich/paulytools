/******************************************************************
 * OPTIMA Library                                                 *
 ******************************************************************
 * Author:						Paul Yushkevich
 *
 * Date:							Apr 1, 1999
 *
 * Description					Multidimensional optimization algorithms
 *									See http://www.cs.unc.edu/~pauly/optima
 *									
 *	Sources:						"Numerical Recepies in C", 
 *									Michaelewitz, "Genetic Algorithms + Data
 *									Structures = Evolutionary Programs"
 *
 * Dependencies:				PY Matrix library, CLAPACK
 ******************************************************************
 * support.cpp
 *	---------
 * Support functions for Optima library 
 ******************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include <mylibs.h>
#include <assert.h>
#include "support.h"

#include <limits>
using namespace std;

// Begin namespace
NAMESPACE_PAULY_START

const double PI = acos(-1.0);

const double randFactor = 1.0 / (1.0 + RAND_MAX);

// N(0,1) generation:
// Wow - Checked this against Excels NORMINV function - 
// we are very close with this crude function!

// Lookup table (goes to 0.5)
// Get Uniform normalized gaussian value
double gaussian(double x) {
	static double gaussianCoeff = 1.0/sqrt(2*PI);
	return gaussianCoeff*exp(-x*x/2);
}

/**
 * Inverse error function 
 */
double dierfc(double y)
{
    double s, t, u, w, x, z;

    z = y;
    if (y > 1) {
        z = 2 - y;
    }
    w = 0.916461398268964 - log(z);
    u = sqrt(w);
    s = (log(u) + 0.488826640273108) / w;
    t = 1 / (u + 0.231729200323405);
    x = u * (1 - s * (s * 0.124610454613712 + 0.5)) - 
        ((((-0.0728846765585675 * t + 0.269999308670029) * t + 
        0.150689047360223) * t + 0.116065025341614) * t + 
        0.499999303439796) * t;
    t = 3.97886080735226 / (x + 3.97886080735226);
    u = t - 0.5;
    s = (((((((((0.00112648096188977922 * u + 
        1.05739299623423047e-4) * u - 0.00351287146129100025) * u - 
        7.71708358954120939e-4) * u + 0.00685649426074558612) * u + 
        0.00339721910367775861) * u - 0.011274916933250487) * u - 
        0.0118598117047771104) * u + 0.0142961988697898018) * u + 
        0.0346494207789099922) * u + 0.00220995927012179067;
    s = ((((((((((((s * u - 0.0743424357241784861) * u - 
        0.105872177941595488) * u + 0.0147297938331485121) * u + 
        0.316847638520135944) * u + 0.713657635868730364) * u + 
        1.05375024970847138) * u + 1.21448730779995237) * u + 
        1.16374581931560831) * u + 0.956464974744799006) * u + 
        0.686265948274097816) * u + 0.434397492331430115) * u + 
        0.244044510593190935) * t - 
        z * exp(x * x - 0.120782237635245222);
    x += s * (x * s + 1);
    if (y > 1) {
        x = -x;
    }
    return x;
}

/**
 * Inverse uniform normal distribution
 */
double norminv(double y) {
	const double _sqrt2 = -sqrt(2.0);
	return _sqrt2 * dierfc(2*y);
}


// Compute a the exponent via lookup tables, in range -5 to 5.  The second parameter
// computeOutOfRange determines whether the values out of range are computed by calling
// exp or whether a constant is returned.  Default is to compute.
double fastExp(double x,bool computeOutOfRange) {
   static bool firstTime = true;
   static const double xMin = -5.0,xMax = 5.0,xStep = 0.001;
   static int tableSize = (int)((xMax-xMin) / xStep) + 1;
   static double *table;
   static double xToIndexMult = (tableSize-1) / (xMax-xMin);

   if(firstTime) {
      // Compute all the values in the table
      table = new double[tableSize];
      int tableIndex = 0;
      for(double xx = xMin;xx < xMax + xStep/2; xx+= xStep) {
         table[tableIndex] = exp(xx);
         tableIndex++;
      }
      firstTime = false;
   }

   if(x >= xMin && x<= xMax) {
      int index = (int)((x - xMin) * xToIndexMult);
      return table[index];
   }

   else if(computeOutOfRange) {
      return exp(x);
   }

   else if(x < xMin) {
      return 0;
   }

   else {
      return 1000000;
   }
}

// Sample the Gaussian distribution via lookup tables.
double getGaussianRnd(double mean,double sigma) {
   static bool firstTime = true;
   static const double xMax = 0.5,xStep = 0.0001;
   static int tableSize = (int)((xMax) / xStep) + 1;
   static double *table;
   static double xToIndexMult = (tableSize-1) / (xMax);

   if(firstTime) {
      // Compute all the values in the table
      table = new double[tableSize];
      double integral = 0.0;
      int tableIndex = 0;
      while(tableIndex < tableSize) {
         table[tableIndex] = integral;
         integral += xStep / gaussian(integral);
         tableIndex++;
      }
      
      firstTime = false;
   }

	double rnd = rand(2.0) - 1.0;
	double offset = (tableSize-1) * fabs(rnd);
	double result = table[(int)offset];
	if(rnd < 0)
		result = -result;
	return sigma*result+mean;
}


// This method is slower than the above because it does not use lookup tables,
// but it is more precise and scientifically sound
double getGaussianRnd2(double mean,double sigma) {
	const double twoPi = 2*PI;
	static double saved = 0.0;
	static bool haveSaved = false;

	if(haveSaved) {
		haveSaved = false;
		return mean + sigma * saved;
	}
	else {
		double u1 = randFactor * (::rand()+1);
		double u2 = randFactor * (::rand()+1);
		double f = sqrt(-2.0*log(u1));
		saved = f*sin(twoPi * u2);
		haveSaved = true;
		double r = mean + sigma * f * cos(2*PI*u2); 

        if(numeric_limits<double>::signaling_NaN() == r) {
            throw "error";
        }
        


		//if(_finite(r))
			return r;
		//else
		//	return mean;
	}
}

// Compute the 'smoothed box' penalty function.  This function
// is similar to the Gaussian, but uses a higher power
// for x in order to penalize values within one std. dev. less
// and ones outside of it more.
double penaltyFunction(double x,double mu,double sigma) {
   static bool firstTime = true;
   static const double xMax = 5.0,xStep = 0.001;
   static int tableSize = (int)((xMax) / xStep) + 1;
   static double *table;
   static double xToIndexMult = (tableSize-1) / (xMax);

   if(firstTime) {
      // Compute all the values in the table
      table = new double[tableSize];
      int tableIndex = 0;
      for(double xx = 0;xx < xMax + xStep/2; xx+= xStep) {
         table[tableIndex] = exp(- xx * xx * xx * xx / 2.0);
         tableIndex++;
      }
      firstTime = false;
   }

   double z = fabs((x-mu) / sigma);

   if(z<= xMax) {
      int index = (int)(z * xToIndexMult);
      return table[index];
   }

   return 0;
}

double penaltyFunction(double x,double mu,double sigma1,double sigma2) {
   double z = fabs((x-mu));
   if(z <= sigma1)
      return 1;
   else
      return penaltyFunction(z-sigma1,0,sigma2-sigma1);
}

// Get a random bit.  The current version assumes that rand_max is 0x7fff which seems to be pretty standard
bool randBit() {
   static short randomShort;
   static short bitIndex = 15;

   if(bitIndex == 15) {
      // Just in case someone compiles this on a screwed up operating system
      assert(2 == sizeof(short));

      // Get a random value
      randomShort = ::rand();

      // Bit index is 0
      bitIndex = 0;
   }

   bool rtn = ((0x0001 & randomShort) == 0x0001) ;
   randomShort >>= 1;
   bitIndex++;

   return rtn;
}



// End namespace
NAMESPACE_PAULY_END


