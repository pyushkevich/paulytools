#include "ispace.h"
#include <optima.h>

#include <fftw3.h>
#include <iostream>
using namespace std;

extern "C" {
  CFFT2DC(void *,void *,int,int,int);
}

/**************************************************************
ARRAY 2D Template - For fast image storage,
no bells or whistles.  Can use any primitive for storage
*************************************************************/
#ifndef M_PI
#define M_PI   3.14159265358979
#endif

/**************************************************************
Scale Space Implementation
************************************************************
void ScaleSpace::computeBlur(int slice) {
   // Spatial scale
   double s = slices[slice].scale;
   
   // Print out which slice we are computing
   printf("computing blur slice %d\n",slice);
   
   // Fourier scale
   double sfx = w / (2*M_PI*s);
   double sfy = h / (2*M_PI*s);
   double denomX = -2.0 * sfx * sfx;
   double denomY = -2.0 * sfy * sfy;
   
   // Helper variables
   int halfx = w / 2;
   int halfy = h / 2;
   int x1,y1;
   
   // First, create a new transform array
   // Allocate the data to perform the transform on
   fftw_complex *work = new fftw_complex [h*w];
   
   int index = 0;
   for(int y=0;y<h;y++) {
      y1 = (y < halfy) ? y : h - y;
      for(int x=0;x<w;x++) {
         x1 = (x < halfx) ? x : w - x;
         
         // Calculate value of Gaussian 
         double value = appExp(x1 * x1 / denomX + y1 * y1 / denomY);
         
         work[index].re = baseFFT[index].re * value;
         work[index].im = baseFFT[index].im * value;
         index++;
      }
   }
   
   // Now, compute the inverse transform
   // First, take a forward fourier transform of our image
   fftwnd_plan p = fftw2d_create_plan(w, h, FFTW_BACKWARD, FFTW_IN_PLACE);
   fftwnd_one(p,work,NULL);
   fftwnd_destroy_plan(p);
   
   // Create a new cimage and copy the results to it
   slices[slice].blur.resize(w,h);
   
   index = 0;
   double mult = 1.0 / (w*h);
   for(y=0;y<h;y++) {
      for(int x=0;x<w;x++) {
         double val = work[index].re*mult;
         slices[slice].blur(x,y) = val;
         index++;
      }
   }
   
   delete work;
}
*/

/*
void ScaleSpace::computeGrad(int slice) {
   // Spatial scale
   double s = slices[slice].scale;
   
   // Print out which slice we are computing
   printf("computing gradient slice %d\n",slice);
   
   

   // First, take a forward fourier transform of our image

   fftw_complex 
      *fxFFT = new fftw_complex [h*w],
      *fyFFT = new fftw_complex [h*w],
      *fFFT = new fftw_complex [h*w];
   
   for(int y=0;y<h;y++) {
      int yc = y < h/2 ? y : y - h;
      for(int x=0;x<w;x++) {
         int xc = x < w/2 ? x : x - w;

         int index = y*w+x;
         fFFT[index].re = exp(-(xc*xc+yc*yc)/(2*s*s))/(2*M_PI*s*s);
         fFFT[index].im = 0;
         fxFFT[index].re = fFFT[index].re * (-xc/(s*s));
         fxFFT[index].im = 0;
         fyFFT[index].re = fFFT[index].re * (-yc/(s*s));
         fyFFT[index].im = 0;
      }
   }

   
   fftwnd_plan p0 = fftw2d_create_plan(w, h, FFTW_FORWARD, FFTW_IN_PLACE);
   fftwnd_one(p0,fFFT,NULL);
   fftwnd_one(p0,fxFFT,NULL);
   fftwnd_one(p0,fyFFT,NULL);
   fftwnd_destroy_plan(p0);
   

   // Fourier scale
   double sfx = w / (2*M_PI*s);
   double sfy = h / (2*M_PI*s);
   double denomX = -2.0 * sfx * sfx;
   double denomY = -2.0 * sfy * sfy;
   
   // Helper variables
   int halfx = w / 2;
   int halfy = h / 2;
   int x1,y1;
   
   // Create a couple arrays, one for horizontal derivative, one for vertical
   fftw_complex *wx = new fftw_complex [h*w];
   fftw_complex *wy = new fftw_complex [h*w];
   
   int index = 0;
   for(y=0;y<h;y++) {
      int yFreq = (y < halfy) ? y : y - h;
      int yGauss = (y < halfy) ? y : h - y;

      for(int x=0;x<w;x++) {
         int xFreq = (x < halfx) ? x : x - w;
         int xGauss = (x < halfx) ? x : w - x;
         
         // Calculate value of Gaussian 
         double value = exp(xFreq * xFreq / denomX + yFreq * yFreq / denomY);
         double fx = value * 2 * M_PI * xFreq / w;
         double fy = value * 2 * M_PI * yFreq / h;
   
         //
         double vf = fFFT[index].re;
         double vfx = fxFFT[index].im;
         double vfy = fyFFT[index].im;
             
         wx[index].re = 0;
         wx[index].im = baseFFT[index].re * fx;
         wy[index].re = 0;
         wy[index].im = baseFFT[index].re * fy;

         index++;
      }
   }
   
   // Now, compute the inverse transform
   // First, take a forward fourier transform of our image
   fftwnd_plan p = fftw2d_create_plan(w, h, FFTW_BACKWARD, FFTW_IN_PLACE);
   fftwnd_one(p,wx,NULL);
   fftwnd_destroy_plan(p);
   
   p = fftw2d_create_plan(w, h, FFTW_BACKWARD, FFTW_IN_PLACE);
   fftwnd_one(p,wy,NULL);
   fftwnd_destroy_plan(p);
   
   // Create a new cimage and copy the results to it
   slices[slice].grad.resize(w,h);
   
   index = 0;
   double mult = 1.0 / (w*h);
   for(y=0;y<h;y++) {
      for(int x=0;x<w;x++) {
         slices[slice].grad(x,y).x = wx[index].re*mult;
         slices[slice].grad(x,y).y = wy[index].re*mult;
         index++;
      }
   }
   
   delete wx;
   delete wy;
}

*/

void ScaleSpace::computeJet(int slice) {
   int x,y;

   // Spatial scale
   float s = sliceScales[slice];
   
   // Print out which slice we are computing
   printf("computing jet at scale %lg\n",s);
  
   // Fourier scale
   float sfx = (float) (w2 / (2*M_PI*s));
   float sfy = (float) (h2 / (2*M_PI*s));

   float sqrt2Pi = (float) sqrt(M_PI*2);

   // Factor to multiply x and y by
   float multX = (float) (1.0 / (-2.0 * sfx * sfx));
   float multY = (float) (1.0 / (-2.0 * sfy * sfy));
   
   // Helper variables
   int halfx = w2 / 2;
   int halfy = h2 / 2;
   
   // First we will compute the exponent term for all x and y seprately
   // because each Gaussian is just a product of horizontal and vertical terms
   float *xCache = new float[w2];
   float *yCache = new float[h2];
   float *xDerCache = new float[w2];
   float *yDerCache = new float[h2];

   for(x=0;x<=halfx;x++) {
      xCache[x] = (float)exp(x*x*multX);
      xDerCache[x] = sqrt2Pi * x / sfx;
   }
   for(x=halfx+1;x<w2;x++) {
      xCache[x] = xCache[w-x];
      xDerCache[x] = -xDerCache[w-x];
   }
   for(y=0;y<=halfy;y++) {
      yCache[y] = (float) exp(y*y*multY);
      yDerCache[y] = sqrt2Pi * y / sfy;
   }
   for(y=halfy+1;y<h2;y++) {
      yCache[y] = yCache[h-y];
      yDerCache[y] = -yDerCache[h-y];
   }

   // Create storage array for fourier transforms
   fftw_complex *jet[6];
   fftw_plan plan[6];

   for(int i=0;i<6;i++) {
      jet[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * w2 * h2);
      plan[i] = fftw_plan_dft_2d(w2,h2,jet[i],jet[i],FFTW_BACKWARD, FFTW_ESTIMATE);
   }
         
   // Now compute all filters as products of their buddies
   int offset = 0;
   for(y = 0;y < h2; y++) {
      for(x = 0;x < w2; x++) {
         float g = xCache[x]*yCache[y];
         float gx = g * xDerCache[x];  // times i 
         float gy = g * yDerCache[y];  // times i
         float gxx = xDerCache[x] * gx;   // times -1
         float gxy = xDerCache[x] * gy;   // times -1
         float gyy = yDerCache[y] * gy;   // times -1

         float reBase = baseFFT[offset][0];
         float imBase = baseFFT[offset][1];

         jet[0][offset][0] = g*reBase;
         jet[0][offset][1] = g*imBase;
         
         jet[1][offset][0] = -gx * imBase;
         jet[1][offset][1] = gx * reBase;
         jet[2][offset][0] = -gy * imBase;
         jet[2][offset][1] = gy * reBase;

         jet[3][offset][0] = -gxx*reBase;
         jet[3][offset][1] = -gxx*imBase;
         jet[4][offset][0] = -gxy*reBase;
         jet[4][offset][1] = -gxy*imBase;
         jet[5][offset][0] = -gyy*reBase;
         jet[5][offset][1] = -gyy*imBase;
  
         offset++;
      }
   }

   delete[] xCache;
   delete[] yCache;
   delete[] xDerCache;
   delete[] yDerCache;
   
   // Now, compute the inverse transform
   // First, take a forward fourier transform of our image
   for(i=0;i<6;i++) {
      fftw_execute(plan[i]); 
   }

    // Create a new cimage and copy the results to it
   // float mult = 1.0f / (w2*h2);
   for(i=0;i<6;i++) {
      jets[slice][i].resize(w,h);

      // Copy from the FFT output to the slice
      for(y=0;y<h;y++) {
         int offset = w2*y;
         for(x=0;x<w;x++) {
            jets[slice][i](x,y) = (float) jet[i][offset][0];
            offset++;
         }
      }

      // Clean up the fft array
      fftw_destroy_plan(plan[i]);
      fftw_free(jet[i]);
   }
}


// Performs interpolation for whole jet or parts of the jet
void ScaleSpace::sampleJet(const Vector2D &p,double s,int first,int last) {
   
   // Check bounds - we are currently wrapping around the image
   if(s < 0.001 || s > 1000 || p.x < 0 || p.x >= w-1 || p.y < 0 || p.y >= h-1) {
      // Out of bounds.  Return 0s
      for(int i=first;i<=last;i++)
         jetlet[i] = 0;
      return;
   }

	// Get slices for the scale
	int sl = findSlice((float)s);
	int su = sl + 1;
	float sScale = 1.0f / (sliceScales[su] - sliceScales[sl]);
	
	float s1 = (float) (sScale * (s - sliceScales[sl]));
	float s2 = (float) (sScale * (sliceScales[su] - s));
   
   // Check for slice existance
   if(jets[sl][0].width()==0)
		computeJet(sl);
   
	if(jets[su][0].width()==0)
		computeJet(su);
   
   // Compute the space bilerp parameters
   int xl = (int)p.x;
   int xu = xl + 1;
   float x1 = (float)p.x - xl;
   float x2 = 1.0f - x1;

   int yl = (int)p.y;
   int yu = yl + 1;
   float y1 = (float)p.y - yl;
   float y2 = 1.0f - y1;

   float f11 = s1*x1;
   float f12 = s1*x2;
   float f21 = s2*x1;
   float f22 = s2*x2;
   float f111 = f11*y1;
   float f112 = f11*y2;
   float f121 = f12*y1;
   float f122 = f12*y2;
   float f211 = f21*y1;
   float f212 = f21*y2;
   float f221 = f22*y1;
   float f222 = f22*y2;
   
   for(int i=first;i<=last;i++) {
      jetlet[i] =
         jets[su][i](xu,yu) * f111 +
         jets[su][i](xu,yl) * f112 + 
         jets[su][i](xl,yu) * f121 + 
         jets[su][i](xl,yl) * f122 +
         jets[sl][i](xu,yu) * f211 +
         jets[sl][i](xu,yl) * f212 +
         jets[sl][i](xl,yu) * f221 +
         jets[sl][i](xl,yl) * f222;
   }
}



int roundUpToPowerOfTwo(int x) {
	int X = x;

	if(--x < 1)
		return 1;

	int log2 = -1;
	while(x != 0) {
		x >>= 1;
		log2++;
	}

	x = 1;
	for(int i=0;i<=log2;i++) {
		x <<= 1;
	}
		
	cout << "round up of " << X << " is " << x;
	return x;
}

void ScaleSpace::computeScaleSpace() 
{

	// Create an image of an appropriate size
	w2 = roundUpToPowerOfTwo(w);
	h2 = roundUpToPowerOfTwo(h);

   // Create storage array for fourier transforms
   if(baseFFT) fftw_free(baseFFT);
   baseFFT = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * w2 * h2);   
   fftw_plan plan = fftw_plan_dft_2d(w2,h2,baseFFT,baseFFT,FFTW_FORWARD, FFTW_ESTIMATE);
   
	// Copy the image into the real component
	int offset = 0;
   for(int j=0;j<h;j++) {
      for(int i=0;i<w;i++) {
		   baseFFT[offset][0] = base(i,j);
			baseFFT[offset][1] = 0;
         ++offset;
		}
	}
	
   // Take a forward fourier transform of our image
   fftw_execute(plan);
   fftw_destroy_plan(plan);
}

// const int ScaleSpace::nSlices = 60;
// const int ScaleSpace::sliceZero = 30;
const double ScaleSpace::logSqrt2 = log(sqrt(2.0));


ScaleSpace::ScaleSpace(cimage &im) {

	int i;

	w = roundUpToPowerOfTwo(cimage_xdim(im));
	h = roundUpToPowerOfTwo(cimage_ydim(im));
   
   // Copy the incoming image to storage
   base.resize(w,h);
   baseFFT=NULL;

   iMax = cimage_get(im,0,0),iMin = cimage_get(im,0,0);
   for(i=0;i<base.width();i++) {
      for(int j=0;j<base.height();j++) {
         if(i < cimage_xdim(im) && j < cimage_xdim(im))
				base(i,j) = cimage_get(im,i,j);
			else 
				base(i,j) = (float) iMin;

         iMax = (iMax < base(i,j)) ? base(i,j) : iMax;
         iMin = (iMin > base(i,j)) ? base(i,j) : iMin;
      }
   }
   
   for(i=0;i<base.width();i++) {
      for(int j=0;j<base.height();j++) {
         base(i,j) = (float) (255.0 * (base(i,j) - iMin) / (iMax-iMin));
      }
   }
   
	// Precompute the scales of all slices
   for(i=0;i<nSlices;i++) {
      sliceScales[i] = (float) sliceToScale(i);
   }
   
	// Initialize the binary searcher
	fSlice = sliceZero;

   // Perform forward FFT
   computeScaleSpace();
}

ScaleSpace::~ScaleSpace() {
   if(baseFFT) fftw_free(baseFFT);
}


/**************************************************************
Exponent table
*************************************************************
double appExp(double x) {
   static ExponentTable *expTable = NULL;
   
   if(expTable==NULL) {
      expTable = new ExponentTable(0.0001);
   }
   
   return expTable->get(x);
}
*/

/**************************************************************
 Traditional image space,
 Blurring is done each time on demand
 **************************************************************/
CImageSpace::CImageSpace(cimage &im) {
   // Copy the incoming image to storage
   image.resize(cimage_xdim(im),cimage_ydim(im));
   
   double iMax = cimage_get(im,0,0),iMin = cimage_get(im,0,0);
   for(int i=0;i<image.width();i++) {
      for(int j=0;j<image.height();j++) {
         image(i,j) = cimage_get(im,i,j);
         iMax = (iMax < image(i,j)) ? image(i,j) : iMax;
         iMin = (iMin > image(i,j)) ? image(i,j) : iMin;
      }
   }

   for(i=0;i<image.width();i++) {
      for(int j=0;j<image.height();j++) {
         image(i,j) = 255.0 * (image(i,j) - iMin) / (iMax-iMin);
      }
   }

   
   // Cache width and height
   w = image.width();
   h = image.height();

   // Allocate exponent array
   xExpArray = new double[w];
}

CImageSpace::~CImageSpace() {
   delete xExpArray;
}

// Get the value of the image blurred at scale s
double CImageSpace::getValue(const Vector2D &p,double s) {
   // Number of standard deviations away that we evaluate.  Let this be 3 or 4
   const double nSD = 3.0;
   int x,y;
   
   // First determine the range of the kernel, the starting pixel and ending pixel
   double range = s*nSD;
   double range2 = range*range;

   int xStart = (int) (p.x - range);
   int xEnd = (int) (p.x + range);
   int yStart = (int) (p.y - range + 0.5);
   int yEnd = (int) (p.y + range - 0.5);

   // Adjust the ranges to the image plane
   xStart = xStart < 0 ? 0 : xStart;
   xEnd = xEnd >= w ? w-1 : xEnd;
   yStart = yStart < 0 ? 0 : yStart;
   yEnd = yEnd >= w ? w-1 : yEnd;

   // Precompute the x-exponent for each column.  First, dx is the distance from the center of kernel
   // to center of each pixel, scaled so that dx*dx = dist^2/2*s^2, 
   double dScale = 1.0 / (sqrt(2.0)*s);
   double dxScaled = ((xStart+0.5)-p.x) * dScale;  
   for(x=xStart;x<=xEnd;x+=1) {
      xExpArray[x] = exp(-dxScaled*dxScaled);
      dxScaled += dScale;
   }
   
   // Now we need to compute scaled dy, because inside this loop we will be computing 
   // the y-component of the exponent
   double dy = (yStart+0.5) - p.y;  
   double dyScaled = dy * dScale;

   // The sum of all kernel weights, should add up to 1/(2*pi*s^2), but may not if the kernel falls off the
   // image
   // double wSum = 0.0;

   // The kernel value
   double kSum = 0.0;
   
   // Ok, go through each row and compute the kernel
   for(y=yStart;y<=yEnd;y+=1) {
      double yExp = exp(-dyScaled*dyScaled);
      double wght;

      // Now compute the range of x that we will be evaluating.  We are restricting the kernel to a circle
      // so the range of x depends on dy
      double xRange = sqrt(range2 - dy*dy);
      xStart = (int) (p.x - xRange + 0.5);
      xEnd = (int) (p.x + xRange - 0.5);
      xStart = xStart < 0 ? 0 : xStart;
      xEnd = xEnd >= w ? w-1 : xEnd;

      // Get a pointer to the image data
      double *iRow = image.row(y);

      // OK, now, we can run through the x array and comput everything
      for(x = xStart;x<=xEnd;x+=1) {
         wght = yExp * xExpArray[x];
         // wSum += wght;
         kSum += wght*iRow[x];
      }

      dy += 1.0;
      dyScaled += dScale;
   }

   // OK, we are now almost finished.  Divide the kSum by wSum.  This way, the kernel response is independant of the
   // area of the kernel.
   double area = 2*M_PI*s*s;
   return kSum / area;
}

// Interpolate the scale space of image gradient
Vector2D CImageSpace::getGradient(const Vector2D &p,double s) {
   // Number of standard deviations away that we evaluate.  Let this be 3 or 4
   const double nSD = 3.0;
   int x,y;
   
   // First determine the range of the kernel, the starting pixel and ending pixel
   double range = s*nSD;
   double range2 = range*range;

   int xStart = (int) (p.x - range);
   int xEnd = (int) (p.x + range);
   int yStart = (int) (p.y - range + 0.5);
   int yEnd = (int) (p.y + range-0.5);

   // Adjust the ranges to the image plane
   xStart = xStart < 0 ? 0 : xStart;
   xEnd = xEnd >= w ? w-1 : xEnd;
   yStart = yStart < 0 ? 0 : yStart;
   yEnd = yEnd >= w ? w-1 : yEnd;

   // Precompute the x-exponent for each column.  First, dx is the distance from the center of kernel
   // to center of each pixel, scaled so that dx*dx = dist^2/2*s^2, 
   double dScale = 1.0 / (sqrt(2.0)*s);
   double dxScaled = ((xStart+0.5)-p.x) * dScale;  
   for(x=xStart;x<=xEnd;x+=1) {
      xExpArray[x] = exp(-dxScaled*dxScaled);
      dxScaled += dScale;
   }
   
   // Now we need to compute scaled dy, because inside this loop we will be computing 
   // the y-component of the exponent
   double dy = (yStart+0.5) - p.y;  
   double dyScaled = dy * dScale;

   // The sum of all kernel weights, should add up to 1/(2*pi*s^2), but may not if the kernel falls off the
   // image
   // double wSumDX = 0.0,wSumDY = 0.0;

   // The kernel value
   double kSumDX = 0.0,kSumDY = 0.0;
   
   // Ok, go through each row and compute the kernel
   for(y=yStart;y<=yEnd;y+=1) {
      double yExp = exp(-dyScaled*dyScaled);
      double wght;

      // Now compute the range of x that we will be evaluating.  We are restricting the kernel to a circle
      // so the range of x depends on dy
      double xRange = sqrt(range2 - dy*dy);
      xStart = (int) (p.x - xRange + 0.5);
      xEnd = (int) (p.x + xRange - 0.5);
      xStart = xStart < 0 ? 0 : xStart;
      xEnd = xEnd >= w ? w-1 : xEnd;

      // Get a pointer to the image data
      double *iRow = image.row(y);

      // OK, now, we can run through the x array and compute everything
      double dx = ((xStart+0.5) - p.x);  
      for(x = xStart;x<=xEnd;x+=1) {
         wght = yExp * xExpArray[x];
         kSumDX += dx*wght*iRow[x];
         kSumDY += dy*wght*iRow[x];

         // Count of xRange - it's distance to the pixel's center
         dx += 1.0;
      }

      dy += 1.0;
      dyScaled += dScale;
   }

   // OK, we are now almost finished.  Derivative gets scaled by cube of sigma
   double sqrt2Pi = sqrt(M_PI*2);
   double scaleMult = 1.0 / (s*s*s*sqrt2Pi);
   return Vector2D(kSumDX*scaleMult,kSumDY*scaleMult);
}

// Compute the one-jet
double CImageSpace::getOneJet(const Vector2D &p,double s,Vector2D &outGradient) {
   // Number of standard deviations away that we evaluate.  Let this be 3 or 4
   const double nSD = 3.0;
   int x,y;
   
   // First determine the range of the kernel, the starting pixel and ending pixel
   double range = s*nSD;
   double range2 = range*range;

   int xStart = (int) (p.x - range);
   int xEnd = (int) (p.x + range);
   int yStart = (int) (p.y - range + 0.5);
   int yEnd = (int) (p.y + range-0.5);

   // Adjust the ranges to the image plane
   xStart = xStart < 0 ? 0 : xStart;
   xEnd = xEnd >= w ? w-1 : xEnd;
   yStart = yStart < 0 ? 0 : yStart;
   yEnd = yEnd >= w ? w-1 : yEnd;

   // Precompute the x-exponent for each column.  First, dx is the distance from the center of kernel
   // to center of each pixel, scaled so that dx*dx = dist^2/2*s^2, 
   double dScale = 1.0 / (sqrt(2.0)*s);
   double dxScaled = ((xStart+0.5)-p.x) * dScale;  
   for(x=xStart;x<=xEnd;x+=1) {
      xExpArray[x] = exp(-dxScaled*dxScaled);
      dxScaled += dScale;
   }
   
   // Now we need to compute scaled dy, because inside this loop we will be computing 
   // the y-component of the exponent
   double dy = (yStart+0.5) - p.y;  
   double dyScaled = dy * dScale;

   // The sum of all kernel weights, should add up to 1/(2*pi*s^2), but may not if the kernel falls off the
   // image
   // double wSumDX = 0.0,wSumDY = 0.0;

   // The kernel value
   double kSum = 0.0;
   double kSumDX = 0.0,kSumDY = 0.0;
   
   // Ok, go through each row and compute the kernel
   for(y=yStart;y<=yEnd;y+=1) {
      double yExp = exp(-dyScaled*dyScaled);
      double wght;

      // Now compute the range of x that we will be evaluating.  We are restricting the kernel to a circle
      // so the range of x depends on dy
      double xRange = sqrt(range2 - dy*dy);
      xStart = (int) (p.x - xRange + 0.5);
      xEnd = (int) (p.x + xRange - 0.5);
      xStart = xStart < 0 ? 0 : xStart;
      xEnd = xEnd >= w ? w-1 : xEnd;

      // Get a pointer to the image data
      double *iRow = image.row(y);

      // OK, now, we can run through the x array and compute everything
      double dx = ((xStart+0.5) - p.x);  
      for(x = xStart;x<=xEnd;x+=1) {
         wght = yExp * xExpArray[x];
         kSumDX += dx*wght*iRow[x];
         kSumDY += dy*wght*iRow[x];
         kSum += wght*iRow[x];

         // Count of xRange - it's distance to the pixel's center
         dx += 1.0;
      }

      dy += 1.0;
      dyScaled += dScale;
   }

   // OK, we are now almost finished.  Derivative gets scaled by cube of sigma
   double sqrt2Pi = sqrt(M_PI*2);
   double scaleMult1 = 1.0 / (s*s*2*M_PI);
   double scaleMult2 = 1.0 / (s*s*s*sqrt2Pi);
   outGradient.x = kSumDX*scaleMult2;
   outGradient.y = kSumDY*scaleMult2;
   return kSum * scaleMult1;
}

// Compute the two-jet
double CImageSpace::getTwoJet(const Vector2D &p,double s,Vector2D &outGradient,Transform2D &outHessian) {
   return 0;
}

// Compute the Hessian
Transform2D CImageSpace::getHessian(const Vector2D &p,double s) {
   return Transform2D();
}





