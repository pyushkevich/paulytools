#ifndef _ISPACE_H_
#define _ISPACE_H_

#include <math.h>
#include <vector2d.h>
#include <assert.h>
#include "array2d.h"

#ifdef _MATRIX_BOUNDS_CHECK
#define dassert(a) assert(a);
#else
#define dassert(a) ;
#endif

// ITK Includes and forward references
#include <vnl/vnl_matrix.h>
#include <vcl_complex.h>
namespace itk {
  template <class TPixel, unsigned int VDim> class Image;
};

class ImageSpace {
public:
   virtual int width() = 0;
   virtual int height() = 0;
   virtual double getPixel(int x,int y) = 0;
   
   // Get the value of the image function at a point.  This is the only one you MUST implement
   virtual double getValue(const Vector2D &p,double s) = 0;
   
   // Compute the one-jet of the image.  Subclass may implement this differently
   // for faster computation
   virtual double getOneJet(const Vector2D &p,double s,Vector2D &outGradient) {
      // Deltas
      double delta = 0.001;
      Vector2D dx(delta,0,0),dy(0,delta,0);

      double fp = getValue(p,s);
      double fx = getValue(p+dx,s);
      double fy = getValue(p+dy,s);

      outGradient = Vector2D((fx-fp)/delta,(fy-fp)/delta,false);
      return fp;
   }

   // Compute the two-jet of the image.  This will only set the upper left matrix of the transform, so
   // make sure to pass in a transform that does not have a shift
   virtual double getTwoJet(const Vector2D &p,double s,Vector2D &outGradient,Transform2D &outHessian) {
      // Deltas
      double delta = 0.001;
      double mult = 1.0 / delta;
      Vector2D dx(delta,0,0),dy(0,delta,0);
      Vector2D gp;

      // Compute the function at six places
      double f00 = getValue(p,s);
      double f10 = getValue(p+dx,s);
      double f20 = getValue(p+dx*2,s);
      double f01 = getValue(p+dy,s);
      double f02 = getValue(p+dy*2,s);
      double f11 = getValue(p+dy+dx,s);

      // Compute the gradient
      outGradient.x = (f10-f00) * mult;
      outGradient.y = - (f01-f00) * mult;

      // Compute the Hessian
      outHessian[0][0] = ((f20-f10)*mult - outGradient.x)*mult;
      outHessian[0][1] = ((f11-f01)*mult - outGradient.x)*mult;
      outHessian[0][1] = ((f11-f10)*mult - outGradient.y)*mult;
      outHessian[1][1] = ((f02-f01)*mult - outGradient.y)*mult;

      return f00;
   }

   // Compute the image gradient at a given location in the image  Subclass may implement this differently
   // for faster computation
   virtual Vector2D getGradient(const Vector2D &p,double s) {
      Vector2D G;
      getOneJet(p,s,G);
      return G;
   }

   // Compute the image gradient at a given location in the image  Subclass may implement this differently
   // for faster computation
   virtual Transform2D getHessian(const Vector2D &p,double s) {
      Transform2D T;
      Vector2D G;
      getTwoJet(p,s,G,T);
      return T;
   }
   
   // The power kernel returns log(|GRAD|)*<normalize(GRAD),n>^power
   virtual double applyPowerKernel(const Vector2D &p,const Vector2D &n,double s,int power) {
      // Compute normalized gradient and gradient magnitude
      Vector2D grad = getGradient(p,s);
      double gMag = grad.normalize();

		// If gMag is zero, then there is no image match
		if(gMag == 0.0)
			return 0;

      // Cosine of the angle between n and grad
      double c = n.dotProduct(grad);

      // Take c to power (use this discrete algorithm, pow costs too much).
      double a=1;
      while(power > 0) {
         a = (power & 1) ? a*c : a;
         c = c*c;
         power >>= 1;
      } 

      // Great.  Now take log of the gradient magnitude times the cosine
      return log(gMag+1)*a;
   }

   // Derivative of Gaussian Kernel
   double applyDOGKernel(const Vector2D &p,const Vector2D &n,double s) {
      return getGradient(p,s).dotProduct(n);
   }

   // Get the optimal orientation of the LPP kernel at a given position
   virtual Vector2D getLPPDirection(const Vector2D &p,double s) {
      // We need a Hessian from the image.
      Transform2D H = getHessian(p,s);
      Vector2D n;

      // l is the smaller eigenvalue of the hessian
      double l = 0.5 * (H[0][0]+H[1][1]-sqrt((H[0][0]-H[1][1])*(H[0][0]-H[1][1]) + 4*H[0][1]*H[1][0]));

      // n is the eigenvector corresponding to l
      n.x = H[0][1]/(l-H[0][0]);
      n.y = 1.0;

      // Return a unit vector
      n.normalize();
      return n;
   }
   
   // Apply the LPP kernel at given position, scale and orientation
   virtual double applyLPPKernel(const Vector2D &p,const Vector2D &n,double s) {
      // We need a Hessian from the image.
      Transform2D H = getHessian(p,s);

      return s*s*H[0][0]*n.x*n.x+2*H[0][1]*n.x*n.y+H[1][1]*n.y*n.y;
   }

   // Apply the LPP kernel at given position, scale and at optimal orientation
   virtual double applyLPPKernel(const Vector2D &p,double s) {
      // We need a Hessian from the image.
      Transform2D H = getHessian(p,s);

      return s*s*(H[0][0]+H[1][1]-sqrt((H[0][0]-H[1][1])*(H[0][0]-H[1][1]) + 4*H[0][1]*H[1][0]));
   }

   // Apply the Fxx+Fyy laplacean kernel
   virtual double applyLaplaceanKernel(const Vector2D &p,double s) {
      // We need a Hessian from the image. (ok, so one operation is extra, big deal)
      Transform2D H = getHessian(p,s);
      return s*s*(H[0][0]+H[1][1]);
   }

   virtual ~ImageSpace() {};

	virtual void getImageBytes(Array2D<unsigned char> &bytes) = 0;

	virtual int getImageMax() = 0;
	virtual int getImageMin() = 0;
};

/*****************************************************************
 Scale Space - Precomputes blurred slices and gradients
  ****************************************************************/

class ScaleSpace : public ImageSpace {
private:
	// Some sizes
	enum {nSlices=60,sliceZero=30};

   // The base image
	Array2D<float> base;
   
	// The transform of the image
  typedef vcl_complex<float> ComplexType;
  typedef vnl_matrix<ComplexType> ComplexMatrixType;
  ComplexMatrixType *baseFFT;

   // Number of stacks that have been computed
   static const double logSqrt2;

	// Array of scales
	float sliceScales[nSlices];

   // Stacks of images
   Array2D<float> jets[nSlices][6];

   // Image dimensions
   int w,h;

   // Scaled up image dimensions
   int w2,h2;

	// Image min/max
	double iMin,iMax;

   // Convert slice number to scale
   double sliceToScale(double slice) {
      return exp(logSqrt2*(slice-sliceZero));
   }

   // void computeBlur(int slice);
   // void computeGrad(int slice);
   
   void computeJet(int slice);
   void sampleJet(const Vector2D &p,double s,int start,int last);
   void computeScaleSpace();

	int fSlice;
	int findSlice(float scale) {
		// Search position is static.  This way we remember where we began the search	
		/*static*/ // int fSlice = sliceZero;
		int lo = 0;
		int hi = nSlices-1;
		
		while(1) {
			if(scale < sliceScales[fSlice]) {
				hi = fSlice;
				fSlice = (lo + fSlice) >> 1;
			}
			else if(scale > sliceScales[fSlice+1]) {
				lo = fSlice;
				fSlice = (hi + fSlice) >> 1;
			}
			else {
				return fSlice;
			}
		}
	}


   // This is used for jet computations
   double jetlet[6];
   
public:
  typedef itk::Image<short,2> ImageType;

   ScaleSpace(ImageType *image);
   ~ScaleSpace();

   // Interpolate the scale space of image intensity
   double getValue(const Vector2D &p,double s) {
      sampleJet(p,s,0,0);
      return jetlet[0];
   }

   // Interpolate the scale space of image gradient
   Vector2D getGradient(const Vector2D &p,double s) {
      sampleJet(p,s,1,2);
      return Vector2D(jetlet[1],jetlet[2]);
   }

   // Compute the one-jet
   double getOneJet(const Vector2D &p,double s,Vector2D &outGradient) {
      sampleJet(p,s,0,2);
      outGradient.x = jetlet[1];
      outGradient.y = jetlet[2];
      return jetlet[0];
   }

   // Compute the two-jet
   double getTwoJet(const Vector2D &p,double s,Vector2D &outGradient,Transform2D &outHessian) {
      sampleJet(p,s,0,5);
      outGradient.x = jetlet[1];
      outGradient.y = jetlet[2];
      outHessian[0][0] = jetlet[3];
      outHessian[0][1] = jetlet[4];
      outHessian[1][0] = jetlet[4];
      outHessian[1][1] = jetlet[5];
      return jetlet[0];
   }

   // Compute the Hessian
   Transform2D getHessian(const Vector2D &p,double s) {
      sampleJet(p,s,3,5);
      Transform2D outHessian;
      outHessian[0][0] = jetlet[3];
      outHessian[0][1] = jetlet[4];
      outHessian[1][0] = jetlet[4];
      outHessian[1][1] = jetlet[5];
      return outHessian;
   }

   int width() {
      return w;
   }
   
   int height() {
      return h;
   }

   double getPixel(int x,int y) {
      return base(x,y);
   }

	int getImageMax() { 
		return (int) iMax; 
	}

	int getImageMin() {
		return (int) iMin;
	}

	void getImageBytes(Array2D<unsigned char> &bytes) {
		bytes.resize(base.width(),base.height());
		double scale = 255.0 / (iMax - iMin);
		for(int y=0;y<base.height();y++) {
			for(int x=0;x<base.width();x++) {
				bytes(x,y) = (unsigned char)((base(x,y) - iMin) * scale);
			}
		}
	}
};


/**************************************************************
 Exponent table
class ExponentTable {
private:
   double *table;
   double step;
   double max;
   int size;

public:
   ExponentTable(double step) {
      max = 4.0;

      this->step = step;
      size = (int)(max/step)+1;
      table = new double[size];

      for(int i=0;i<size;i++)
         table[i] = -1.0;
   }

   ~ExponentTable() {
      delete table;
   }

   double get(double x) {
      if(x > 0)
         return exp(x);
      else if(-x >= max)
         return 0.0;
      
      int i = (int) (size * (-x) / max);
      if(table[i] < 0.0) {
         table[i] = exp(-i*max/size);
      }

      return table[i];
   }
};

double appExp(double x);
 *************************************************************/
class CImageSpace : public ImageSpace {
private:
   // The base image, a 2D array
   Array2D<double> image;

   // Dimensions, copied for faster access
   int w,h;

	// Image min/max
	double iMin,iMax;

   // The exponent array for columns in the kernel box.  This is recomputed each time at a new
   // scale, but I didn't want to reallocate the memory each time
   double *xExpArray;
   
public:
   CImageSpace(itk::Image<short,2> *image);
   ~CImageSpace();

   // Get thje value of the image blurred at scale s
   double getValue(const Vector2D &p,double s);

   // Interpolate the scale space of image gradient
   Vector2D getGradient(const Vector2D &p,double s);

   // Compute the one-jet
   double getOneJet(const Vector2D &p,double s,Vector2D &outGradient);

   // Compute the two-jet
   double getTwoJet(const Vector2D &p,double s,Vector2D &outGradient,Transform2D &outHessian);

   // Compute the Hessian
   Transform2D getHessian(const Vector2D &p,double s);

   int width() {
      return w;
   }
   
   int height() {
      return h;
   }

   double getPixel(int x,int y) {
      return image(x,y);
   }

	int getImageMax() { 
		return (int) iMax; 
	}

	int getImageMin() {
		return (int) iMin;
	}

	void getImageBytes(Array2D<unsigned char> &bytes) {
		bytes.resize(image.width(),image.height());
		double scale = 255.0 / (iMax - iMin);
		for(int y=0;y<image.height();y++) {
			for(int x=0;x<image.width();x++) {
				bytes(x,y) = (unsigned char)((image(x,y) - iMin) * scale);
			}
		}
	}
};

#endif

