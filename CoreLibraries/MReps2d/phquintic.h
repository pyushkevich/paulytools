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
 * phquintic.h
 *	-----------
 * Code for 5th order pyhtagorean hodographs
 ******************************************************************/
#ifndef _PH_QUINTIC_H_
#define _PH_QUINTIC_H_

#ifdef WIN32
#pragma warning(disable:4275)
#endif

#include <cmath>
#include <complex>
#include <vector>



#ifndef M_PI
#define M_PI 3.141592653590
#define M_PI_2  1.570796326795
#endif

#include <mylibs.h>
#include <vector2d.h>

// Begin namespace
NAMESPACE_PAULY_START

typedef std::complex<double> cmplx;

// A Bezier curve
class BezierCurve {
private:
	int degree;
  std::vector <double> x,y;
	std::vector <double> ti,ti1;

	int nChooseK(int n,int k);
public:
	BezierCurve(int degree);
	~BezierCurve();

	void setControlPoint(int number,double x,double y);
	void getInterpolation(double t,double &x_out,double &y_out);
};

// All points are in the complex plane
class PHQuintic {
private:
	// Starting and ending points
	cmplx p0,p5;

	// Starting and ending tangent vectors
	cmplx T0,T5;

	// All the things that parameterize our curve
	cmplx a,b,c,k,a1,a2,a3,b1,b2,b3;

	// Total bending energy of the curve
	double E;

	// Interpolate between the two points
	void interpolate(cmplx d0,cmplx d1);

	double getBendingEnergy(double t);

	double getBendingEnergy(double t,cmplx a,cmplx b,cmplx k,
									cmplx a1,cmplx a2,cmplx a3,
									cmplx b1,cmplx b2,cmplx b3);

	void getInterpolation(double t,cmplx &r,cmplx &rPrime);
	void makeBezier();

public:
	// Bezier curve associated with this curve
	BezierCurve bezier;

	double getArcLength(double t0,double t1);
	double getBendingEnergy(double t0,double t1);

	// Get curve at a point (and normals if needed too)
	void getInterpolation(double t,double &xt,double &yt);
	void getInterpolation(double t,double &xt,double &yt,double &nxt,double &nyt);

   // Compute curvature at a point on the quintic.
   double getCurvature(double t);

	// Create a new curve without specifying any information.  This creates a line between the
   // two points.
	PHQuintic();

   // Create a new curve that has endpoints at 0 and 1 on the x axis with provided
   // normal vectors.  This is the new way to operate with PHQuinitcs that allows
   // curvature control, and should be used along with the CBoundlet class.
   PHQuintic(double nx1,double ny1,double nx2,double ny2);

   // Set endpoint vectors of the curve
   void setEndVectors(double nx1,double ny1,double nx2,double ny2);
};


class CBoundlet {
private:
   // This describes the curve
   PHQuintic curve;
   
   // Vectors that describe normals of the ends of the boundlet when it
   // is mapped onto a unit interval.
   Vector2D vStart,vEnd;

   // A transform that maps the boundlet into real space
   Transform2D transform,inverseTransform;

   // This multiplyer is applied to the normals
   double nMult;

   // The arc length
   double arcLength;

public:
   // Create a dummy boundlet
   CBoundlet();
	virtual ~CBoundlet() {};

   // Get the interpolation of the boundary.  The vectors that are input as references will
   // be returned as 3 element vectors with 1 in the third position for vector x and
   // 0 in the third position for vector n.
   void getInterpolation(double t,Vector2D &p,Vector2D &n);

   // Compute the boundlet given absolute position of the two boundary
   // sites.  The normal vector direction is the desired normal to the
   // boundary and the length is the desired radius of curvature
   void compute(const Vector2D &xStart,const Vector2D &xEnd,const Vector2D &nStart,const Vector2D &nEnd);

   // Return the arc length of this piece of curve.
   double getArcLength() {
      return arcLength;
   }
};


// End namespace
NAMESPACE_PAULY_END

#endif //_PH_QUINTIC_H_
