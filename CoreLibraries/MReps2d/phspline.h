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
* phspline.h
*	-----------
* Closed continuous splines consisting of quintics
******************************************************************/
#ifndef _PH_SPLINE_H_
#define _PH_SPLINE_H_


#ifdef _WIN32
#pragma warning(disable:4275)
#endif

#include <complex>
#include <vector>

using namespace std;

#ifndef M_PI
#define M_PI 3.141592653590
#define M_PI_2  1.570796326795
#endif

#include <mylibs.h>
#include "phquintic.h"
#include "bspline.h"
// #include <cimage/include/cimage.h>
#include <vector2d.h>
#include <optima.h>

// Begin namespace
NAMESPACE_PAULY_START

typedef vector<PHQuintic> PHQuinticList;


/**
* A vertex of a polygon - not too useful outside of this code
*/
struct PolygonVertex {
	// Position and normal
	Vector2D x,n;

	PolygonVertex(const Vector2D &inX,const Vector2D &inN) : 
	x(inX),n(inN) {}

};
typedef vector<PolygonVertex> PVList;

typedef vector<CBoundlet> CBList;

/**
* A boundary representation via a polygon (list of points and normals)
* (for image ops, the points are expected to be in a unit square)
*/
class PolygonBoundary {
public:
	PVList vertices;

	// A blank constructor that creates an empty list
	PolygonBoundary() {};

	// Load the boundary representation from a text file (2 or 4 doubles per line),
	// depending on whether the normals are needed
	// (Scale all points by scaleFactor is an available option)
	// Returns negative value of failure
	int loadPoints(char *fname,bool hasNormals,double scaleFactorX=1.0,double scaleFactorY=1.0);

	// Load the boundary representation from a text file (2 or 4 doubles per line),
	// depending on whether the normals are needed
	// (Scale all points by scaleFactor is an available option)
	// Returns negative value of failure
	int savePoints(char *fname,double scaleFactorX=1.0,double scaleFactorY=1.0);

	// Compute normals numerically for points that don't have normals (by approximating
	// the tangent direction)
	void approximateNormals(bool clockwise=true);

	// Add a polygon to the boudnary (sequential construction)
	void addVertex(const Vector2D &x,const Vector2D &n) {
		vertices.insert(vertices.end(),PolygonVertex(x,n));
	}   
};


/**
* An abstract boundary representation
*/
class ContinuousBoundary : public Function {
public:
	// Get the position along the spline at given t (t ranges from 0 to 1) with same
	// point corresponding to t=0 and t=1
	virtual void getInterpolation(double t,Vector2D &x) = 0;
	virtual void getInterpolation(double t,Vector2D &x,Vector2D &n) = 0;

	virtual void getInterpolationJet(double t,Vector2D &x,Vector2D &x1,Vector2D &x2) = 0;

	// Get maximum range of t
	virtual double tmax() = 0;

	// Construct a spline from a list of point coordinates and normals 
	virtual void buildFromPB(PolygonBoundary &pb,int numSamples) = 0;   

	// Evaluate a function that reflects how medial pair s,t is
	double mediality(double s,double t);

	// Compute the mediality jet (up to 2nd der)
	Vector medialityJet(double s,double t);

	// Evaluate function required to make this a 'Function' subclass
	double evaluate(const Vector &x) {
		double x0 = mediality(x(0),x(1));
		// return x0;
		double x1s = (mediality(x(0)+0.00001,x(1)) - x0) / 0.00001;
		double x1t = (mediality(x(0),fmod(x(1)+0.00001,tmax())) - x0) / 0.00001;

		//if(x1s > x1t || x1s < 0 || x1t < 0)
		//	return 1;
		// else 
		return x0;

		//Vector v = medialityJet(x(0),x(1));
		//return v(1) - x1s;

		// return sqrt(x1s*x1s+x1t*x1t);
		//return sqrt(x1s*x1s+x1t*x1t) - sqrt(v(1)*v(1)+v(2)*v(2));
	}

	// Get medial point formed by two boundary positions
	void getMedialPoint(double s,double t,Vector2D &v,double &r,double &axialAngle,double &objectAngle);

	// Trace the core given a starting search point s,t
	int traceCore(double s,double t,vector<Vector> &out) ;

	// A very rough test to see if the spline contains the circle centered at x,y
	// with radius r.
	bool containsCircle(const Vector2D &x,double r,int samples);

	// Create a polygon boundary that corrensponds to the spline
	void constructPolygonBoundary(PolygonBoundary &pb,int numPolygons); 
};


/**
* A boundary representation using Fourier series
*/
class FourierBoundary : public ContinuousBoundary {
public:
	// Positive and negative coefficients
	vector <cmplx> pos,neg;

	// Scaling/Translation matrix
	Transform2D T;

	// Constructor
	FourierBoundary();

	// Get the position along the spline at given t (t ranges from 0 to 1) with same
	// point corresponding to t=0 and t=1
	void getInterpolation(double t,Vector2D &x);
	void getInterpolation(double t,Vector2D &x,Vector2D &n);
	void getInterpolationJet(double t,Vector2D &x,Vector2D &x1,Vector2D &x2);

	// Get maximum range of t
	double tmax() {
		return 1.0;
	}

	// Load the boundary representation from a text file 
	// that contains fourier coefficients in a list 
	// (first line is the order of highest coefficient, and the consequent lines
	// are coefficients of x,y, alternating, of order 1,-1,2,-2,...
	int load(char *fname,int format);
	int save(char *fname,int format);

	// Formats
	const static int ONEPERLINE;
	const static int MATHEMATICA;
	const static int FOURPERLINE;

	// Set the translation and scaling factors
	void setTransform(const Vector2D &dx,const Vector2D &sx);

	// Construct a spline from a list of point coordinates and normals 
	void buildFromPB(PolygonBoundary &pb,int numSamples);

	// This method returns a vector containing the coefficients of the fourier representation
	Vector getFourierCoeff() const;

	// This loads the fourier boundary from a vector containing the representation
	void setFourierCoeff(const Vector &coeff);
};


/**
* A continous closed spline consisting of a number of PHSplines
*/
class PHSpline : public ContinuousBoundary {
public:
	CBList segments;

	// A blank constructor
	PHSpline();

	// Construct a spline from a list of point coordinates and normals 
	void buildFromPB(PolygonBoundary &pb,int numSplines);
	void buildFromPB(PolygonBoundary &pb);

	// Get the position along the spline at given t (t ranges from 0 to 1) with same
	// point corresponding to t=0 and t=1
	void getInterpolation(double t,Vector2D &x);
	void getInterpolation(double t,Vector2D &x,Vector2D &n);
	void getInterpolationJet(double t,Vector2D &x,Vector2D &x1,Vector2D &x2);

	// Get a list of points of optimal curvature.
	// int getCurvatureMaxima(vector<double> &points,vector<double> &curv);

	// Get maximum range of t
	double tmax() {
		return segments.size();
	}

	int load(char *fname);
	int save(char *fname);
};



/**
* A continous closed b-spline
*/
class CBSpline : public ContinuousBoundary {
private:
public:
	// Pointer to the spline object
	BSpline1D *spline;

	// A blank constructor
	CBSpline() {
		spline = NULL;
	}

	// Destructor
	~CBSpline() {
		reset();
	}

	// Restart
	void reset() {
		if(spline) delete spline;
		spline = NULL;
	}

	// Construct a spline from a list of point coordinates and normals 
	void buildFromPB(PolygonBoundary &pb,int numSplines);
	void buildFromPB(PolygonBoundary &pb);

	// Get the position along the spline at given t (t ranges from 0 to 1) with same
	// point corresponding to t=0 and t=1
	void getInterpolation(double t,Vector2D &x);
	void getInterpolation(double t,Vector2D &x,Vector2D &n);
	void getInterpolationJet(double t,Vector2D &x,Vector2D &x1,Vector2D &x2);

	// Get a list of points of optimal curvature.
	// int getCurvatureMaxima(vector<double> &points,vector<double> &curv);

	// Get maximum range of t
	double tmax() {
		return 1.0;
	}

	int load(char *fname);
	int save(char *fname);
};



// End namespace
NAMESPACE_PAULY_END

#endif //_PH_SPLINE_H_
