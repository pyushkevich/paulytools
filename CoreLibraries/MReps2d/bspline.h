#ifndef _BSPLINE_H_
#define _BSPLINE_H_

#include <SMLVec4f.h>

#include <cmath>
#include <vector>
#include <matrix.h>

#include <vnl/vnl_matrix.h>

using namespace std;

/**
 * Spline knot definition
 * This class can be shared by 1 and 2 dimensional splines
 */
class KnotVector
  {
public:
  /**
   * m: Number of control points - 1
   * z: Number of zero-knots at each end (2-4)
   */
  KnotVector(int m,int z);
  ~KnotVector();

  /**
   * This method evaluates the first three derivatives of the basis functions
   * i: Knot index
   * u: Value of the parameter
   * W: Array of four vectors that will receive the derivatives
   * This function depends on the knot sequence.
   * When there are only two zero knots at ends, the 
   * function is matrix multiplication. 
   */
  void basisJet(int i,float u,MySMLVec4f *W) {
    // void (KnotVector ::* bfnBJ)(int,float,SMLVec4f*) = kv.basisJet;
    (this->*basisJetFN)(i,u,W);
  }

  /**
   * Return the u or v value for a given knot index
   */
  float getParmAtKnot(int knotIndex) {
    return knots[knotIndex];
  }

  /**
   * Get knot value for a given t
   */
  int getKnotAtParm(float u);

  /**
   * Create a knot index sequence for an ordered sequence of u values
   */
  void getKnotIndexSequence(vector<float> &in,vector<int> &out);

private:
  // Number of knots, spline size
  int nk, m;

  // A vector of knots
  float *knots;   

  // Matrix multiplication basis jet function (z=2)
  void basisJetMtx(int i,float u,MySMLVec4f *W);

  // Recursive basis jet function (z<>2)
  void basisJetRcr(int i,float u,MySMLVec4f *W);

  // Basis jet function pointer
  void (KnotVector ::* basisJetFN)(int i,float u,MySMLVec4f *W);

  // Precomputed data needed for matrix basis jet computation
  MySMLVec4f *UB;
  MySMLVec4f *jetMul;

  };

/**
 * 1D Control point
 */
class ControlPoint1D
  {
public:
  float x;
  MySMLVec4f nbr;
  };

/**
 * This represents a one dimensional B-Spline
 */
class BSpline1D
{
public:
  typedef vnl_matrix<double> MatrixType;

  /**
   * m: Number of control points - 1
   * k: Number of dimensions at each point
   * z: Number of zero-knots at each end (2-4)
   */
  BSpline1D(int m,int k,int z=4);

  /**
   * Destroy
   */
  virtual ~BSpline1D();

  /**
   * Spline interpolation
   * i:   Control point number
   * W:   Basis vector computed by basisJet
   * o:   Order of result (x,x',x'',x''')
   * d1:  Starting dimension (x,y,z)
   * d2:  Ending dimension (x,y,z)
   * out: Array of length d2-d1+1
   */
  void interpolatePoint(int i,MySMLVec4f *W,int o,int d1,int d2,float *out);

  /**
   * Set control point
   * i: index of control point
   * d: dimension of the control point (x,y,z, etc)
   */
  virtual void setControl(int i,int d,float val);

  /**
   * Get control point
   * i: index of control point
   * d: dimension of the control point (x,y,z, etc)
   */
  virtual float getControl(int i,int d) {
    return P[d][i].x;
  }

  unsigned int numControls() { return m;}

  /**
   * Fit spline to point data, interpolating the last two points
   * Q: Matrix whose rows are points that must be interpolated
   */
  void fitToPoints(const MatrixType &Q);

  // Knot sequence associated with the spline
  KnotVector kv;

protected:  
  // Spline size, number of dimensions
  int m,k;

  // An array of control points (indexed by dimension, position)
  ControlPoint1D **P;
  };

class SimpleBSplineCurve
  {
public:
  SimpleBSplineCurve(unsigned int nControlPoints, unsigned int nDimensionsPerPoint)
  : spline(nControlPoints-1,nDimensionsPerPoint)
  {
    d = nDimensionsPerPoint;
  }

  // Set a control point
  void SetControlPoint(unsigned int iIndex, unsigned int iDim, double value)
  {
    spline.setControl(iIndex,iDim,value);
  }

  // Get a control point
  double GetControlPoint(unsigned int iIndex, unsigned int iDim)
  {
    return spline.getControl(iIndex,iDim);
  }

  // Get the number of control points
  unsigned int GetNumberOfControlPoints()
  {
    return spline.numControls();
  }

  // Get the size of the control grid
  unsigned int GetInterpolationGridSize() { return grid.size();}

  // Set uniform interpolation grid
  void SetUniformInterpolationGrid(unsigned int nPoints)
  {
    double *values = new double[nPoints];
    double step = 1.0 / (nPoints - 1);
    double val = 0.0;
    for (unsigned int i=0;i<nPoints;i++)
      {
      values[i] = val < 1.0 ? val : 1.0;
      val += step;
      }    
    SetInterpolationGrid(nPoints,values);
    delete values;
  }

  // Set the values where the spline is interpolated
  void SetInterpolationGrid(unsigned int nPoints, const double *values)
  {
    grid.resize(nPoints);
    for (unsigned int i=0;i<nPoints;i++)
      {
      grid[i].u = values[i];
      grid[i].k = spline.kv.getKnotAtParm(grid[i].u);
      spline.kv.basisJet(grid[i].k,grid[i].u,&grid[i].W);
      }
  }

  // Interpolate a grid point
  void InterpolateGrid(unsigned int iGrid, float *valArray)
  {
    spline.interpolatePoint(grid[iGrid].k-3,&grid[iGrid].W,0,0,d-1,valArray);
  }

  // Fit the spline to points
  void FitToPoints(BSpline1D::MatrixType &Q)
  {
    spline.fitToPoints(Q);
  }


private:
  struct Node
    {
    double u;
    unsigned int k;
    MySMLVec4f W;
    };

  BSpline1D spline;
  vector<Node> grid;
  unsigned int d;


  };





#endif
