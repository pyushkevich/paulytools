/************************************************************************
 * COMP 257 Final Project
 * Differential Geometry of Implicit Surfaces
 * Author:		Paul Yushkevich
 * Module:
 * Last Update: Dec 13, 1998
 *
 * Description:
 *
 *
 *************************************************************************/

#ifndef _BLOB_MODELLER_H_
#define _BLOB_MODELLER_H_

#include <vnl/vnl_matrix_fixed.h>

#include "implicit.h"
#include <Registry/registry.h>
#include "plib.h"

class ISurface;

/************************************************************************
 * A structure we will use for holding local geometry information
 ************************************************************************/
class PointData {
public:
  // Position of the point
  double x[3];

  // Gradient
  double Di[3];

  // Hessian
  double Dij[3][3];

  // Frame vectors
  double f[3][3];

  // Principal curvatures
  double kappa1,kappa2;

  // Gaussian curvature, Mean curvature
  double K,H;

  // Gradient magnitude
  double gradMag;

  // Asymptotic direction angle
  double thetaAsymptote;

  // Certain binary properties
  bool umbillic;

  PointData(){
    x[0] = x[1] = x[2] = 0;
  };

  PointData(double x,double y,double z) {
    this->x[0] = x;
    this->x[1] = y;
    this->x[2] = z;
  }

  // Computes local shape properties, limiting them just to kappas
  // if fastkappa is true
  void compute(ISurface *surf,bool fastKappa=false);

  // Return a vector that's a rotation of f1princ by theta in tangent plane
  void getVectorForAngle(double theta,double vector[3]);
};

/************************************************************************
 * Class for remembering the triangles in the model
 ************************************************************************/
class Triangle {
public:
  PointData *p1,*p2,*p3;

  Triangle (PointData *p1,PointData *p2,PointData *p3) {
    this->p1 = p1;
    this->p2 = p2;
    this->p3 = p3;
  }
};


/**
 * This is a parent class for ovoids, models and various other
 * surfaces.  It maintains a GL display list, can display itself,
 * can be arbitrarily rotated, translated and scaled
 */
class ISurface {
private:
  // Current instance of this class
  static ISurface *current;
  static double function(double x,double y,double z);
  static int polygonCallback(int v1,int v2,int v3,IMP_VERTICES vertices);

  // Destroy the list of triangles and points
  void destroyPointData();

protected:
  // Display list for this ovoid's surface
  int dl;

  // Affine transform matrix
  typedef vnl_matrix_fixed<double,4,4> Matrix4x4;
  Matrix4x4 R,Rinv;

  double computeFu(double x[3],double u[3],double F);
  double computeFuv(double x[3],double u[3],double v[3],double F,double Fu,double Fv);
  double computeFuu(double x[3],double u[3],double F,double Fu);

  // Size of tetrahedra used
  double tetraSize;

  // Build the GL display list
  virtual void buildDisplayList();

  // internal function that returns bounding size for current object
  virtual double getBoundCubeSize() {
    return 100;
  }

public:
  virtual double getFunction(double x,double y,double z) = 0;

  ISurface();
  ~ISurface();

  // Scale translation and euler angles of rotation for this ovoid
  // These are not used in computations, but are provided for
  // quering
  double Sx,Sy,Sz,Tx,Ty,Tz,Rx,Ry,Rz;

  // Set the translation, rotation in degrees about axis and translation
  void setAffineTransform(double sx,double sy,double sz,
    double tx,double ty,double tz,
    double rx,double ry,double rz);

  //virtual void computeNormal(double x,double y,double z,double normal[3],bool normalize=true);
  virtual double computeJet(double x,double y,double z,double Di[3],double Dij[3][3]);

  // void computePointData(PointData &data);

  // Find itersection of surface witha line segment,  if no intersection
  // returns -1
  int rootTrap(double x1,double y1,double z1,double x2,double y2,double z2,
    double &x,double &y,double &z,double treshold);

  // gl METHODS

  virtual void recompute(bool buildList = true);

  void display();

  // Color
  float color[4];
  void setColor(float r,float g,float b,float a);

  // These flags describe the current state of the surface
  bool computed;
  bool displayed;
  bool solid;

  void invalidate() {
    computed = false;
  }

  void drawWireframe() {
    solid = false;
  }

  void drawSolid() {
    solid = true;
  }

  void setDisplayed(bool b) {
    displayed = b;
    if(displayed && !computed) 
      recompute();
  }

  void setTetraSize(double ts) {
    tetraSize = ts;
    invalidate();
  }
  double getTetraSize() {
    return tetraSize;
  }

  // Registry routines
  void read(Registry *reg);
  void write(Registry *reg);


  // All the triangles in the surface
  PList triangles;

  // All the points in the surface
  PList points;

  // The bounding cube dimension
  double boundingCubeSize;


};

class Ovoid : public ISurface {
private:
  // Alpha, beta and strength of the ovoid
  double a,b,s;

public:
  // Constructor
  Ovoid(double alpha,double beta,double strength=1);

  // Used by the model
  double getFieldComponent(double x,double y,double z);

  // Used to render the ovoid
  double getFunction(double x,double y,double z);

  // Change alpha,beta
  void setParms(double alpha,double beta);
  double getAlpha(),getBeta();

  // Compute normal analytically
  virtual void computeNormal(double x,double y,double z,double normal[3],bool normalize=true);

  // Same, but use matrix in computation
  void computePartials(double x,double y,double z,double normal[]);

  // Compute function,derivatives and second paritals
  double computeJet(double x,double y,double z,double Di[3],double Dij[3][3]);

  // Registry routines
  void read(Registry *reg);
  void write(Registry *reg);
};


class Model : public ISurface {
protected:
  void buildDisplayList();

  void colorCode(PointData *pd);
  void drawColorCoded();
  void drawPrincField();
  void drawAsymptoticField();

public:
  PList *ovoids;

  // A mode enumerator for telling what mode we are running in
  enum DisplayModes {
    plainSurface,gcCodedSurface,princField,asymptField
  } mode;


  Model();
  ~Model();

  double getFunction(double x,double y,double z);

  // Compute normal analytically
  virtual void computeNormal(double x,double y,double z,double normal[3],bool normalize=true);

  // Compute function,derivatives and second paritals
  double computeJet(double x,double y,double z,double Di[3],double Dij[3][3]);

  // Registry routines
  void read(Registry *reg);
  void write(Registry *reg);

  // My own recompute routine
  void recompute(bool buildList=true);

  // Way to set the mode
  void setMode(DisplayModes mode);

};



class ParabolicSurface : public ISurface {
  ISurface *surf;
public:
  ParabolicSurface(ISurface *s);
  double getFunction(double x,double y,double z);

  double getBoundCubeSize() {
    return surf->boundingCubeSize*1.2;
  }


};

class RidgeSurface : public ISurface {
  ISurface *surf;
public:
  RidgeSurface(ISurface *s);
  double getFunction(double x,double y,double z);

  double getBoundCubeSize() {
    return surf->boundingCubeSize*1.5;
  }



};


class CausticSurface : public ISurface {
  ISurface *surf;
public:
  CausticSurface(ISurface *s) {
    surf = s;
    setColor(0.6,0,0.6,1);
  }

  double getFunction(double x,double y,double z) {
    // Not reallya caustic surface
    return 0;
  }

  void recompute(bool buildList=true) {
    computed = true;
    if(!surf->computed)
      surf->recompute();
    if(buildList)
      buildDisplayList();
  }

  double getBoundCubeSize() {
    return surf->boundingCubeSize*4;
  }

  void buildDisplayList();
};

class TestSurface:public ISurface {
public:
  double getFunction(double x,double y,double z) {
    // Not reallya caustic surface
    return x*x+y*y+z*z - 1;
  }
};


#endif
