#ifndef _VECTOR_2D_H_
#define _VECTOR_2D_H_

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#include "matrix.h"

class Vector2D;

class Vector2D {
public:
   // Everyone can access the x and y components of this vector, but w is protected
   double x,y,w;

   Vector2D() : x(0),y(0),w(0) {
   }
 
   // Let w=1 for a point and w=0 for a vector
   Vector2D(double inX,double inY,double inW=0.0) : x(inX),y(inY),w(inW) {
   }

   // Addition
   Vector2D operator +(const Vector2D &v) const {
      return Vector2D(x+v.x,y+v.y,w+v.w);
   }

   void operator +=(const Vector2D &v) {
      x += v.x;
      y += v.y;
      w += v.w;
   }

   // Substraction
   Vector2D operator -(const Vector2D &v) const {
      return Vector2D(x-v.x,y-v.y,w-v.w);
   }

   void operator -=(const Vector2D &v) {
      x -= v.x;
      y -= v.y;
      w -= v.w;
   }

   // Negation
   Vector2D operator -() const {
      return Vector2D(-x,-y,-w);
   }

   // Constant multiplication
   Vector2D operator *(const double k) const {
      return Vector2D(x*k,y*k,w*k);
   }

   void operator *=(const double k) {
      x *= k;
      y *= k;
      w *= k;
   }

   // Constant multiplication
   Vector2D operator /(const double k) const {
      return Vector2D(x/k,y/k,w/k);
   }

   void operator /=(const double k) {
      x /= k;
      y /= k;
      w /= k;
   }

   // Dot product
   double dotProduct(const Vector2D &v) const {
      return x*v.x + y*v.y;
   }

   // Vector length
   double twoNorm() const {
      return sqrt(x*x+y*y);
   }

   // One norm
   double oneNorm() const {
      return fabs(x) + fabs(y);
   }

   // Infinity norm
   double infinityNorm() const {
      return (fabs(x) > fabs(y)) ? fabs(x) : fabs(y);
   }

   // Distance to another vector
   double distanceTo(const Vector2D &v) const {
      double dx = v.x - x;
      double dy = v.y - y;

      return sqrt(dx*dx+dy*dy);
   }

   // Distance to another vector
   double distanceToSqr(const Vector2D &v) const {
      double dx = v.x - x;
      double dy = v.y - y;

      return dx*dx+dy*dy;
   }

   // Normalize the vector, return its length before normalization
   double normalize() {
      double length = twoNorm();
      x /= length;
      y /= length;
      w /= length;
      return length;
   }

   // Scale the vector
   Vector2D scaledBy(double sx,double sy) const {
      return Vector2D(x*sx,y*sy,w);
   }

   // Get a vector normal to this vector
   Vector2D getNormal() const {
      return Vector2D(-y,x,w);
   }

   // Get the slope of the vector, in radians
   double getSlope() const {
      return atan2(y,x);
   }

   // Get the slope of the vector, in radians
   double getSlopeDeg() const {
      return 180 * atan2(y,x) / M_PI;
   }

	// Returns a vector whose x,y,z components are the minimum of components of v1, v2
   static Vector2D minimal(const Vector2D &v1,const Vector2D &v2) {
      Vector2D v;
      v.x = (v1.x < v2.x) ? v1.x : v2.x;
      v.y = (v1.y < v2.y) ? v1.y : v2.y;
      v.w = (v1.w < v2.w) ? v1.w : v2.w;
      return v;
   }


	// Returns a vector whose x,y,z components are the maximum of components of v1, v2
   static Vector2D maximal(const Vector2D &v1,const Vector2D &v2) {
      Vector2D v;
      v.x = (v1.x > v2.x) ? v1.x : v2.x;
      v.y = (v1.y > v2.y) ? v1.y : v2.y;
      v.w = (v1.w > v2.w) ? v1.w : v2.w;
      return v;
   }

   friend class Transform2D;
};

class Transform2D {
private:
   double t[3][3];
public:
   
   // Constructor
   Transform2D() {
      t[0][0]=t[1][1]=t[2][2]=1.0;
      t[0][1]=t[0][2]=t[1][0]=0.0;
      t[1][2]=t[2][0]=t[2][1]=0.0;
   }

   // Construction helpers
   static Transform2D rotation(double angleRadians) {
      Transform2D T;
      T.t[0][0] = T.t[1][1] = cos(angleRadians);
      T.t[0][1] = sin(angleRadians);
      T.t[1][0] = -T.t[0][1];
      return T;
   }

   static Transform2D translation(double dx,double dy) {
      Transform2D T;
      T.t[0][2] = dx;
      T.t[1][2] = dy;
      return T;      
   }

   static Transform2D scaling(double xScale,double yScale) {
      Transform2D T;
      T.t[0][0] = xScale;
      T.t[1][1] = yScale;
      return T;      
   }

   // This creates a transform that moves,rotates and scales the coordinate system.  The 
   // transform of point [0,0,1] is origin and transform of [1,0,1] is origin+unit
   static Transform2D frameChange(const Vector2D &origin,const Vector2D &unit) {
      Transform2D T;
      T.t[0][0] = unit.x;
      T.t[0][1] = -unit.y;
      T.t[1][0] = unit.y;
      T.t[1][1] = unit.x;
      T.t[0][2] = origin.x;
      T.t[1][2] = origin.y;
      return T;      
   }


   // Member access
   double* operator[](int i) {
      return t[i];
   }

   // Transforms can only be stacked together (multiplied)
   Transform2D operator *(const Transform2D &T) const {
      Transform2D S;
      for(int i=0;i<3;i++) {
         for(int j=0;j<3;j++) {
            S.t[i][j] = 0.0;
            for(int k=0;k<3;k++) {
               S.t[i][j] += t[i][k]*T.t[k][j];
            }
         }
      }

      return S;
   }

   // Transforms may be applied to vectors
   Vector2D operator *(const Vector2D &v) const {
      return Vector2D(
         t[0][0]*v.x+t[0][1]*v.y+t[0][2]*v.w,
         t[1][0]*v.x+t[1][1]*v.y+t[1][2]*v.w,
         t[2][0]*v.x+t[2][1]*v.y+t[2][2]*v.w);
   }

   // Transforms have inverses
   Transform2D inverse() const {
      // Inverse of a transform is easy.  Since the bottom row is 0 0 1, we need to compute the 
      // inverse of the inside 2x2 matrix using cramer's rule and the rest is cakewalk.
      Transform2D inv;

      // Zero determinants are bad.  But for standard transforms they do not happen.
      double det = t[0][0]*t[1][1]-t[0][1]*t[1][0];

      inv.t[0][0] = t[1][1]/det;
      inv.t[0][1] = -t[0][1]/det;
      inv.t[1][0] = -t[1][0]/det;
      inv.t[1][1] = t[0][0]/det;

      // Now compute the shifts.
      inv.t[0][2] = (t[0][1]*t[1][2]-t[1][1]*t[0][2]) / det;
      inv.t[1][2] = (t[1][0]*t[0][2]-t[0][0]*t[1][2]) / det;

      return inv;
   }

};

#endif // _VECTOR_2D_H_
