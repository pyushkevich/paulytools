#ifndef _VECTOR_3D_H_
#define _VECTOR_3D_H_

class Vector3D {
public:
   // Everyone can access the x and y components of this vector, but w is protected
   double x,y,z,w;

   Vector3D() : x(0),y(0),z(0),w(0) {
   }
 
   // Let w=1 for a point and w=0 for a vector
   Vector3D(double inX,double inY,double inZ,double inW=0.0) : x(inX),y(inY),z(inZ),w(inW) {
   }

   // Addition
   Vector3D operator +(const Vector3D &v) const {
      return Vector3D(x+v.x,y+v.y,z+v.z,w+v.w);
   }

   void operator +=(const Vector3D &v) {
      x += v.x;
      y += v.y;
      z += v.z;
      w += v.w;
   }

   // Substraction
   Vector3D operator -(const Vector3D &v) const {
      return Vector3D(x-v.x,y-v.y,z-v.z,w-v.w);
   }

   void operator -=(const Vector3D &v) {
      x -= v.x;
      y -= v.y;
      z -= v.z;
      w -= v.w;
   }

   // Negation
   Vector3D operator -() const {
      return Vector3D(-x,-y,-z,-w);
   }

   // Constant multiplication
   Vector3D operator *(const double k) const {
      return Vector3D(x*k,y*k,z*k,w*k);
   }

   void operator *=(const double k) {
      x *= k;
      y *= k;
      z *= k;
      w *= k;
   }

   // Constant multiplication
   Vector3D operator /(const double k) const {
      return Vector3D(x/k,y/k,z/k,w/k);
   }

   void operator /=(const double k) {
      x /= k;
      y /= k;
      z /= k;
      w /= k;
   }

   // Dot product
   double dotProduct(const Vector3D &v) const {
      return x*v.x + y*v.y + z*v.z;
   }

   // Vector length
   double twoNorm() const {
      return sqrt(x*x+y*y+z*z);
   }

   // One norm
   double oneNorm() const {
      return fabs(x) + fabs(y) + fabs(z);
   }

   // Infinity norm
   double infinityNorm() const {
      return max(max(fabs(x),fabs(y)),fabs(z));
      // return (fabs(x) > fabs(y)) ? fabs(x) : fabs(y);
   }

   // Distance to another vector
   double distanceTo(const Vector3D &v) const {
      double dx = v.x - x;
      double dy = v.y - y;
      double dz = v.z - z;

      return sqrt(dx*dx+dy*dy+dz*dz);
   }

   // Normalize the vector, return its length before normalization
   double normalize() {
      double length = twoNorm();
      x /= length;
      y /= length;
      z /= length;
      w /= length;
      return length;
   }

   // Scale the vector
   Vector3D scaledBy(double sx,double sy,double sz) const {
      return Vector3D(x*sx,y*sy,z*sz,w);
   }

   // Get a vector normal to this vector
   Vector3D crossProduct(const Vector3D &v) const {
      Vector3D r;
      r.x = y*v.z - z*v.y;
      r.y = z*v.x - x*v.z;
      r.z = x*v.y - y*v.x;

      return r;
   }

   friend class Transform3D;
};
/*
class Transform3D {
private:
   double t[4][4];
public:
   
   // Constructor
   Transform3D() {
      for(int i=0;i<4;i++)
         for(int j=0;j<4;j++)
            t[i][j] = (i==j) ? 1.0 : 0.0;
   }

   // Rotation about an axis
   static Transform2D rotation(double angleRadians,Vector3D axis) {
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
   static Transform2D frameChange(const Vector3D &origin,const Vector3D &unit) {
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
   Vector3D operator *(const Vector3D &v) const {
      return Vector3D(
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
*/
#endif _VECTOR_3D_H_