/************************************************************************
 * COMP 257 Final Project
 * Differential Geometry of Implicit Surfaces
 * Author:	Paul Yushkevich
 * Module:	blobmodl.cpp	
 * Last Update: Dec 13, 1998
 *
 * Description:
 * Implementation for classes dealing with implicit surfaces, 
 * curvature computation, models, ovoids and singular surfaces.
 *************************************************************************/
#include "blobmodl.h"
#include <FL/glut.h>
#include <FL/fl_draw.h>
#include "glcode.h"
#include <vnl/vnl_inverse.h>

extern Model *model;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/************************************************************************
 * Vector and matrix functions.  I have my own C++ library for that, but 
 * these are faster since they are fixed size
 ************************************************************************/
inline double kNorm(double x,double y,double k) {
	return pow(pow(fabs(x),k)+pow(fabs(y),k),1/k);
}

inline double dotProduct(double u[3],double v[3]) {
	return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
}

inline void crossProduct(double u[3],double v[3],double uxv[3]) {
	uxv[0] = u[1]*v[2] - u[2]*v[1];
	uxv[1] = u[2]*v[0] - u[0]*v[2];
	uxv[2] = u[0]*v[1] - u[1]*v[0];
}


void multiplyMatrixByVector(double H[3][3],double u[3],double Hu[3]) {
	for(int i=0;i<3;i++) {
		Hu[i]=0;
		for(int j=0;j<3;j++) {
			Hu[i] += H[i][j]*u[j];
		}
	}
}

double multiplyMatrixByTwoVectors(double H[3][3],double u[3],double v[3]) {
	double Hu[3];
	multiplyMatrixByVector(H,u,Hu);
	return dotProduct(Hu,v);
}

inline void scaleVector(double u[3],double s,double su[3]) {
	su[0] = s*u[0];
	su[1] = s*u[1];
	su[2] = s*u[2];
}

inline void vectorSum(double u[3],double v[3],double uplusv[3]) {
	uplusv[0] = u[0]+v[0];
	uplusv[1] = u[1]+v[1];
	uplusv[2] = u[2]+v[2];
}

inline double normVector(double n[3]) {
	double length = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	n[0] /= length;
	n[1] /= length;
	n[2] /= length;
	return length;
}

/************************************************************************
 * This function computes local shape properties
 * by method outlined in O'Neil V.4 
 ************************************************************************/
void PointData::compute(ISurface *surf,bool fastKappa) {
	// V,W are vectors in the tangent plane
	double V[3],W[3];
	
	// Compute the gradient and hessian
	surf->computeJet(x[0],x[1],x[2],Di,Dij);				

	if(Di[0]==0 && Di[1]==0) {
		// A degenerate case when gradient is in the Z direction
		V[0] = 1;V[1] = 0;V[2] = 0;
	}
	else {
		// Pick a vector in the tangent plane
		V[0] = -Di[1];
		V[1] = Di[0];
		V[2] = 0;
	}

	// W is the cross product of the gradient with tangent vector 
	crossProduct(Di,V,W);	

	// Find length of G
	double Gmag2 = dotProduct(Di,Di);
	gradMag = sqrt(Gmag2);

	// Scale V and W in such a way that V x W = Di
	double Vmag = sqrt(dotProduct(V,V));
	double Wmag = sqrt(dotProduct(W,W));
	scaleVector(V,gradMag/Vmag,V);
	scaleVector(W,1.0/Wmag,W);

	// Compute covariant derivatives of the gradient in V,W directions
	// as product of Hessian with V,W
	double NablaVofG[3],NablaWofG[3];
	multiplyMatrixByVector(Dij,V,NablaVofG);
	multiplyMatrixByVector(Dij,W,NablaWofG);

	// Compute K,H according to formulae on p. 217
	double NablaV_x_NablaW[3],NablaV_x_W[3],V_x_NablaW[3],sum[3];
	crossProduct(NablaVofG,NablaWofG,NablaV_x_NablaW);
	crossProduct(NablaVofG,W,NablaV_x_W);
	crossProduct(V,NablaWofG,V_x_NablaW);
	vectorSum(NablaV_x_W,V_x_NablaW,sum);

	// Gaussian curvature
	K = dotProduct(Di,NablaV_x_NablaW)/(Gmag2*Gmag2);

	// I drop the minus because my my normals are pointing the opposite of what they
	// should.
	H = dotProduct(Di,sum)/(2*Gmag2*gradMag);

	// Compute the kappas from K and H.  Absolute value prevents numerical error
	// artifacts.
	kappa1 = H+sqrt(fabs(H*H-K));
	kappa2 = H-sqrt(fabs(H*H-K));

	// Find the principal directions unless the user wants fast kappa computation
	if(fastKappa)
		return;

	// The rest of computation assume that length of V is one.  That scales covariant
	// derivative in V direction scale by same amount.  That's why we divide a,b by
	// powers of gradient magnitude
	double a = dotProduct(V,NablaVofG) / Gmag2;
	double b = dotProduct(W,NablaVofG) / gradMag;
	double c = dotProduct(W,NablaWofG);

	// Find the angle separating the principal directions from V,W.  This computation
	// follows from Excercise 22 on p. 223
	if(a==c && b==0) {
		// Umbillic - I set f0,f1 to zeros to denote the fact.
		f[0][0] = f[0][1] = f[0][2] = 0;
		f[1][0] = f[1][1] = f[1][2] = 0;
		umbillic = true;
	}
	else {
		// Theta is the angle between V and principal direction 1
		double theta = atan(2*b/(a-c));		

		double t1[3],t2[3];
		scaleVector(V,cos(theta/2),t1);
		scaleVector(W,sin(theta/2),t2);
		vectorSum(t1,t2,f[0]);

		scaleVector(V,-sin(theta/2),t1);
		scaleVector(W,cos(theta/2),t2);
		vectorSum(t1,t2,f[1]);

		umbillic = false;

		// There is also an asymptotic direction angle
		if(K < 0)
			thetaAsymptote = atan(sqrt(-kappa1/kappa2));
	}

	// Don't forger the last frame element
	scaleVector(Di,1.0/gradMag,f[2]);
}

void PointData::getVectorForAngle(double theta,double vector[3]) {
	static double t1[3],t2[3];
	scaleVector(f[0],cos(theta),t1);
	scaleVector(f[1],sin(theta),t2);
	vectorSum(t1,t2,vector);
}


/************************************************************************
 * ISurface Attributes and Methods
 ************************************************************************/
// Currently rendering ISurface instance
ISurface *ISurface::current;

// Used by polygonize to get function values
double ISurface::function(double x,double y,double z) {
	return ISurface::current->getFunction(x,y,z);
}

/**
 * Called by polygonize for every triangle in the mesh
 * stores triangles and points in internal data structures
 */
int ISurface::polygonCallback(int v1,int v2,int v3,IMP_VERTICES vv) {
	int vindex[3] = {v1,v2,v3};

	for(int i=0;i<3;i++) {
		int v = vindex[i];
		if(ISurface::current->points.get(v)==NULL) {
			PointData *p = new PointData(vv.ptr[v].position.x,
												  vv.ptr[v].position.y,
												  vv.ptr[v].position.z);
			p->f[2][0] = vv.ptr[v].normal.x;
			p->f[2][1] = vv.ptr[v].normal.y;
			p->f[2][2] = vv.ptr[v].normal.z;
			ISurface::current->points.set(v,p);

			// Compute the bounding cube
			double b = ISurface::current->boundingCubeSize;
			b = (fabs(p->x[0])>b) ? fabs(p->x[0]) : b;
			b = (fabs(p->x[1])>b) ? fabs(p->x[1]) : b;
			b = (fabs(p->x[2])>b) ? fabs(p->x[2]) : b;
			ISurface::current->boundingCubeSize = b;
		}
	}

	Triangle *t = new Triangle((PointData *)ISurface::current->points.get(v1),
										(PointData *)ISurface::current->points.get(v2),
										(PointData *)ISurface::current->points.get(v3));

	ISurface::current->triangles.append(t);

	return 1;
}

/**
 * ISurface initializer
 */
ISurface::ISurface() : R(0.0), Rinv(0.0), points(1000,1000), triangles(1000,1000)
{
	setAffineTransform(1,1,1,0,0,0,0,0,0);
	dl = glGenLists(1);
	setColor(1,0,0,1);
	solid = true;
	tetraSize = 0.1;

	computed = false;
	displayed = true;
	solid = false;
}
	
ISurface::~ISurface() {
	destroyPointData();
	glDeleteLists(dl,1);
}

/**
 * Set the translation, rotation in degrees about axis and translation
 */
void ISurface::setAffineTransform(double sx,double sy,double sz,
											 double tx,double ty,double tz,
											 double rx,double ry,double rz)
											 
{
	Sx=sx;Sy=sy;Sz=sz;
	Tx=tx;Ty=ty;Tz=tz;
	Rx=rx;Ry=ry;Rz=rz;
	
	Matrix4x4 T(0.0), S(0.0), Rx(0.0), Ry(0.0), Rz(0.0);
	double s,c;
	
	T.fill_diagonal(1.0);
	T(0,3) = tx;
	T(1,3) = ty;
	T(2,3) = tz;
	
	S.fill_diagonal(1.0);
	S(0,0) = sx;
	S(1,1) = sy;
	S(2,2) = sz;
	
	Rx.fill_diagonal(1.0);
	s = sin(M_PI*rx/180);
	c = cos(M_PI*rx/180);
	Rx(1,1) = c;
	Rx(2,2) = c;
	Rx(1,2) = -s;
	Rx(2,1) = s;
	
	Ry.fill_diagonal(1.0);
	s = sin(M_PI*ry/180);
	c = cos(M_PI*ry/180);
	Ry(0,0) = c;
	Ry(2,2) = c;
	Ry(0,2) = -s;
	Ry(2,0) = s;
	
	Rz.fill_diagonal(1.0);
	s = sin(M_PI*rz/180);
	c = cos(M_PI*rz/180);
	Rz(0,0) = c;
	Rz(1,1) = c;
	Rz(0,1) = -s;
	Rz(1,0) = s;
	
	R = T*Rx*Ry*Rz*S;
  Rinv = vnl_inverse(R);
}

/**
 * Junk the data used for storing triangles
 */
void ISurface::destroyPointData() {
  int v;
  for(v=0;v<points.getSize();v++) 
    {
    PointData *p = (PointData *)points.get(v);
    if(p)
      delete p;
    points.set(v,NULL);
    }

  for(v=0;v<triangles.getSize();v++) 
    {
    Triangle *t = (Triangle *)triangles.get(v);
    if(t)
      delete t;
    triangles.set(v,NULL);
    }

  points.clear();
  triangles.clear();
}


/**
 * Recompute the mesh representation of the surface
 */
void ISurface::recompute(bool buildList) {
	fl_cursor(FL_CURSOR_WAIT);

	// Set the curent ovoid pointer to myself
	current = this;
	
	destroyPointData();

	boundingCubeSize = 0;

	polygonize(ISurface::function,
				  tetraSize,
				  (int)(getBoundCubeSize()/tetraSize),
				  0,0,0,
				  ISurface::polygonCallback,
				  TET); 

	if(buildList)
		buildDisplayList();

	fl_cursor(FL_CURSOR_DEFAULT);
	
	computed = true;
}

/**
 * Build GL list formed by triangles
 */
void ISurface::buildDisplayList() {

	glNewList(dl,GL_COMPILE);
	
	//glEnable(GL_CULL_FACE);
	//glCullFace(GL_BACK);
	
	glBegin(GL_TRIANGLES);
	
	// Create a bunch of normals and vertices
	for(int i=0;i<triangles.getSize();i++) {
		Triangle *t = (Triangle *)triangles.get(i);
		
		glNormal3d(t->p1->f[2][0],t->p1->f[2][1],t->p1->f[2][2]);
		glVertex3d(t->p1->x[0],t->p1->x[1],t->p1->x[2]);
		glNormal3d(t->p2->f[2][0],t->p2->f[2][1],t->p2->f[2][2]);
		glVertex3d(t->p2->x[0],t->p2->x[1],t->p2->x[2]);
		glNormal3d(t->p3->f[2][0],t->p3->f[2][1],t->p3->f[2][2]);
		glVertex3d(t->p3->x[0],t->p3->x[1],t->p3->x[2]);
	}
	
	glEnd();
	
	glEndList();

	
}

/**
 * Called each time GL displays
 */
void ISurface::display() {
	if(!computed || !displayed)
		return;

	// Set up the transformation
	glPushMatrix();
	glMultMatrixd( R.data_block() );

	if(solid) {
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	} 
	else {
		glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
	}

	glColor4fv(color);
	glCallList(dl);
	glPopMatrix();
}

/**
 * Set material color for this surface
 */
void ISurface::setColor(float r,float g,float b,float a) 
{
	color[0]=r;
	color[1]=g;
	color[2]=b;
	color[3]=a;
	if(computed & displayed)
		buildDisplayList();
}

const double eps = 0.0001;

/**
 * Compute derivative numerically in u direction of a function, given value
 * of the function at x
 */
double ISurface::computeFu(double x[3],double u[3],double F) {
	double rtn = (getFunction(x[0]+eps*u[0],x[1]+eps*u[1],x[2]+eps*u[2]) - F) / eps;
	return rtn;
}

/**
 * Compute Fuv given partial first derivatives and function value
 * of the function at x
 */
double ISurface::computeFuv(double x[3],double u[3],double v[3],double F,double Fu,double Fv) {
	double fuv = getFunction(x[0]+eps*u[0]+eps*v[0],
								    x[1]+eps*u[1]+eps*v[1],
									 x[2]+eps*u[2]+eps*v[2]);
	double rtn = ((fuv - F) / eps - (Fu + Fv)) / eps;
	return rtn;
}

/**
 * Compute Fuu given partial first derivatives and function value
 * of the function at x
 */
double ISurface::computeFuu(double x[3],double u[3],double F,double Fu) {
	double rtn = ((getFunction(x[0]+2*eps*u[0],x[1]+2*eps*u[1],x[2]+2*eps*u[2]) - F) / eps - 2*Fu)/eps;
	return rtn;
}

/**
 * Compute function value, gradient and Hessian at a point numerically
 */
double ISurface::computeJet(double x,double y,double z,double Di[3],double Dij[3][3]) {
	static double ux[3] = {1,0,0};
	static double uy[3] = {0,1,0};
	static double uz[3] = {0,0,1};
	static double X[3];

	X[0] = x;
	X[1] = y;
	X[2] = z;

	double F;
		
	F = getFunction(x,y,z);
	Di[0] = computeFu(X,ux,F);
	Di[1] = computeFu(X,uy,F);
	Di[2] = computeFu(X,uz,F);

	Dij[0][0] = computeFuu(X,ux,F,Di[0]);
	Dij[1][1] = computeFuu(X,uy,F,Di[1]);
	Dij[2][2] = computeFuu(X,uz,F,Di[2]);

	Dij[0][1] = Dij[1][0] = computeFuv(X,ux,uy,F,Di[0],Di[1]);
	Dij[0][2] = Dij[2][0] = computeFuv(X,ux,uz,F,Di[0],Di[2]);
	Dij[1][2] = Dij[2][1] = computeFuv(X,uy,uz,F,Di[1],Di[2]);

	return F;
}

/**
 * Write surface to a registry object
 */
void ISurface::write(Registry *reg) {
	reg->setDoubleValue("color.red",color[0]);
	reg->setDoubleValue("color.green",color[1]);
	reg->setDoubleValue("color.blue",color[2]);

	reg->setDoubleValue("translate.x",Tx);
	reg->setDoubleValue("translate.y",Ty);
	reg->setDoubleValue("translate.z",Tz);

	reg->setDoubleValue("rotate.x",Rx);
	reg->setDoubleValue("rotate.y",Ry);
	reg->setDoubleValue("rotate.z",Rz);

	reg->setDoubleValue("scale.x",Sx);
	reg->setDoubleValue("scale.y",Sy);
	reg->setDoubleValue("scale.z",Sz);

	reg->setDoubleValue("tetrahedronSize",tetraSize);
}


/**
 * Read surface from a registry object
 */
void ISurface::read(Registry *reg) {
	color[0] = reg->getDoubleValue("color.red",1);
	color[1] = reg->getDoubleValue("color.green",1);
	color[2] = reg->getDoubleValue("color.blue",1);

	Tx = reg->getDoubleValue("translate.x",0);
	Ty = reg->getDoubleValue("translate.y",0);
	Tz = reg->getDoubleValue("translate.z",0);

	Rx = reg->getDoubleValue("rotate.x",0);
	Ry = reg->getDoubleValue("rotate.y",0);
	Rz = reg->getDoubleValue("rotate.z",0);

	Sx = reg->getDoubleValue("scale.x",1);
	Sy = reg->getDoubleValue("scale.y",1);
	Sz = reg->getDoubleValue("scale.z",1);

	tetraSize = reg->getDoubleValue("tetrahedronSize",0.2);
	setAffineTransform(Sx,Sy,Sz,Tx,Ty,Tz,Rx,Ry,Rz);

	invalidate();
}



/**
 * Find zero crossing of the surface with line segment
 * returns -1 if no crossing
 */
int ISurface::rootTrap(double x1,double y1,double z1,
				  double x2,double y2,double z2,
				  double &x,double &y,double &z,
				  double treshold) 
{
	double f1 = getFunction(x1,y1,z1);
	double f2 = getFunction(x2,y2,z2);

	if((f1 < 0 && f2 < 0) || (f1 > 0 && f2 > 0)) {
		return -1;
	}

	while(true) {
		x = (x2+x1)/2;
		y = (y2+y1)/2;
		z = (z2+z1)/2;

		double f = getFunction(x,y,z);

		if(fabs(f) < treshold)
			return 0;

		if((f < 0 && f1 < 0) || (f > 0 && f1 > 0)) {
			x1 = x;y1 = y;z1 = z;
			f1 = f;
		}
		else {
			x2 = x;y2 = y;z2 = z;
			f2 = f;
		}
	}
}


/************************************************************************
 * Ovoid Methods
 ************************************************************************/
Ovoid::Ovoid(double alpha,double beta,double strength) :
a(alpha),b(beta),s(strength),ISurface()
{
	tetraSize = 0.1;
	recompute();
	setColor(0.6f,0.6f,0.6f,1.0f);
	solid=true;
}

/**
 * Get the component of the field contributed to by this ovoid
 */ 
inline double Ovoid::getFieldComponent(double x,double y,double z) {
	
	static double xt[3];
	
	for(int i=0;i<3;i++)
		xt[i] = x*Rinv(i,0) + y*Rinv(i,1) + z*Rinv(i,2) + Rinv(i,3);
	
	// Compute the norm function
	double n = kNorm(kNorm(xt[0],xt[1],a),xt[2],b); 
	
	// Return strength / norm.  
	return s/(n);
}

/**
 * Get the value of implicit function defining the ovoid
 */ 
inline double Ovoid::getFunction(double x,double y,double z) {
	
	// Compute the norm function
	double n = kNorm(kNorm(x,y,a),z,b); 
	
	// Surface at 1
	return n - 1;
}

/**
 * Set the value of alpha, beta parameters
 */ 
void Ovoid::setParms(double alpha,double beta) {
	a = alpha;
	b = beta;
	
	recompute();
}

double Ovoid::getAlpha(){
	return a;
}

double Ovoid::getBeta(){
	return b;
}

/**
 * Compute function, gradient and Hessian analytically
 * The expression for the derivative of this function is a bloody mess.
 * It works.
 */ 
double Ovoid::computeJet(double x,double y,double z,double Di[3],double Dij[3][3]) {	
	double g,gx,gy,gz,gxx,gxy,gyy,gyz,gzz,gxz;
	double u,v,w;
	
	u = x*Rinv(0,0) + y*Rinv(0,1) + z*Rinv(0,2) + Rinv(0,3);
	v = x*Rinv(1,0) + y*Rinv(1,1) + z*Rinv(1,2) + Rinv(1,3);
	w = x*Rinv(2,0) + y*Rinv(2,1) + z*Rinv(2,2) + Rinv(2,3);

	// Powers of u,v,w
	double ua_2 = pow(fabs(u),a-2);
	double va_2 = pow(fabs(v),a-2);
	double wb_2 = pow(fabs(w),b-2);
	double ua_1 = ua_2*u;
	double va_1 = va_2*v;
	double wb_1 = wb_2*w;
	double ua = ua_1*u;
	double va = va_1*v;
	double wb = wb_1*w;

	double ux = Rinv(0,0),vx = Rinv(1,0),wx = Rinv(2,0);
	double uy = Rinv(0,1),vy = Rinv(1,1),wy = Rinv(2,1);
	double uz = Rinv(0,2),vz = Rinv(1,2),wz = Rinv(2,2);
	
	double uava = ua+va;
	double T1_2 = pow(uava,b/a-2);
	double T1_1 = T1_2*uava;
	double T1 = T1_1*uava;

	double uvw = T1+wb;
	double T2_2 = pow(T1+wb,1/b-2);
	double T2_1 = T2_2*uvw;
	double T2 = T2_1*uvw;

	double T1x = (T1_1)*(ua_1*ux + va_1*vx) + wb_1*wx;
	double T1y = (T1_1)*(ua_1*uy + va_1*vy) + wb_1*wy;
	double T1z = (T1_1)*(ua_1*uz + va_1*vz) + wb_1*wz;

	// Function evaluation
	g = T2;

	// First derivatives
	gx = (T2_1) * T1x;
	gy = (T2_1) * T1y;
	gz = (T2_1) * T1z;

	// Second derivatives
	gxy = (1-b)*T2_2*T1x*T1y + T2_1*((b-a)*T1_2*(ua_1*uy+va_1*vy)*(ua_1*ux+va_1*vx) +
			(a-1)*T1_1*(ua_2*ux*uy+va_2*vx*vy)+(b-1)*wb_2*wx*wy);

	gxz = (1-b)*T2_2*T1x*T1z + T2_1*((b-a)*T1_2*(ua_1*uz+va_1*vz)*(ua_1*ux+va_1*vx) +
			(a-1)*T1_1*(ua_2*ux*uz+va_2*vx*vz)+(b-1)*wb_2*wx*wz);

	gyz = (1-b)*T2_2*T1y*T1z + T2_1*((b-a)*T1_2*(ua_1*uz+va_1*vz)*(ua_1*uy+va_1*vy) +
			(a-1)*T1_1*(ua_2*uy*uz+va_2*vy*vz)+(b-1)*wb_2*wy*wz);

	gxx = (1-b)*T2_2*T1x*T1x + T2_1*((b-a)*T1_2*(ua_1*ux+va_1*vx)*(ua_1*ux+va_1*vx) +
			(a-1)*T1_1*(ua_2*ux*ux+va_2*vx*vx)+(b-1)*wb_2*wx*wx);
	
	gyy = (1-b)*T2_2*T1y*T1y + T2_1*((b-a)*T1_2*(ua_1*uy+va_1*vy)*(ua_1*uy+va_1*vy) +
			(a-1)*T1_1*(ua_2*uy*uy+va_2*vy*vy)+(b-1)*wb_2*wy*wy);
	
	gzz = (1-b)*T2_2*T1z*T1z + T2_1*((b-a)*T1_2*(ua_1*uz+va_1*vz)*(ua_1*uz+va_1*vz) +
			(a-1)*T1_1*(ua_2*uz*uz+va_2*vz*vz)+(b-1)*wb_2*wz*wz);
	
	// Function f we return is a reciprocal of g
	double g_2 = 1.0 / (g*g);
	double g_3 = g_2 / g;

	Di[0] = -gx*g_2;
	Di[1] = -gy*g_2;
	Di[2] = -gz*g_2;

	Dij[0][0] = 2*gx*gx*g_3 - gxx*g_2;
	Dij[1][1] = 2*gy*gy*g_3 - gyy*g_2;
	Dij[2][2] = 2*gz*gz*g_3 - gzz*g_2;

	Dij[0][1] = Dij[1][0] = 2*gx*gy*g_3 - gxy*g_2;
	Dij[0][2] = Dij[2][0] = 2*gx*gz*g_3 - gxz*g_2;
	Dij[1][2] = Dij[2][1] = 2*gy*gz*g_3 - gyz*g_2;

	// Return the function value
	return 1.0 / g;
}


/*
0.591056

0.90896
0.37315
0.185881

-0.101161       -0.0695591      -0.0346296
-0.0695591      0.0398442       -0.0142082
-0.0346296      -0.0142082      0.0613135
  */

/**
 * Compute the first derivative. 
 */
void Ovoid::computePartials(double x,double y,double z,double normal[]) {
	double u[3];
	for(int i=0;i<3;i++)
		u[i] = x*Rinv(i,0) + y*Rinv(i,1) + z*Rinv(i,2) + Rinv(i,3);
	
	// Compute the normal of this ovoid
	double ua1 = pow(u[0],a-1);
	double va1 = pow(u[1],a-1);
	double wb1 = pow(u[2],b-1);

	double uava = u[0]*ua1+u[1]*va1;
	double uavaTerm1 = pow(uava,b/a-1);

	double firstTerm = - pow(uava*uavaTerm1 + u[2]*wb1,-1/b-1);

	normal[0] = firstTerm * (uavaTerm1 * (ua1*Rinv(0,0)+va1*Rinv(1,0)) + wb1*Rinv(2,0));
	normal[1] = firstTerm * (uavaTerm1 * (ua1*Rinv(0,1)+va1*Rinv(1,1)) + wb1*Rinv(2,1));
	normal[2] = firstTerm * (uavaTerm1 * (ua1*Rinv(0,2)+va1*Rinv(1,2)) + wb1*Rinv(2,2));

	// Compare to numeric result
	// ((ISurface *)this)->computeNormal(x,y,z,normal,normalize);
}

/**
 * Compute the first derivative. 
 */
void Ovoid::computeNormal(double x,double y,double z,double normal[],bool normalize) {
	// Compute the normal of this ovoid
	double xyTerm = pow(x,a)+pow(y,a);
	double firstTerm = pow(pow(xyTerm,b/a)+pow(z,b),1/b-1);
	double secondTerm = pow(xyTerm,b/a-1);

	normal[0] = -1 * firstTerm * secondTerm * x;
	normal[1] = -1 * firstTerm * secondTerm * y;
	normal[2] = -1 * firstTerm * z;

	if(normalize) {
		normVector(normal);
	}
}

/**
 * Write to Registry object
 */
void Ovoid::write(Registry *reg) {
	reg->setDoubleValue("alpha",a);
	reg->setDoubleValue("beta",b);

	this->ISurface::write(reg);
}

/**
 * Read from Registry object
 */
void Ovoid::read(Registry *reg) {
	a = reg->getDoubleValue("alpha",2);
	b = reg->getDoubleValue("beta",2);

	this->ISurface::read(reg);
	recompute();
}


/************************************************************************
 * Model Methods
 ************************************************************************/
Model::Model() : ISurface() 
{
	tetraSize = 0.1;
	//polygonFunction = isurfCallbackColored;
	ovoids = new PList();

	setColor(1,0.6,0,1);

	mode = plainSurface;
}

Model::~Model() {
	delete ovoids;
}

/**
 * Get implicit function value
 */
inline double Model::getFunction(double x,double y,double z) 
{
	double fSum = 0;
	Ovoid *ovoid;
	
	for(int i=0;i<ovoids->getSize();i++)  {
		ovoid = (Ovoid *)ovoids->get(i);
		fSum += ovoid->getFieldComponent(x,y,z);
	}
	
	return 1.0 - fSum;
}

/**
 * Set the display mode for this model
 */
void Model::setMode(DisplayModes mode) {
	this->mode = mode;
	buildDisplayList();
	glutPostRedisplay();
}

/**
 * Choose color depending on Gaussian curvatuire
 */
void Model::colorCode(PointData *pd) {
	if(pd->kappa1*pd->kappa2 > 0)
		glColor3d(color[0],color[1],color[2]);
	else if (pd->kappa1*pd->kappa2 < 0)
		glColor3d(color[2],color[0],color[1]);
	else
		glColor3d(1,1,1);
}


/**
 * Draw in Gaussian Curvature mode
 */
void Model::drawColorCoded() {
	glNewList(dl,GL_COMPILE);
	
	//glEnable(GL_CULL_FACE);
	//glCullFace(GL_BACK);
	
	glBegin(GL_TRIANGLES);
	
	// Create a bunch of normals and vertices
	for(int i=0;i<triangles.getSize();i++) {
		
		Triangle *t = (Triangle *)triangles.get(i);
		
		colorCode(t->p1);
		glNormal3d(t->p1->f[2][0],t->p1->f[2][1],t->p1->f[2][2]);
		glVertex3d(t->p1->x[0],t->p1->x[1],t->p1->x[2]);

		colorCode(t->p2);
		glNormal3d(t->p2->f[2][0],t->p2->f[2][1],t->p2->f[2][2]);
		glVertex3d(t->p2->x[0],t->p2->x[1],t->p2->x[2]);

		colorCode(t->p2);
		glNormal3d(t->p3->f[2][0],t->p3->f[2][1],t->p3->f[2][2]);
		glVertex3d(t->p3->x[0],t->p3->x[1],t->p3->x[2]);
	}
	
	glEnd();
	
	glEndList();
}

/**
 * Draw in principal field mode
 */
void Model::drawPrincField() {

	glNewList(dl,GL_COMPILE);
	
	//glEnable(GL_CULL_FACE);
	//glCullFace(GL_BACK);
	
	glBegin(GL_LINES);
	
	// Create a bunch of normals and vertices
	for(int i=0;i<points.getSize();i++) {
		
		PointData *p = (PointData *)points.get(i);

		double vm = tetraSize/3;
		
		glNormal3d(p->f[2][0],p->f[2][1],p->f[2][2]);
		
		if(!p->umbillic) {	
			glColor3d(1,0,0);
			glVertex3d(p->x[0],p->x[1],p->x[2]);
			glVertex3d(p->x[0]+vm*p->f[0][0],p->x[1]+vm*p->f[0][1],p->x[2]+vm*p->f[0][2]);

			glColor3d(1,1,0);
			glVertex3d(p->x[0],p->x[1],p->x[2]);
			glVertex3d(p->x[0]+vm*p->f[1][0],p->x[1]+vm*p->f[1][1],p->x[2]+vm*p->f[1][2]);
		}

		glColor3d(0,0,1);
		glVertex3d(p->x[0],p->x[1],p->x[2]);
		glVertex3d(p->x[0]+vm*p->f[2][0],p->x[1]+vm*p->f[2][1],p->x[2]+vm*p->f[2][2]);
	}
	
	glEnd();
	
	glEndList();
}


/**
 * Draw in asymptotic field mode
 */
void Model::drawAsymptoticField() {

	glNewList(dl,GL_COMPILE);
	
	//glEnable(GL_CULL_FACE);
	//glCullFace(GL_BACK);
	
	glBegin(GL_LINES);
	
	// Create a bunch of normals and vertices
	for(int i=0;i<points.getSize();i++) {
		
		PointData *p = (PointData *)points.get(i);

		double vm = tetraSize/3;

		// Asymptotic directions are just for hyperbolics
		if(p->K > 0)
			continue;
		
		// Find two rotations of principal axis
		double ad1[3],ad2[3];
		p->getVectorForAngle(p->thetaAsymptote,ad1);
		p->getVectorForAngle(-p->thetaAsymptote,ad2);
		
		glNormal3d(p->f[2][0],p->f[2][1],p->f[2][2]);
				
		glColor3d(1,0,0);
		glVertex3d(p->x[0]-vm*ad1[0],p->x[1]-vm*ad1[1],p->x[2]-vm*ad1[2]);
		glVertex3d(p->x[0]+vm*ad1[0],p->x[1]+vm*ad1[1],p->x[2]+vm*ad1[2]);

		glColor3d(1,1,0);
		glVertex3d(p->x[0]-vm*ad2[0],p->x[1]-vm*ad2[1],p->x[2]-vm*ad2[2]);
		glVertex3d(p->x[0]+vm*ad2[0],p->x[1]+vm*ad2[1],p->x[2]+vm*ad2[2]);
	}
	
	glEnd();
	
	glEndList();
}


void Model::buildDisplayList() {
	switch(mode) {
	case plainSurface:
		this->ISurface::buildDisplayList();
		break;

	case gcCodedSurface:
		drawColorCoded();
		break;

	case asymptField:
		drawAsymptoticField();
		break;
		
	case princField:
		drawPrincField();
		break;
	}

	
}

/**
 * My own version of recompute.  Gets mesh from ISurface::recompute
 * and then computes local geometry of each point
 */
void Model::recompute(bool buildList) {
	this->ISurface::recompute(false);

	for(int v=0;v<points.getSize();v++) {
		PointData *pd = (PointData *)points.get(v);
		pd->compute(this);
	}

	if(buildList)
		buildDisplayList();
}

void Model::computeNormal(double x,double y,double z,double normal[],bool normalize) {
	// Compute the normal of each ovoid
	Ovoid *ovoid;
	
	normal[0] = normal[1] = normal[2] = 0.0;

	for(int i=0;i<ovoids->getSize();i++)  {
		double n[3];
		
		ovoid = (Ovoid *)ovoids->get(i);
		ovoid->computePartials(x,y,z,n);
		
		normal[0] += n[0];
		normal[1] += n[1];
		normal[2] += n[2];
	}	

	if(normalize) {
		normVector(normal);
	}
}

/**
 * Compute all derivatives analytically
 */
double Model::computeJet(double x,double y,double z,double Di[3],double Dij[3][3]) {
	// Compute the normal of each ovoid
	Ovoid *ovoid;
	
	double f = 0;
	for(int r=0;r<3;r++) {
		Di[r] = 0.0;
		for(int c=0;c<3;c++) {
			Dij[r][c] = 0.0;
		}
	}

	for(int i=0;i<ovoids->getSize();i++)  {
		double di[3];
		double dij[3][3];
		
		ovoid = (Ovoid *)ovoids->get(i);
		f += ovoid->computeJet(x,y,z,di,dij);
		
		for(int r=0;r<3;r++) {
			Di[r] -= di[r];
			for(int c=0;c<3;c++) {
				Dij[r][c] -= dij[r][c];
			}
		}
	}	

	return 1-f;
}

/**
 * Write to Registry object
 */
void Model::write(Registry *reg) {
	reg->setIntValue("ovoidCount",ovoids->getSize());
	for(int i=0;i<ovoids->getSize();i++) {
		Ovoid *o = (Ovoid *)ovoids->get(i);
		Registry *oreg = &reg->getSubFolder("ovoid[%d]",i);
		o->write(oreg);
	}

	this->ISurface::write(reg);
}

/**
 * Read from Registry object
 */
void Model::read(Registry *reg) {
	ovoids->clear();
	int nOvoids = reg->getIntValue("ovoidCount",0);

	for(int o=0;o<nOvoids;o++) {
		Ovoid *ovoid = new Ovoid(2,2);
		Registry *oreg = &reg->getSubFolder("ovoid[%d]",o);
		ovoid->read(oreg);
		ovoids->append(ovoid);
	}

	this->ISurface::read(reg);
}

/************************************************************************
 * Parabolic Surface methods
 ************************************************************************/
ParabolicSurface::ParabolicSurface(ISurface *s) : 
surf(s),ISurface() 
{
	tetraSize = 0.2;
	setColor(0,0.6,0,1);
}

/**
 * Get implicit function value
 */
double ParabolicSurface::getFunction(double x,double y,double z)
{
	PointData pd;
	pd.x[0] = x;
	pd.x[1] = y;
	pd.x[2] = z;

	// Compute point properties, but only kappa1 and kappa2
	pd.compute(surf,true);

	// This function has zeroes at parabolic points
	return pd.kappa1*pd.kappa2;
}

/************************************************************************
 * Ridge Surface Methods
 ************************************************************************/
RidgeSurface::RidgeSurface(ISurface *s) : 
surf(s),ISurface() 
{
	tetraSize = 0.2;
	setColor(0.6,0,0,1);
}

/**
 * Get implicit function value
 * Big problem - I don't know which kappa corresponds to ehich
 * principal direction, so until I do this doesn't work right
 */
double RidgeSurface::getFunction(double x,double y,double z)
{
	double epsilon = 0.00001;

	PointData p(x,y,z);

	// Compute point properties.
	p.compute(surf);

	// Take a step in principal direction 1
	PointData p1(x+epsilon*p.f[1][0],
					 y+epsilon*p.f[1][1],
					 z+epsilon*p.f[1][2]);
	p1.compute(surf,true);

	// Take a step in principal direction 2
	
	PointData p2(x+epsilon*p.f[0][0],
					 y+epsilon*p.f[0][1],
					 z+epsilon*p.f[0][2]);
	p2.compute(surf,true);
	

	// Compute derivative of kappa1,kappa2
	double dk1 = (p1.kappa1-p.kappa1)/epsilon;
	double dk2 = (p2.kappa2-p.kappa2)/epsilon;

	// This function has zeroes at parabolic points
	return dk1;//*dk2;
}

/**
 * Draw caustic surface
 */
void CausticSurface::buildDisplayList() {
	glNewList(dl,GL_COMPILE);
	
	//glEnable(GL_CULL_FACE);
	//glCullFace(GL_BACK);
	
	glBegin(GL_POINTS);
	
	// Create a bunch of normals and vertices
	for(int i=0;i<surf->triangles.getSize();i++) {
		Triangle *t = (Triangle *)surf->triangles.get(i);
		
		if(fabs(t->p1->kappa1) > 0.1 && fabs(t->p2->kappa1) > 0.1 && fabs(t->p3->kappa1) > 0.1) {
			glColor3d(color[0],color[1],color[2]);
			
			PointData &p = *(t->p1);
			glNormal3d(p.f[2][0] / p.kappa1,p.f[2][1] / p.kappa1,p.f[2][2] / p.kappa1);
			glVertex3d(p.x[0]-p.f[2][0] / p.kappa1,
				p.x[1]-p.f[2][1] / p.kappa1,
				p.x[2]-p.f[2][2] / p.kappa1);
			
			p = *(t->p2);
			glNormal3d(p.f[2][0] / p.kappa1,p.f[2][1] / p.kappa1,p.f[2][2] / p.kappa1);
			glVertex3d(p.x[0]-p.f[2][0] / p.kappa1,
				p.x[1]-p.f[2][1] / p.kappa1,
				p.x[2]-p.f[2][2] / p.kappa1);
			
			p = *(t->p3);
			glNormal3d(p.f[2][0] / p.kappa1,p.f[2][1] / p.kappa1,p.f[2][2] / p.kappa1);
			glVertex3d(p.x[0]-p.f[2][0] / p.kappa1,
				p.x[1]-p.f[2][1] / p.kappa1,
				p.x[2]-p.f[2][2] / p.kappa1);
		}
		
		if(fabs(t->p1->kappa2) > 0.1 && fabs(t->p2->kappa2) > 0.1 && fabs(t->p3->kappa2) > 0.1) {
			glColor3d(color[2],color[0],color[1]);
			
			PointData &p = *(t->p1);
			glNormal3d(p.f[2][0] / p.kappa2,p.f[2][1] / p.kappa2,p.f[2][2] / p.kappa2);
			glVertex3d(p.x[0]-p.f[2][0] / p.kappa2,
				p.x[1]-p.f[2][1] / p.kappa2,
				p.x[2]-p.f[2][2] / p.kappa2);
			
			p = *(t->p2);
			glNormal3d(p.f[2][0] / p.kappa2,p.f[2][1] / p.kappa2,p.f[2][2] / p.kappa2);
			glVertex3d(p.x[0]-p.f[2][0] / p.kappa2,
				p.x[1]-p.f[2][1] / p.kappa2,
				p.x[2]-p.f[2][2] / p.kappa2);
			
			p = *(t->p3);
			glNormal3d(p.f[2][0] / p.kappa2,p.f[2][1] / p.kappa2,p.f[2][2] / p.kappa2);
			glVertex3d(p.x[0]-p.f[2][0] / p.kappa2,
				p.x[1]-p.f[2][1] / p.kappa2,
				p.x[2]-p.f[2][2] / p.kappa2);
		}
	}
	
	glEnd();
	
	glEndList();
}
