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
#include <math.h>
#include "phquintic.h"
#include <BrentLinearMethod.h>
#include <ConjugateGradientMethod.h>

// Begin namespace
NAMESPACE_PAULY_START

using namespace std;

// Windows is screwed up about using the standard C++ libraries
inline cmplx conjg(cmplx &x) {
	cmplx y(x.real(),-x.imag());
	return y;
}

// Pascal's trianlge
int pascalTriangle[9][9] = {
{1,0,0,0,0,0,0,0,0},
{1,1,0,0,0,0,0,0,0},
{1,2,1,0,0,0,0,0,0},
{1,3,3,1,0,0,0,0,0},
{1,4,6,4,1,0,0,0,0},
{1,5,10,10,5,1,0,0,0},
{1,6,15,20,15,6,1,0,0},
{1,7,21,35,35,21,7,1,0},
{1,8,28,56,70,56,28,8,1}};

/***********************************************************************
 * Bezier Curve interpolation
 ***********************************************************************/
BezierCurve::BezierCurve(int degree) {
	this->degree = degree;

   x.resize(degree+1);
	y.resize(degree+1);
	ti.resize(degree+1);
	ti1.resize(degree+1);

   ti[0] = 1;
	ti1[0] = 1;
}

BezierCurve::~BezierCurve() {
}

void BezierCurve::setControlPoint(int number,double x,double y) {
	this->x[number] = x;
	this->y[number] = y;
}

int BezierCurve::nChooseK(int n,int k) {
	if(n < 9)
		return(pascalTriangle[n][k]);
	
	return 0;
}

void BezierCurve::getInterpolation(double t,double &x_out,double &y_out) {
	x_out = y_out = 0;
	int i;
	
	for(i=1;i<=degree;i++) {
		ti[i] = ti[i-1]*t;
		ti1[i] = ti1[i-1]*(1-t);
	}
	
	for(i=0;i<=degree;i++) {
		double mult = nChooseK(degree,i)*ti[i]*ti1[degree-i];
		x_out += x[i]*mult;
		y_out += y[i]*mult;
	}
}	


/***********************************************************************
 * Pythagorean Quintic Computation
 ***********************************************************************/
// Interpolate between the two points
void PHQuintic::interpolate(cmplx d0,cmplx d1) {
    int j;
    // Two rho values
    cmplx rho[2];
    rho[0] = sqrt(d0/d1);
    rho[1] = -rho[0];

    // Four alpha values
    cmplx alpha[4];
    for(j=0;j<2;j++) {
	cmplx discriminant = sqrt(9.0*(1.0+rho[j])*(1.0+rho[j])-4.0*(6.0*rho[j]*rho[j]+2.0*rho[j]+6.0-30.0/d1));
	alpha[j*2] = (3.0+3.0*rho[j]-discriminant)/2.0;
	alpha[j*2+1] = (3.0+3.0*rho[j]+discriminant)/2.0;
    }

    // Eight mu values, four a,b,k,c values
    cmplx mu1[4],mu2[4],A[4],B[4],K[4],C[4],A1[4],A2[4],A3[4],B1[4],B2[4],B3[4];
    double E[4];

    for(j=0;j<4;j++) {
	// Compute curves
	static cmplx i(0.0,1.0);
	cmplx a,b,c,k;
	cmplx a_2,b_2,a_4,b_4;

	mu1[j] = (alpha[j]-sqrt(alpha[j]*alpha[j]-4.0*rho[j/2]))/2.0;
	mu2[j] = (alpha[j]+sqrt(alpha[j]*alpha[j]-4.0*rho[j/2]))/2.0;

      if(mu1[j]==-1.0 || mu2[j]==-1.0) {
         E[j] = 1000000.0;
         continue;
      }

		a = A[j] = mu1[j]/(mu1[j]+1.0);
		b = B[j] = mu2[j]/(mu2[j]+1.0);

		a_2 = a*a;
		b_2 = b*b;
		a_4 = a_2*a_2;

		k = K[j] = d0/(a_2*b_2);
		c = C[j] = k*(a*a_4 - 5.0*b*a_4 + 10.0*a*a_2*b_2) / 30.0;

		// Conjugates
		cmplx a_ = conjg(a);
		cmplx b_ = conjg(b);
		
		// Compute bending energies
		cmplx a3 = A3[j] = i/(8.0*a.imag()*(a-b)*(a-b_));
		cmplx b3 = B3[j] =  i/(8.0*b.imag()*(a-b)*(a_-b));
		cmplx a2 = A2[j] = a3*(1.5*i/a.imag() + 1.0/(a-b) - 3.0/(a-b_));
		cmplx b2 = B2[j] = b3*(1.5*i/b.imag() - 1.0/(a-b) + 3.0/(a_-b));
		cmplx a1 = 1.5*a2*i/a.imag();
		a1 += a3*(0.75/(a.imag()*a.imag())-2.0/((a-b)*(a-b)) +
			6.0/((a - b_)*(a - b_)) - (1.0-2.0*b.imag()/a.imag())/((a-b)*(a-b_)));
		A1[j] = a1;
		
		cmplx b1 = 1.5*b2*i/b.imag();
		b1 += b3*(0.75/(b.imag()*b.imag())-2.0/((a-b)*(a-b)) + 
			6.0/((a_ - b)*(a_ - b)) - (1.0-2.0*a.imag()/b.imag())/((a-b)*(a_-b)));
		B1[j] = b1;
		
		E[j] = getBendingEnergy(1,a,b,k,a1,a2,a3,b1,b2,b3)-getBendingEnergy(0,a,b,k,a1,a2,a3,b1,b2,b3);
	}
	
	// Choose the curve with least bending energy
	int min = 0;
	for(j=1;j<4;j++) {
		if(fabs(E[j]) < fabs(E[min]))
			min = j;
	}
	
	// Remember the parameters of that curve
	a = A[min];b = B[min];c = C[min];k = K[min];
	a1 = A1[min];a2 = A2[min];a3 = A3[min];
	b1 = B1[min];b2 = B2[min];b3 = B3[min];
	this->E = E[min];
	
	// Make a bezier curve
	makeBezier();
}

double PHQuintic::getBendingEnergy(double t) {
	return getBendingEnergy(t,a,b,k,a1,a2,a3,b1,b2,b3);
}

double PHQuintic::getBendingEnergy(double t,cmplx a,cmplx b,cmplx k,
											  cmplx a1,cmplx a2,cmplx a3,
											  cmplx b1,cmplx b2,cmplx b3) 
{
	double e = 2.0*a1.real()*log(abs(t-a)) + 2.0*b1.real()*log(abs(t-b)) -
		2.0*a1.imag()*arg(t-a) - 2.0*b1.imag()*arg(t-b);
	
	cmplx tmp = 2.0*a2/(t-a) + 2.0*b2/(t-b) + a3/((t-a)*(t-a)) + b3/((t-b)*(t-b));
	e -= tmp.real();
	e *= 4.0 / norm(k);
	
	return e;
}

void PHQuintic::getInterpolation(double t,cmplx &r,cmplx &rPrime) {
	
	cmplx t_a = t-a;
	cmplx t_a2 = t_a*t_a;
	cmplx t_a4 = t_a2*t_a2;
	cmplx t_b = t-b;
	cmplx t_b2 = t_b*t_b;
	
	r = k*(t_a4*t_a-5.0*t_a4*t_b+10.0*t_a2*t_a*t_b2)/30.0+c;
	rPrime = k*(t_a2*t_b2);
}

double PHQuintic::getCurvature(double t) {
   
   static double epsilon = 0.0000001;
   double t2 = t+epsilon;

   double x1,y1,x2,y2,nx1,ny1,nx2,ny2;
   getInterpolation(t,x1,y1,nx1,ny1);
   getInterpolation(t2,x2,y2,nx2,ny2);

   double l1 = sqrt(nx1*nx1+ny1*ny1);
   double l2 = sqrt(nx2*nx2+ny2*ny2);
   nx1/=l1;ny1/=l1;
   nx2/=l1;ny2/=l1;


   double ds = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
   double dt = sqrt((nx1-nx2)*(nx1-nx2)+(ny1-ny2)*(ny1-ny2));
   return dt/ds;
   

/*
   // Precompute some terms
   cmplx t_a = t-a;
	cmplx t_a2 = t_a*t_a;
	cmplx t_b = t-b;
	cmplx t_b2 = t_b*t_b;

   // Compute the first derivative
   cmplx rPrime = k*(t_a2*t_b2);

   // Compute the second derivative
   cmplx rSecond = 2.0*k*(t_a*t_b2+t_a2*t_b);

   // Compute curvature
   double kappa = fabs(rPrime.real()*rSecond.imag()-rPrime.imag()*rSecond.real());
	// double kappa = fabs(rPrime.real()*rSecond.real()
	//	+rPrime.imag()*rSecond.imag());

   double denom = sqrt(rPrime.real()*rPrime.real()+rPrime.imag()*rPrime.imag());
   denom = denom*denom*denom;
   kappa /= denom;
   
   // Compute ratio of their norms
   return kappa;
*/
}

void PHQuintic::makeBezier() {
	cmplx p[6];
	cmplx term1 = (a+b-2.0*a*b);
	
	p[0] = 0.0;
	p[1] = k*a*a*b*b;
	p[2] = p[1] - 0.5*k*a*b*term1;
	p[3] = p[2] + (1.0/6.0)*k*(term1*term1 + 2.0*a*(1.0-a)*b*(1.0-b));
	p[4] = p[3] - 0.5*k*(1.0-a)*(1.0-b)*term1;
	p[5] = p[4] + k*(1.0-a)*(1.0-a)*(1.0-b)*(1.0-b);
	
	for(int i=0;i<6;i++) {
		cmplx x = p0 + (p5-p0)*p[i]/5.0;
		bezier.setControlPoint(i,x.real(),x.imag());
	}
}

double PHQuintic::getArcLength(double t0,double t1) {
	return 1;
}

double PHQuintic::getBendingEnergy(double t0,double t1) {
	return getBendingEnergy(t1) - getBendingEnergy(t0);
}

// Get curve at a point (and normals if needed too)
void PHQuintic::getInterpolation(double t,double &xt,double &yt) {
	cmplx r,rPrime;
	getInterpolation(t,r,rPrime);

	// Removed - quintic always on unit interval
   // r = p0 + r*(p5-p0);
	
	xt = r.real();
	yt = r.imag();
}

void PHQuintic::getInterpolation(double t,double &xt,double &yt,double &nxt,double &nyt) {
	cmplx r,rPrime;
	getInterpolation(t,r,rPrime);
	
	// Removed - quintic always on unit interval
   // r = p0 + (p5-p0)*r;
	xt = r.real();
	yt = r.imag();

   // Removed - quintic always on unit interval
	// rPrime = (p5-p0)*rPrime;
	// nxt = rPrime.imag();
	// nyt = rPrime.real();
   nyt = rPrime.imag();
   nxt = rPrime.real();
}

PHQuintic::PHQuintic() : p0(0.0,0.0),p5(1.0,0.0),bezier(5)
{
   T0 = cmplx(1.0,0.0);
   T5 = cmplx(1.0,0.0);
   interpolate(T0,T5);
}

PHQuintic::PHQuintic(double nx1,double ny1,double nx2,double ny2) : 
p0(0.0,0.0),p5(1.0,0.0),
bezier(5)
{
   T0 = cmplx(nx1,ny1);
   T5 = cmplx(nx2,ny2);
   
   interpolate(T0,T5);
}

void PHQuintic::setEndVectors(double nx1,double ny1,double nx2,double ny2) {
   T0 = cmplx(nx1,ny1);
   T5 = cmplx(nx2,ny2);
   
   interpolate(T0,T5);
}

/***********************************************************************
 * An optimization problem for maximizing curvature of an object
 ***********************************************************************/
class PHCurvaturePbm : public Function {
   PHQuintic &phq;
public:
   PHCurvaturePbm(PHQuintic &quintic) : phq(quintic) {};

   double evaluate(const Vector &v) {
      
		if(v(0) < 0.0) {
         return -phq.getCurvature(0.0) + fabs(v(0));
      }
      else if(v(0) > 1.0) {
         return -phq.getCurvature(1.0) + v(0)-1.0;
      }
      else
         return -phq.getCurvature(v(0));
		
		// return -phq.getCurvature(v(0));

   }
};



/***********************************************************************
 * A PHQuinitic with defined end points.
 ***********************************************************************/
CBoundlet::CBoundlet() : vStart(0,0,true), vEnd(0,0,true) {
}

void CBoundlet::compute(const Vector2D &xStart,const Vector2D &xEnd,
                        const Vector2D &nStart,const Vector2D &nEnd) 
{
   

   // We will compute the transform that takes a vector in real space and maps it
   // onto the normalized space.  This vector maps (xEnd-xStart) to (1,0)   
   transform = Transform2D::frameChange(xStart,xEnd-xStart);
   inverseTransform = transform.inverse();

   // Test the inverse
   Transform2D T = transform * inverseTransform;
      
   // Now let's map all the normal vectors to our unit space
   Vector2D vs = inverseTransform * nStart;
   Vector2D ve = inverseTransform * nEnd;

   // Now let's compare the old vectors to the new vectors.  If the difference in one-norm
   // is less than epsilon, the quintic stands
   double dv = (vStart-vs).oneNorm() + (vEnd-ve).oneNorm();

   // OK, store the new vectors and go home
   vStart = vs;
   vEnd = ve;
   
   if(dv > 0.00001) {
      // Now, remember that the user passes in normal vectors, but what we need is tangent vectors,
      // and in the right direction.
      
      vs = Vector2D(vs.y,-vs.x,true);
      ve = Vector2D(ve.y,-ve.x,true);
/*
      if(vs.x < 0) {
         nMult = -1;
         vs = -vs;
         ve = -ve;
      }
      else {
         nMult = 1;
      }
*/
      nMult = 1;
      curve.setEndVectors(vs.x,vs.y,ve.x,ve.y);      
   }

   arcLength = curve.getArcLength(0,1) * xStart.distanceTo(xEnd);
}

// Get the interpolation of the boundary
void CBoundlet::getInterpolation(double t,Vector2D &p,Vector2D &n) {
   Vector2D P(0,0,true),N(0,0,false);
   
   // Get interpolation from ph curve
   curve.getInterpolation(t,P.x,P.y,N.x,N.y);

   // Lets transform these vectors into the real world
   p = transform * P;
   n = (transform * N)*nMult;


}




/***********************************************************************
 * Pythagorean Quintic Test Program
 ***********************************************************************/
/*
void main(void) {
	PHQuintic q(0,0,1,0,60*M_PI/180,-18*M_PI/180);
	for(double t=0;t<=1.0;t+=0.01) {
		double x,y,nx,ny;
		q.getInterpolation(t,x,y,nx,ny);
		printf("%lg\t%lg\t%lg\n",t,x,y);
	}
	printf("BE = %lg\n",q.getBendingEnergy(0,1));
}
*/

// End namespace
NAMESPACE_PAULY_END

