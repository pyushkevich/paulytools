/******************************************************************
 * OPTIMA Library                                                 *
 ******************************************************************
 * Author:						Paul Yushkevich
 *
 * Date:							Apr 1, 1999
 *
 * Description					Multidimensional optimization algorithms
 *									See http://www.cs.unc.edu/~pauly/optima
 *									
 *	Sources:						"Numerical Recepies in C", 
 *									Michaelewitz, "Genetic Algorithms + Data
 *									Structures = Evolutionary Programs"
 *
 * Dependencies:				PY Matrix library, CLAPACK
 ******************************************************************
 * BrentLinearMethod.cpp
 *	---------------------
 * This method from NRC-10.2 allows us to optimize along a vector in
 * n-space.
 ******************************************************************/


#include <math.h>
#include <BrentLinearMethod.h>

// Begin namespace
NAMESPACE_PAULY_START

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define TOL 2.0e-4
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);


inline double SIGN(double a,double b) { return (b<0) ? -a : a;}
inline double FMAX(double a,double b) { return (a<b) ? b : a;}

bool MinBrakRoutine::isFinished() {
	return done || (fb <= fc);
}

void MinBrakRoutine::performIteration() {
	// Skip if finished
	if(isFinished())
		return;

	// Keep returning here until we bracket.
	// Compute u by parabolic extrapolation from a; b; c. TINY is used to 
	// prevent any possible division by zero.
	r=(bx-ax)*(fb-fc); 
	q=(bx-cx)*(fb-fa);
	u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
	ulim=(bx)+GLIMIT*(cx-bx);
	
	// We won't go farther than this. Test various possibilities:
	if ((bx-u)*(u-cx) > 0.0) { 
		// Parabolic u is between b and c: try it.
		fu=func(u);
		if (fu < fc) { 
			// Got a minimum between b and c.
			ax=(bx);
			bx=u;
			fa=(fb);
			fb=fu;
			done = true;
			return;
		} 
		else if (fu > fb) { 
			// Got a minimum between between a and u.
			cx=u;
			fc=fu;
			done = true;
		}
		
		u=cx+GOLD*(cx-bx); 
		// Parabolic t was no use. Use default mag-nication. 
		fu=func(u);
	} else if ((cx-u)*(u-ulim) > 0.0) { 
		//Parabolic t is between c and its allowed limit. 
		fu=func(u);
		if (fu < fc) {
			SHFT(bx,cx,u,cx+GOLD*(cx-bx));
			SHFT(fb,fc,fu,func(u));
		}
	} else if ((u-ulim)*(ulim-cx) >= 0.0) { 
		//Limit parabolic u to maximumallowed value. 
		u=ulim;
		fu=func(u);
	} else { 
		// Reject parabolic u, use default magnica-tion. 
		u=(cx)+GOLD*(cx-bx);
		fu=func(u);
	}
	SHFT(ax,bx,cx,u); // Eliminate oldest point and continue.
	SHFT(fa,fb,fc,fu);
}

inline double MinBrakRoutine::func(double value) {
	// Function must be of one dimension
	vec(0) = value;
	return problem->evaluate(vec);
}

MinBrakRoutine::MinBrakRoutine(Function *problem,double _ax, double _bx) 
{
	this->problem = problem;
	ax = _ax;
	bx = _bx;

	done = false;
	vec.setSize(1);

	fa=func(ax);
	fb=func(bx);
	if (fb > fa) { 
		//Switch roles of a and b so that we can go
		// downhill in the direction from a to b. 
		SHFT(dum,ax,bx,dum);
		SHFT(dum,fb,fa,dum);
	}
	
	// First guess for c.
	cx=(bx)+GOLD*(bx-ax); 
	fc=func(cx);
}

LineMinRoutine::LineMinRoutine(Function *problem) {
	this->problem = problem;
	vec.setSize(1);
}

inline double LineMinRoutine::func(double value) {
	vec(0) = value;
	return problem->evaluate(vec);
}

inline double BrentRoutine::getFAtMinimum() {
	return fx;
}

inline double BrentRoutine::getMinimum() {
	return x;
}

void BrentRoutine::performIteration() {
	if(++iter > ITMAX) {
		// Should never happen
		dassert(0);
		xmin=x; 
		fxmin = fx;
		done = true;
		return;
	}

	xm=0.5*(a+b);
	tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
	
	// Test for done here.
	// printf("T1: %.12g\t%.12g\n",x-xm,b-a);
	if (fabs(x-xm) <= (tol2-0.5*(b-a))) {	
		xmin=x;
		fxmin = fx;
		done = true;
		return;
	}
	if (fabs(e) > tol1) {				
		// Construct a trial parabolic t
		r=(x-w)*(fx-fv);
		q=(x-v)*(fx-fw);
		p=(x-v)*q-(x-w)*r;
		q=2.0*(q-r);
		if (q > 0.0) p = -p;
		q=fabs(q);
		etemp=e;
		e=d;
		if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
			// printf("G");
		}			

		// The above conditions determine the acceptability of the parabolic t. Here we
		// take the golden section step into the larger of the two segments.
		else {
			// Take the parabolic step.
			d=p/q;							
			u=x+d;
			if (u-a < tol2 || b-u < tol2)
				d=SIGN(tol1,xm-x);
			// printf("P");
		}
	} 
	else {
		// printf("g");
		d=CGOLD*(e=(x >= xm ? a-x : b-x));
	}
	u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
	fu=func(u);
	// This is the one function evaluation per iteration.
	if (fu <= fx) {						
		// printf("+\n");
		// Now decide what to do with our function evaluation. 
		if (u >= x) a=x; else b=x;
		SHFT(v,w,x,u);						
		SHFT(fv,fw,fx,fu);
	} 
	else {
		// printf("-\n");
		//Housekeeping follows:
		if (u < x) a=u; else b=u;
		if (fu <= fw || w == x) {
			v=w;
			w=u;
			fv=fw;
			fw=fu;
		} else if (fu <= fv || v == x || v == w) {
			v=u;
			fv=fu;
		}
	} //Done with housekeeping. Back for another iteration. 	

	// printf("Brent: %16.12g  %16.12g  %16.12g  %16.12g\n",a,x,b,fx);
}

bool BrentRoutine::isFinished() {
	return done;
}

BrentRoutine::BrentRoutine(Function *problem,double ax, double bx, double cx,double tol)
: LineMinRoutine(problem) 
{
	this->ax = ax;
	this->bx = bx;
	this->cx = cx;
	this->tol = tol;
	e=0.0;

	// Main program loop.
	iter=0;
	done = false;

	// Initializations...
	a=(ax < cx ? ax : cx);					
	b=(ax > cx ? ax : cx);
	x=w=v=bx;									
	fw=fv=fx=func(x);	
}

const double GoldenRoutine::R = 0.61803399;
const double GoldenRoutine::C = 1.0 - GoldenRoutine::R;

inline double GoldenRoutine::getFAtMinimum() {
	return (f1 < f2) ? f1 : f2;
}

inline double GoldenRoutine::getMinimum() {
	return (f1 < f2) ? x1 : x2;
}

void GoldenRoutine::performIteration() {

	if(fabs(x3-x0) <= tol*(fabs(x1)+fabs(x2))) {
		// End condition
		done = true;
		return;
	}

	if(f2 < f1) {
		SHFT3(x0,x1,x2,R*x1+C*x3)
		SHFT2(f1,f2,func(x2))
	}
	else {
		SHFT3(x3,x2,x1,R*x2+C*x0)
		SHFT2(f2,f1,func(x1))
	}

	// printf("Gold: %16.12g  %16.12g  %16.12g  %16.12g\n",x1,f1,x2,f2);
}

bool GoldenRoutine::isFinished() {
	return done;
}

GoldenRoutine::GoldenRoutine(Function *problem,double _ax, double _bx, double _cx,double tol) 
: LineMinRoutine(problem) 
{
	ax = _ax;
	bx = _bx;
	cx = _cx;
	this->tol = tol;

	// Main program loop.
	done = false;

	// Initializations...
	x0 = ax;
	x3 = cx;
	if(fabs(cx-bx) > fabs(bx-ax)) {
		x1 = bx;
		x2 = bx+C*(cx-bx);
	}
	else {
		x2 = bx;
		x1 = bx-C*(bx-ax);
	}
	f1 = func(x1);
	f2 = func(x2);
}






double Directional1DFunction::evaluate(const Vector &V) {
	double u = V(0);
	X1 = X + N * u;
	return fnd->evaluate(X1);
}

Directional1DFunction::Directional1DFunction(Function *fnd,const Vector &X,const Vector &N) {
	this->X = X;
	this->N = N;
	this->fnd = fnd;
}

double BrentLinearMethod::getFAtMinimum() {
	return linmin ? linmin->getFAtMinimum() : mnbrak->getFAtMinimum();
}

Vector BrentLinearMethod::getMinDirection() {
	return f1d.N * (linmin ? linmin->getMinimum() : mnbrak->getMinimum());
}

Vector BrentLinearMethod::getMinimum() {
	return f1d.X + getMinDirection();
}

BrentLinearMethod::~BrentLinearMethod() {
	if(linmin)
		delete linmin;
	delete mnbrak;
}

void BrentLinearMethod::performIteration() {
	if(!mnbrak->isFinished()) {
		// printf("M ");
		mnbrak->performIteration();
	}
	else if(linmin == NULL) {
		// printf("Brent\n");
		// linmin = new BrentRoutine(&f1d,mnbrak->ax,mnbrak->bx,mnbrak->cx,TOL);
		linmin = new GoldenRoutine(&f1d,mnbrak->ax,mnbrak->bx,mnbrak->cx,TOL);
		//printf("%16.12g\t%16.12g\n",mnbrak->ax,mnbrak->fa);
		//printf("%16.12g\t%16.12g\n",mnbrak->cx,mnbrak->fc);
		//printf("%16.12g\t%16.12g\n",brent->x,brent->fx);
	}
	else if(!linmin->isFinished()) {
		//printf("b ");
		linmin->performIteration();
		//printf("%16.12g\t%16.12g\n",brent->x,brent->fx);
	}
	else {
		done = true;
	}
}

BrentLinearMethod::BrentLinearMethod(Function *problem,const Vector &x,const Vector &n)
: f1d(problem,x,n)
{
	// value = 0;	
	mnbrak = new MinBrakRoutine(&f1d,0.0,1.0);
	linmin = NULL;
	done = false;
/*
	double xx = 1.0,xmin,fx,fb,fa,bx,ax = 0.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb);
	value = brent(ax,xx,bx,TOL,&xmin);
	n *= xmin;
	x += n;*/
}


/*******************************************************************
 Solve a problem of finiding s that optimizes f(P+sn) 
 ******************************************************************
BrentLinearMethod::BrentLinearMethod(Function &problem,Vector inX,Vector inN) :  
p(problem) 
{
	x = inX;
	xt = x;
	n = inN;
	value = 0;
}

void BrentLinearMethod::run() {
	double xx = 1.0,xmin,fx,fb,fa,bx,ax = 0.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb);
	value = brent(ax,xx,bx,TOL,&xmin);
	n *= xmin;
	x += n;
}

// Given a function func, and given distinct initial points ax and bx, this routine searches in
// the downhill direction (defned by the function as evaluated at the initial points) and returns
// new points ax, bx, cx that bracket a minimum of the function. Also returned are the function
// values at the three points, fa, fb, and fc.
// Here GOLD is the default ratio by which successive intervals are magnifed; 
// GLIMIT is the maximum magnifcation allowed for a parabolic-ft step.
void BrentLinearMethod::mnbrak(double *ax, double *bx, double *cx, 
									  double *fa, double *fb, double *fc)
{
	double ulim,u,r,q,fu,dum;
	*fa=func(*ax);
	*fb=func(*bx);
	if (*fb > *fa) { 
		//Switch roles of a and b so that we can go
		// downhill in the direction from a to b. 
		SHFT(dum,*ax,*bx,dum);
		SHFT(dum,*fb,*fa,dum);
	}
	
	// First guess for c.
	*cx=(*bx)+GOLD*(*bx-*ax); 
	*fc=func(*cx);
	
	while (*fb > *fc) { 
		// Keep returning here until we bracket.
		// Compute u by parabolic extrapolation from a; b; c. TINY is used to 
		// prevent any possible division by zero.
		r=(*bx-*ax)*(*fb-*fc); 
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		
		// We won't go farther than this. Test various possibilities:
		if ((*bx-u)*(u-*cx) > 0.0) { 
			// Parabolic u is between b and c: try it.
			fu=func(u);
			if (fu < *fc) { 
				// Got a minimum between b and c.
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} 
			else if (fu > *fb) { 
				// Got a minimum between between a and u.
				*cx=u;
				*fc=fu;
				return;
			}
			
			u=(*cx)+GOLD*(*cx-*bx); 
			// Parabolic t was no use. Use default mag-nication. 
			fu=func(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) { 
			//Parabolic t is between c and its allowed limit. 
			fu=func(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx));
				SHFT(*fb,*fc,fu,func(u));
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) { 
			//Limit parabolic u to maximumallowed value. 
			u=ulim;
			fu=func(u);
		} else { 
			// Reject parabolic u, use default magnica-tion. 
			u=(*cx)+GOLD*(*cx-*bx);
			fu=func(u);
		}
		SHFT(*ax,*bx,*cx,u); // Eliminate oldest point and continue.
		SHFT(*fa,*fb,*fc,fu);
	}
}

// Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
// between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
// the minimum to a fractional precision of about tol using Brent's method. The abscissa of
// the minimum is returned as xmin, and the minimum function value is returned as brent, the
// returned function value.
double BrentLinearMethod::brent(double ax, double bx, double cx,
										double tol, double *xmin)
{
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;								
	// This will be the distance moved on the step before last.
	// a and b must be in ascending order,
	// but input abscissas need not be. 
	
	// Initializations...
	a=(ax < cx ? ax : cx);					
	b=(ax > cx ? ax : cx);
	x=w=v=bx;									
	fw=fv=fx=func(x);
	
	// Main program loop.
	for (iter=1;iter<=ITMAX;iter++) {	
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		
		// Test for done here.
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {	
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {				
			// Construct a trial parabolic t
			r=(x-w)*(fx-fv);
		q=(x-v)*(fx-fw);
		p=(x-v)*q-(x-w)*r;
		q=2.0*(q-r);
		if (q > 0.0) p = -p;
		q=fabs(q);
		etemp=e;
		e=d;
		if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		
		// The above conditions determine the acceptability of the parabolic t. Here we
		// take the golden section step into the larger of the two segments.
		else {
			// Take the parabolic step.
			d=p/q;							
			u=x+d;
			if (u-a < tol2 || b-u < tol2)
				d=SIGN(tol1,xm-x);
		}
		} 
		else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=func(u);
		// This is the one function evaluation per iteration.
		if (fu <= fx) {						
			// Now decide what to do with our function evaluation. 
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u);						
			SHFT(fv,fw,fx,fu);
		} 
		else {
			//Housekeeping follows:
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		} //Done with housekeeping. Back for another iteration. 
	}

	dassert(0);
	*xmin=x; //Never get here.
	return fx;
}
*/
/*******************************************************************
 Solve a problem of finiding s that optimizes f(P+sn) with f'
 ******************************************************************/
BrentDerivativeMethod::BrentDerivativeMethod(DifferentiableFunction &problem,Vector inX,Vector inN) :  
p(problem) 
{
	x = inX;
	n = inN;
	xt = x;
	value = 0;
}


void BrentDerivativeMethod::run() {
	double xx = 1.0,xmin,fx,fb,fa,bx,ax = 0.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb);
	value = dbrent(ax,xx,bx,TOL,&xmin);
	n *= xmin;
	x += n;
}

// Given a function func, and given distinct initial points ax and bx, this routine searches in
// the downhill direction (defned by the function as evaluated at the initial points) and returns
// new points ax, bx, cx that bracket a minimum of the function. Also returned are the function
// values at the three points, fa, fb, and fc.
// Here GOLD is the default ratio by which successive intervals are magnifed; 
// GLIMIT is the maximum magnifcation allowed for a parabolic-ft step.
void BrentDerivativeMethod::mnbrak(double *ax, double *bx, double *cx, 
									  double *fa, double *fb, double *fc)
{
	double ulim,u,r,q,fu,dum;
	*fa=func(*ax);
	*fb=func(*bx);
	if (*fb > *fa) { 
		//Switch roles of a and b so that we can go
		// downhill in the direction from a to b. 
		SHFT(dum,*ax,*bx,dum);
		SHFT(dum,*fb,*fa,dum);
	}
	
	// First guess for c.
	*cx=(*bx)+GOLD*(*bx-*ax); 
	*fc=func(*cx);
	
	while (*fb > *fc) { 
		// Keep returning here until we bracket.
		// Compute u by parabolic extrapolation from a; b; c. TINY is used to 
		// prevent any possible division by zero.
		r=(*bx-*ax)*(*fb-*fc); 
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		
		// We won't go farther than this. Test various possibilities:
		if ((*bx-u)*(u-*cx) > 0.0) { 
			// Parabolic u is between b and c: try it.
			fu=func(u);
			if (fu < *fc) { 
				// Got a minimum between b and c.
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} 
			else if (fu > *fb) { 
				// Got a minimum between between a and u.
				*cx=u;
				*fc=fu;
				return;
			}
			
			u=(*cx)+GOLD*(*cx-*bx); 
			// Parabolic t was no use. Use default mag-nication. 
			fu=func(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) { 
			//Parabolic t is between c and its allowed limit. 
			fu=func(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx));
				SHFT(*fb,*fc,fu,func(u));
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) { 
			//Limit parabolic u to maximumallowed value. 
			u=ulim;
			fu=func(u);
		} else { 
			// Reject parabolic u, use default magnica-tion. 
			u=(*cx)+GOLD*(*cx-*bx);
			fu=func(u);
		}
		SHFT(*ax,*bx,*cx,u); // Eliminate oldest point and continue.
		SHFT(*fa,*fb,*fc,fu);
	}
}


// Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
// between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
// the minimum to a fractional precision of about tol using Brent's method. The abscissa of
// the minimum is returned as xmin, and the minimum function value is returned as brent, the
// returned function value.
double BrentDerivativeMethod::dbrent(double ax, double bx, double cx,
												 double tol, double *xmin)
{
	// Will be used as ags for whether proposed steps are acceptable or not.
	int iter,ok1,ok2;  
	double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
	double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
	
	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=dfunc(x,dx);
	dw=dv=dx; 
	// All our housekeeping chores are doubled by the necessity of moving
	// derivative values around as well as function values.
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol1=tol*fabs(x)+ZEPS;
		tol2=2.0*tol1;
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			// Initialize these d's to an out-of-bracket value. 
			d1=2.0*(b-a); 
			d2=d1;
			// Secant method with one point.  Which of these two estimates of d shall we take? 
			// We will insist that they be within the bracket, and on the side pointed to 
			// by the derivative at x:
			if (dw != dx) d1=(w-x)*dx/(dx-dw);
			
			if (dv != dx) d2=(v-x)*dx/(dx-dv); 
			
			u1=x+d1;
			u2=x+d2;
			ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
			ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
			// Movement on the step before last.
			olde=e; 
			e=d;

			if (ok1 || ok2) { 
				// Take only an acceptable d, and if both are acceptable, then take
				// the smallest one.
				if (ok1 && ok2)
					d=(fabs(d1) < fabs(d2) ? d1 : d2);
				else if (ok1)
					d=d1;
				else
					d=d2;
				if (fabs(d) <= fabs(0.5*olde)) {
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				} 
				else { 
					// Bisect, not golden section.
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
					// Decide which segment by the sign of the derivative.
				}
			} 
			else {
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
			}
		} 
		else {
			d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}

		// Optimization (paul)
		double _dfu,_fu;
		_fu = dfunc(u,_dfu);

		if (fabs(d) >= tol1) {
			u=x+d;
			fu=_fu;
		} 
		else {
			u=x+SIGN(tol1,d);
			fu=_fu;
			if (fu > fx) { 
				// If the minimum step in the downhill direction takes us uphill, then
				// we are done.
				*xmin=x;
				return fx;
			}
		}

		// Now all the housekeeping, sigh.
		du=_dfu; 
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			MOV3(v,fv,dv, w,fw,dw)
				MOV3(w,fw,dw, x,fx,dx)
				MOV3(x,fx,dx, u,fu,du)
		} 
		else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				MOV3(v,fv,dv, w,fw,dw)
					MOV3(w,fw,dw, u,fu,du)
			} 
			else if (fu < fv || v == x || v == w) {
				MOV3(v,fv,dv, u,fu,du)
			}
		}
	}

	// printf("Too many iterations in routine dbrent");
	dassert(0);
	return 0.0; // Never get here.
}

const double BrentRootFinder::ROOT_ERROR = -4.56e123;

BrentRootFinder::BrentRootFinder(Function1D &fn) : f(fn) {
}

double BrentRootFinder::zbrent(double x1,double x2,double tol) {
	int iter;
	double a = x1,b = x2,c = x2,d,e,min1,min2;
	double fa = f.evaluate(a),fb = f.evaluate(b),fc,p,q,r,s,tol1,xm;

	if(fabs(fa) < tol) return a;		
	if(fabs(fb) < tol) return b;

	if((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
		return ROOT_ERROR;
	}

	fc = fb;
	for(iter = 1;iter <= ITMAX;iter++) {
		if((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c = a;
			fc = fa;
			e=d=b-a;
		}
		if(fabs(fc) < fabs(fb)) {
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		tol1 = 2.0 * ZEPS * fabs(b)+0.5*tol;
		xm = 0.5 * (c-b);

		if(fabs(xm) < tol1 || fb == 0.0)
			return b;

		if(fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s = fb / fa;
			if(a == c) {
				p = 2.0 * xm * s;
				q = 1.0 - s;
			}
			else {
				q = fa / fc;
				r = fb / fc;
				p = s * (2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q = (q - 1.0)*(r-1.0)*(s - 1.0);
			}
			if(p > 0.0)
				q = -q;
			p = fabs(p);
			min1 = 3.0 * xm *q - fabs(tol1 * q);
			min2 = fabs(e * q);
			if(2.0 * p < (min1 < min2 ? min1 : min2)) {
				e = d;
				d = p/q;
			}
			else {
				d = xm;
				e = d;
			}
		} 
		else {
			d = xm;
			e = d;
		}
		a = b;
		fa = fb;
		if(fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		fb = f.evaluate(b);
	}

	// Too many iterations
	return ROOT_ERROR;
}

// End namespace
NAMESPACE_PAULY_END
