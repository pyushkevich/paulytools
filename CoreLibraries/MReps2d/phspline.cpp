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
#include <math.h>
#include "phspline.h"
#include <BrentLinearMethod.h>
#include <EvolutionaryStrategy.h>
#include <ConjugateGradientMethod.h>
#include <SimplexMethod.h>

#include <algorithm>
#include <iostream>

// Begin namespace
NAMESPACE_PAULY_START

/***********************************************************************
 * Abstract representation with polygons
 ***********************************************************************/
// Is a circle fully inside the object?  (A rough test)
bool ContinuousBoundary::containsCircle(const Vector2D &x,double r,int samples) {
   double tStep = tmax() / samples;
   double t = 0;
   double r2 = r*r;

   for(int i=0;i<samples;i++) {
      Vector2D xt;
      t += tStep;
      getInterpolation(t,xt);

      if(x.distanceToSqr(xt) < r2) {
         return false;
      }
   }

   return true;
}

// A comment to test CVS
double ContinuousBoundary::mediality(double s,double t) {
	Vector2D x1,x2,t1,t2,n1,n2;

   getInterpolation(s,x1,n1);
	n1.normalize();
	t1 = n1.getNormal();

	getInterpolation(t,x2,n2);
	n2.normalize();
	t2 = n2.getNormal();
	
	Vector2D dx = x2 - x1, dn = n1 - n2;
	//dx.normalize();
	//dn.normalize();

	double angle = acos(fabs(dx.dotProduct(dn) / (dx.twoNorm() * dn.twoNorm())));
	angle *= 180 /M_PI;
	return angle;
}

Vector ContinuousBoundary::medialityJet(double s,double t) {
	Vector jet(6);
	Vector2D xs,xt,x1s,x1t,x2s,x2t;

   getInterpolationJet(s,xs,x1s,x2s);
	getInterpolationJet(t,xt,x1t,x2t);

	// Lengths and curvatures
	double ls = x1s.twoNorm(), lt = x1t.twoNorm();
	double ks = (x2s.x*x1s.y - x2s.y*x1s.x) / (ls*ls*ls);
	double kt = (x2t.x*x1t.y - x2t.y*x1t.x) / (lt*lt*lt);
	
	// Tangents
	Vector2D Ts = x1s / ls, Tt = x1t / lt;
	Vector2D T1s = Vector2D(x1s.y * ks,-x1s.x * ks);
	Vector2D T1t = Vector2D(x1t.y * kt,-x1t.x * kt);
	
	// Function
	jet(0) = (xs-xt).dotProduct(Ts-Tt);

	// Gradient
	jet(1) = x1s.dotProduct(Ts-Tt) + T1s.dotProduct(xs-xt);
	jet(2) = -x1t.dotProduct(Ts-Tt) - T1t.dotProduct(xs-xt);

	// Second derivative


	
	return jet;
} 

void ContinuousBoundary::getMedialPoint(double s,double t,Vector2D &v,double &r,double &axialAngle,double &objectAngle) {
	Vector2D X1,X2,t1,t2,n1,n2;

   getInterpolation(s,X1,n1);
	n1.normalize();
	t1 = n1.getNormal();

	getInterpolation(t,X2,n2);
	n2.normalize();
	t2 = n2.getNormal();

	double x1 = X1.x;
	double y1 = X1.y;
	double x2 = X2.x;
	double y2 = X2.y;

	double u1 = n1.x;
	double v1 = n1.y;
	double u2 = n2.x;
	double v2 = n2.y;

	double r1,r2;
	double den = v1 * u2 - u1 * v2;
	if(fabs(den) < 0.0001) {
		r1 = r2 = X1.distanceTo(X2) / 2;
	}
	else {
		r1 = (u2 * (y1-y2) - v2 * (x1-x2)) / den;
		r2 = (u1 * (y1-y2) - v1 * (x1-x2)) / den;
	}

	v = (X1 - n1 * r1 + X2 - n2 * r2) * 0.5;
	r = (r1 + r2 ) * 0.5;

   double slope1 = (X1 - v).getSlope();
   double slope2 = (X2 - v).getSlope();
   
	axialAngle = -(slope1 + slope2) / 2.0;
   objectAngle = -(slope1 - slope2) / 2.0;
}


// A function that is built on top of the spline but returns absolute value of
// the objective function
class CoreTrackPbm : public Function {
   ContinuousBoundary &boundary;
public:
   CoreTrackPbm(ContinuousBoundary &b) : boundary(b) {
   }

   double evaluate(const Vector &x) {
      return fabs(boundary.mediality(x(0),x(1)));
   }
};

class CoreTrackCirclePbm : public Function {
   ContinuousBoundary &boundary;
   Vector center;
   Vector nose;
   Vector normal;
public:
   CoreTrackCirclePbm(ContinuousBoundary &b,Vector center,Vector nose) : boundary(b) {
      this->nose = nose;
      this->center = center;
      normal = Vector(2,-nose(1),nose(0)); 
   }

   Vector getPoint(double x) {
      double theta = 1.5*atan(x);
      return center + cos(theta)*nose - sin(theta)*normal;
   }
   
   double evaluate(const Vector &x) {
      double theta = 2*atan(x(0));
      Vector p =getPoint(x(0));
      return fabs(boundary.mediality(p(0),p(1)));
   }
};


// Trace core from starting position (s,t) storing intermediate points into the
// parameter vector.
int ContinuousBoundary::traceCore(double s,double t,vector<Vector> &out) {
   const double epsilon = 1.0e-4;
   const double stepSize = 1.0e-3;
   const double brentStep = 1.0e-4;
   const double initBrentStep = 1.0e-3;

   // Make a numerical function  for the problem
   CoreTrackPbm p(*this);
   NumericalFunction nf(p,epsilon);

   // First, find an initial point on the zero-set of the function.  We assume that
   // the point given to us is fairly close to us
   Vector x(2,s,t);
   Vector grad(2);
   nf.computeOneJet(x,grad);
   grad.normalize();

   // We look for starting point along the gradient vector
   BrentLinearMethod blm(&nf,x,initBrentStep*grad);
   while(!blm.isFinished())
       blm.performIteration();
   x = blm.getMinimum();

   // Look at the value of the function there.  If it is not zero, we've messed up and could
   // not find a good starting point
   if(blm.getFAtMinimum() > epsilon) {
      cerr << "Could not find starting point along gradient from initial point\n";
      return -1;
   }

   // Initialize our vector with the starting point
   out.clear();
   out.push_back(x);
	Vector xStart = x;

   // We use two vectors through the loop.  x is the last point we located on the track.
   // t is the vector between last x and x.
   
   // A zero-crossing has been found.  Ve now have a start for our 'core'.  Track in each direction.
   for(int dir=-1;dir<=1;dir+=2) {
      Vector t(2,dir*grad(1),-dir*grad(0)); 
	
		// Last position at which we stopped
		Vector xLast = xStart;

      // Run along the core in each direction
      while(1) {
         bool found = false;

         // Find minimum along the gradient at the new position
         for(int i=0;i<10;i++) {
            double ss = (i==0) ? stepSize : fabs(getGaussianRnd(0,stepSize));
            CoreTrackCirclePbm ctc(*this,x+ss*t,ss*t);   
            BrentLinearMethod blm(&ctc,Vector(1,0.0),Vector(1,1.0));            
			while(!blm.isFinished())
				blm.performIteration();

            // If the value is non zero, we've got a problem
            if(blm.getFAtMinimum() < epsilon) {
			   t = ctc.getPoint(blm.getMinimum()(0)) - x;
               t.normalize();
               x = ctc.getPoint(blm.getMinimum()(0));

               found = true;               
               cout << ".";
               break;
            }

            if(i==0) {
               cout << "\n (" << x(0) << "," << x(1) << ")\n";
            }
            cout << "?";            
         }

/*
         // Find minimum along the gradient at the new position
         for(int i=0;i<4;i++) {
            double ss = (i==0) ? stepSize : fabs(getGaussianRnd(0,stepSize));
            BrentLinearMethod blm(nf,x+ss*t,ss*Vector(2,-t(1),t(0)));
            blm.run();

            // If the value is non zero, we've got a problem
            if(blm.value < epsilon) {
               t = blm.x - x;
               t.normalize();
               x = blm.x;

               found = true;               
               cout << ".";
               break;
            }

            if(i==0) {
               cout << "\n (" << x(0) << "," << x(1) << ")\n";
            }
            cout << "?";
            
         }

         if(!found) {
            // Now try something different... Look around the circle
            CoreTrackCirclePbm ctc(*this,x+stepSize*t,stepSize*t);   
            for(double k=-1.5;k<=1.5;k+=0.5) {
               if(k==0.0) continue;
               BrentLinearMethod blm(ctc,Vector(1,0.0),Vector(1,k));
               blm.run();
               
               // If the value is non zero, we've got a problem
               if(blm.value < epsilon) {
                  t = ctc.getPoint(blm.x(0)) - x;
                  t.normalize();
                  x = ctc.getPoint(blm.x(0));
                  
                  found = true;
                  cout << "!";
                  break;
               }
               
               cout << "@";
            }
         }
*/
         if(!found) {            
            cerr << "\nLost track of the core\n";
            break;
         }

         // Append x to the vector
         out.push_back(x);

			// Constrain s,t parameters to unit box
			bool flip[2] = {false,false};
			for(int p=0;p<2;p++) {
				while(x(p) >= tmax()) {
					x(p) -= tmax();
					flip[p] = true;
				}
				while(x(p) < 0) {
					x(p) += tmax();
					flip[p] = true;
				}
			}
			
			xLast.t().print();
         x.t().print();

			// Termination condition:
			// If last x was above the s=t line and current x is at the s=t line, terminate
			if((xLast(0) - xLast(1)) * (x(0)-x(1)) <= 0) {
				if(flip[0] || flip[1]) {
					cerr << "\nCrossed tmax() boundary!\n";
				}
				else {
					cerr << "\nEnd has been reached\n";
					break;
				}
			}

         // If s and t are less than an epsilon away, we say that an end has been reached
         //if(fabs(fmod(x(0),tmax()) - fmod(x(1),tmax())) < epsilon) {
         //   cerr << "\nEnd has been reached\n";
         //   break;
         //}

			xLast = x;
      }

      // Reverse the out vector if dir==-1
      if(dir == -1) {
         x = out[0];
         reverse(out.begin(),out.end());
      }
   }

   return 0;
}

// Create a polygon boundary that corrensponds to the spline
void ContinuousBoundary::constructPolygonBoundary(PolygonBoundary &pb,int numPolygons) {
   pb.vertices.clear();

   double tStep = tmax() / numPolygons;
   double t = 0;

   for(int i=0;i<numPolygons;i++) {
      Vector2D x,n;
      getInterpolation(t,x,n);
      pb.addVertex(x,n);

      t+=tStep;
   }
}


/***********************************************************************
 * Boundary representation with polygons
 ***********************************************************************/
// Load the boundary representation from a text file (2 or 4 doubles per line),
// depending on whether the normals are needed
// (Scale all points by scaleFactor is an available option)
// Returns negative value of failure
int PolygonBoundary::loadPoints(char *fname,bool hasNormals,double scaleFactorX,double scaleFactorY) {
   FILE *f = fopen(fname,"rt");
   if(!f)
      return -1;
   
   vertices.clear();
   
   char buffer[256];
   while(!feof(f)) {
      double v[4];
      Vector2D x,n;
      
      fgets(buffer,255,f);
      char *p = strtok(buffer," \r\n\t");
      int i = -1;
      while(p && i<4) {
         v[++i] = atof(p);
         p = strtok(NULL," \r\n\t");
      }
      
      x.x = v[0];x.y = v[1];
      if(hasNormals) 
         n.x = v[2];n.y = v[3];

      vertices.insert(vertices.end(),PolygonVertex(x.scaledBy(scaleFactorX,scaleFactorY),n.scaledBy(scaleFactorX,scaleFactorY)));
   }

   return 0;
}


// Save the boundary representation from a text file (2 or 4 doubles per line),
int PolygonBoundary::savePoints(char *fname,double scaleFactorX,double scaleFactorY) {
   FILE *f = fopen(fname,"wt");
   if(!f)
      return -1;
   
   for(int i=0;i<vertices.size();i++) {
      fprintf(f,"%lg %lg %lg %lg\n",scaleFactorX * vertices[i].x.x,
                                    scaleFactorY * vertices[i].x.y,
                                    scaleFactorX * vertices[i].n.x,
                                    scaleFactorY * vertices[i].n.y);
   }

   return 0;
}

// Compute normals numerically for points that don't have normals (by approximating
// the tangent direction)
void PolygonBoundary::approximateNormals(bool clockWise) {
   double dAlpha = (clockWise) ? M_PI_2 : -M_PI_2;
   int n = vertices.size();

   for(int i=0;i<n;i++) {
      PolygonVertex v = vertices[i];
      PolygonVertex p = vertices[(i-1) % n];
      PolygonVertex s = vertices[(i+1) % n];

      
      // Find the tangent angle
      Vector2D tng = s.x - p.x;
      
      // Computed normal
      Vector2D nrm = (clockWise) ? tng.getNormal() : -tng.getNormal();
      nrm.normalize();

      // Assign the angle if no normal exisits
      if(v.n.x == 0 && v.n.y == 0) {
         v.n = nrm;
      }
      else {
         // Make sure that normal is in the right direction
         if(nrm.dotProduct(v.n) < 0) {
            v.n =  -v.n;
         }
      }
   }
}

/***********************************************************************
 * A Fourier spline representation
 ***********************************************************************/
FourierBoundary::FourierBoundary() {
   T;
}

void FourierBoundary::setTransform(const Vector2D &dx,const Vector2D &sx) {
   T = Transform2D::frameChange(dx,sx);
}

// File formats for fourier boundaries
const int FourierBoundary::ONEPERLINE = 0;
const int FourierBoundary::MATHEMATICA = 1;
const int FourierBoundary::FOURPERLINE = 2;


// Load the boundary representation from a text file 
// that contains fourier coefficients in a list 
// (first line is the order of highest coefficient, and the consequent lines
// are coefficients of x,y, alternating, of order 1,-1,2,-2,...
// It's the for format used by Sean Ho
int FourierBoundary::load(char *fname,int format) {
   char buffer[256];
	FILE *f = fopen(fname,"rt");
   if(!f)
      return -1;

	pos.clear();
	neg.clear();

   // Insert zeros for zero order coefficients
   if(format == ONEPERLINE) {
		// Skip the first two lines (zero order coefficients are always zero)
		fgets(buffer,256,f);
		fgets(buffer,256,f);
		pos.push_back(cmplx(0.0,0.0));
		neg.push_back(cmplx(0.0,0.0));
   
		// Read coefficients
		while(!feof(f)) {
			double xp,yp,xn,yn;
			fgets(buffer,256,f);xp = atof(buffer);
			fgets(buffer,256,f);yp = atof(buffer);
			fgets(buffer,256,f);xn = atof(buffer);
			fgets(buffer,256,f);yn = atof(buffer);
			
			pos.push_back(cmplx(xp,yp));
			neg.push_back(cmplx(xn,yn));
		}
	}
	else if(format == MATHEMATICA) {
		// Skip the first line
		fgets(buffer,256,f);
		pos.push_back(cmplx(0.0,0.0));
		neg.push_back(cmplx(0.0,0.0));
		
		// Read coefficients
		while(!feof(f)) {
			double xp,yp,xn,yn;
			fgets(buffer,256,f);
			
			if(strlen(buffer) > 4) {
				sscanf(buffer,"{%lg,%lg,%lg,%lg}",&xp,&yp,&xn,&yn);
				neg.push_back(cmplx(xp,xn));
				pos.push_back(cmplx(yp,yn));
			}
		}	

		// Reverse the order
		//reverse(pos.begin()+1,pos.end());
		//reverse(neg.begin()+1,neg.end());
	}
	else if(format == FOURPERLINE) {
		// Skip the first line
		fgets(buffer,256,f);
		fgets(buffer,256,f);
		fgets(buffer,256,f);
		
		pos.push_back(0);
		neg.push_back(0);

		// Read coefficients
		while(!feof(f)) {
			fgets(buffer,256,f);
			if(strlen(buffer) > 20) {				
				double w = atof(strtok(buffer," \t\n"));
				double x = atof(strtok(NULL," \t\n"));
				double y = atof(strtok(NULL," \t\n"));
				double z = atof(strtok(NULL," \t\n"));

				pos.push_back(cmplx(w+z,-y+x)/2.0);
				neg.push_back(cmplx(w-z,-y-x)/2.0);
			}
		}	

		// Reverse the order
		//reverse(pos.begin()+1,pos.end());
		//reverse(neg.begin()+1,neg.end());
	}


   fclose(f);
   return 0;
};

// Load the boundary representation from a text file 
// that contains fourier coefficients in a list 
// (first line is the order of highest coefficient, and the consequent lines
// are coefficients of x,y, alternating, of order 1,-1,2,-2,...
// It's the for format used by Sean Ho
int FourierBoundary::save(char *fname,int format) {
   FILE *f = fopen(fname,"wt");
   if(!f)
      return -1;

   if(format == ONEPERLINE) {
		// Skip the first two lines (zero order coefficients are always zero)
		fprintf(f,"0\n0\n");
      
		// Read coefficients
		for(int i=1;i<pos.size();i++) {
			fprintf(f,"%lg\n%lg\n",real(pos[i]),imag(pos[i]));
			fprintf(f,"%lg\n%lg\n",real(neg[i]),imag(neg[i]));
		}
	}
	else if(format == MATHEMATICA) {
		// Skip the first two lines (zero order coefficients are always zero)
		fprintf(f,"{\n");
      
		// Read coefficients
		for(int i=1;i<pos.size();i++) {
			fprintf(f,"{%lg,%lg,%lg,%lg}",real(pos[i]),imag(pos[i]),real(neg[i]),imag(neg[i]));
			if(i<pos.size()-1) 
				fprintf(f,",");
			fprintf(f,"\n");		
		}
	}

   fclose(f);
   return 0;
};

void FourierBoundary::getInterpolation(double t,Vector2D &xOut,Vector2D &nOut) {
   // Order is known now
   int order = pos.size()-1;

   // We need i
   static cmplx i(0,1);

   cmplx xt,xtPrime; 
   for(int k=0;k<=order;k++) {
      cmplx M = i*(2.0*M_PI*k);
      xt += pos[k]*exp(M*t) + neg[k]*exp(-M*t);
      xtPrime += M*pos[k]*exp(M*t) - M*neg[k]*exp(-M*t);
   }

   // Fill vectors
   Vector2D x(xt.real(),xt.imag(),1);
   Vector2D n(xtPrime.real(),xtPrime.imag(),0);

   // Multiply by transform
   xOut = T * x;
   nOut = (T * n).getNormal();
	nOut.normalize();
}

void FourierBoundary::getInterpolation(double t,Vector2D &xOut) {
   // Order is known now
   int order = pos.size()-1;

   // We need i
   static cmplx i(0,1);

   cmplx xt,xtPrime; 
   for(int k=1;k<=order;k++) {
      cmplx M = i*(2.0*M_PI*k);
      xt += pos[k]*exp(M*t) + neg[k]*exp(-M*t);
   }

   // Fill vectors
   Vector2D x(xt.real(),xt.imag(),1);

   // Multiply by transform
   xOut = T * x;
}

void FourierBoundary::getInterpolationJet(double t,Vector2D &x,Vector2D &x1,Vector2D &x2) {
   // Order is known now
   int order = pos.size()-1;

   // We need i
   static cmplx i(0,1);
	
	cmplx z,z1,z2;

   for(int k=0;k<=order;k++) {
      cmplx M = i*(2.0*M_PI*k);
		cmplx E1 = pos[k]*exp(M*t), E2 = neg[k]*exp(-M*t);

      z += E1 + E2;
      z1 += M*(E1 - E2);
		z2 += M*M*(E1 + E2);
   }

   // Fill vectors
   x  = T * Vector2D(z.real(), z.imag() ,1);
	x1 = T * Vector2D(z1.real(),z1.imag(),0);
	x2 = T * Vector2D(z2.real(),z2.imag(),0);
}

// Construct a spline from a list of point coordinates and normals 
void FourierBoundary::buildFromPB(PolygonBoundary &pb,int numSamples) {
   // Number of coefficients computed
   int N = numSamples;
   pos.resize(N+1);
   neg.resize(N+1);

   // We need i
   static cmplx i(0,1);

   int M = pb.vertices.size();
   double L = 0;
   
   vector<cmplx> ukp,ukn;
   for(int k=0;k<M;k++) {
      L += pb.vertices[k].x.distanceTo(pb.vertices[(k+1)%M].x);
   }

   for(int n=1;n<=N;n++) {
      pos[n] = 0;
      neg[n] = 0;
      double s = 0.0;

      for(int k=0;k<M;k++) {
         double t = 1.0 * k / M;
         cmplx zk1(pb.vertices[(k+1)%M].x.x,pb.vertices[(k+1)%M].x.y);
         cmplx zk(pb.vertices[k].x.x,pb.vertices[k].x.y);
         cmplx dzk = zk1-zk;
         cmplx zkPrime = (dzk * L) / (2*M_PI*abs(dzk));
         
         cmplx uk_p = exp(-i*2.0*M_PI*s*(1.0*n)/L);
         cmplx uk_n = exp(i*2.0*M_PI*s*(1.0*n)/L);
         
         s += abs(dzk);
         
         cmplx uk1_p = exp(-i*2.0*M_PI*s*(1.0*n)/L);                  
         cmplx uk1_n = exp(i*2.0*M_PI*s*(1.0*n)/L);                  

         pos[n] += (1.0/(n*n*M_PI))*zkPrime*(uk1_p-uk_p);
         neg[n] += (1.0/(n*n*M_PI))*zkPrime*(uk1_n-uk_n);
      }
   }   
}


// This method returns a vector containing the coefficients of the fourier representation
Vector FourierBoundary::getFourierCoeff() const {
	int n = neg.size()-1;
	Vector v(n*4);
	int idx = 0;
	for(int i=1;i<=n;i++) {
		v(idx++) = neg[i].real();
		v(idx++) = neg[i].imag();
		v(idx++) = pos[i].real();
		v(idx++) = pos[i].imag();
	}
	return v;
}

// This loads the fourier boundary from a vector containing the representation
void FourierBoundary::setFourierCoeff(const Vector &coeff) {
	// Number of coefficients computed
   int n = coeff.size()/4;
   pos.resize(n+1);
   neg.resize(n+1);

	// Push the zeros
	pos[0] = 0;neg[0] = 0;

	// Load the coefficients
	int idx = 0;
	for(int i=1;i<=n;i++) {
		neg[i] = cmplx(coeff(idx++),coeff(idx++));
		pos[i] = cmplx(coeff(idx++),coeff(idx++));
	}
}

/***********************************************************************
 * A Continous Closed Spline
 ***********************************************************************/
// A blank constructor
PHSpline::PHSpline() {
};

void PHSpline::buildFromPB(PolygonBoundary &pb,int numSamples) {
   // Run through all the polygons.  This is a dumb method that just ignores
   // polygons between the polygons hat correspond to actual boundary points
   double polsPerSpline = ((double)pb.vertices.size()) / numSamples;

   // Make sure all normals are available all
   // pb.approximateNormals();

   // Remove segments
   segments.clear();
   
   // Add splines to the vector
   for(int i=0;i<numSamples;i++) {
      // Get vertex, next vertex
      int nextIndex = (i==numSamples-1) ? 0 : (int)((i+1)*polsPerSpline);

      PolygonVertex &pv = pb.vertices[(int) (i*polsPerSpline)];
      PolygonVertex &npv = pb.vertices[nextIndex];

      // Construct a quintic
      //segments.insert(segments.end(),PHQuintic(-0.5,-0.5,0.5,0.5,M_PI/2,M_PI/2));
       
      // Scale normals to reasonable size
      double dx = (npv.x-pv.x).twoNorm();
      pv.n.normalize();
      npv.n.normalize();
      pv.n *= dx / 4;
      npv.n *= dx / 4;

      segments.insert(segments.end(),CBoundlet());
      segments[segments.size()-1].compute(pv.x,npv.x,pv.n,npv.n);
   }
}

void PHSpline::buildFromPB(PolygonBoundary &pb) {
   // Remove segments
   segments.clear();
   
   // Add splines to the vector
   for(int i=0;i<pb.vertices.size();i++) {
      // Get vertex, next vertex
      int nextIndex = (i==pb.vertices.size()-1) ? 0 : i+1;
      PolygonVertex &pv = pb.vertices[i];
      PolygonVertex &npv = pb.vertices[nextIndex];

      // Construct a quintic
      //segments.insert(segments.end(),PHQuintic(-0.5,-0.5,0.5,0.5,M_PI/2,M_PI/2));
       
      // Scale normals to reasonable size
      double dx = (npv.x-pv.x).twoNorm();
      pv.n.normalize();
      npv.n.normalize();
      pv.n *= dx / 4;
      npv.n *= dx / 4;

      segments.insert(segments.end(),CBoundlet());
      segments[segments.size()-1].compute(pv.x,npv.x,pv.n,npv.n);
   }
}

// Get the position along the spline at given t (t ranges from 0 to tmax()) with same
// point corresponding to t=0 and t=tmax()
void PHSpline::getInterpolation(double t,Vector2D &x) {
   Vector2D n;
   
   t = fmod(t,segments.size());
   t = t < 0 ? t + segments.size() : t;

   segments[(int)floor(t)].getInterpolation(t-floor(t),x,n);
}

// Get the position along the spline at given t (t ranges from 0 to tmax()) with same
// point corresponding to t=0 and t=tmax()
void PHSpline::getInterpolation(double t,Vector2D &x,Vector2D &n) {  
   t = fmod(t,segments.size());
   t = t < 0 ? t + segments.size() : t;

   segments[(int)floor(t)].getInterpolation(t-floor(t),x,n);
   n = n.getNormal();
}

void PHSpline::getInterpolationJet(double t,Vector2D &x,Vector2D &x1,Vector2D &x2) {
	
}


int PHSpline::load(char *fname) {
   FILE *f = fopen(fname,"rt");
   if(!f)
      return -1;
   
   char buffer[256];
   
   segments.clear();
   
   fgets(buffer,255,f);
   fgets(buffer,255,f);
   int nPts = atoi(buffer);
   
   for(int i=0;i<nPts;i++) {
      double v[4],w[4];
            
      fgets(buffer,255,f);
      char *p = strtok(buffer," \r\n\t");
      int j = -1;
      while(p && j<4) {
         v[++j] = atof(p);
         p = strtok(NULL," \r\n\t");
      }

      fgets(buffer,255,f);
      p = strtok(buffer," \r\n\t");
      j = -1;
      while(p && j<4) {
         w[++j] = atof(p);
         p = strtok(NULL," \r\n\t");
      }

      Vector2D x1(v[0],v[1]),n1(v[2],v[3]);
      Vector2D x2(w[0],w[1]),n2(w[2],w[3]);
      
      CBoundlet cb;
      segments.push_back(cb);
      segments[i].compute(x1,x2,n1,n2);        
   }
   
   fclose(f);
   
   return 0;
}


// Save the boundary representation from a text file (2 or 4 doubles per line),
int PHSpline::save(char *fname) {
   FILE *f = fopen(fname,"wt");
   if(!f)
      return -1;
   
   fprintf(f,"#Spline Shape Description\n%d\n",segments.size());
   for(int i=0;i<segments.size();i++) {
      Vector2D x,n;
      segments[i].getInterpolation(0,x,n);
      fprintf(f,"%lg %lg %lg %lg \n",x.x,x.y,-n.y,n.x);
      segments[i].getInterpolation(1,x,n);
      fprintf(f,"%lg %lg %lg %lg \n",x.x,x.y,-n.y,n.x);
   }
   
   fclose(f);
   
   return 0;
}

/*
// Find all local maxima of curvature and store them in a pair of vectors
int PHSpline::getCurvatureMaxima(vector<double> &points,vector<double> &curv) {
   vector<double> p;
   vector<double> c;

   for(int i=0;i<segments.size();i++) {
   
      double tMax = 0;
      double vMax = 0;
      for(double t=0;t<1.0;t+=0.001) {
         double v = segments[i].getCurvature(t);
         if(v > vMax) {
            vMax = v;
            tMax = t;
         }
      }

      p.push_back(tMax);
      c.push_back(-vMax);
   }
*/
/*
      Vector x;
      double val;
      
      // Create a linear minimizn problem
      PHCurvaturePbm pbm(segments[i]);
      // BrentLinearMethod blm(pbm,Vector(1,0.5),Vector(1,1.0));
      // blm.run();
      // x = blm.x;
      // val = blm.value;

      ConjugateGradientMethod cgm(NumericalFunction(pbm),Vector(1,0.5));
      while(!cgm.isFinished()) {
         cgm.performIteration();
      }
      x = cgm.getBestEverX();
      val = cgm.getBestEverValue();


      // OK, now get the value and new x
      double T = x(0); // blm.x(0);

      p.push_back(T);
      c.push_back(val);
   }
*/
/*
	//fclose(f);


   // Clear output arrays
   points.clear();
   curv.clear();

   // Now remove all the boundary maxima
   for(i=0;i<segments.size();i++) {
      if(p[i] >= 1.0 || p[i] <= 0.0) 
			continue;

      points.push_back(p[i] + i);
      curv.push_back(c[i]);
   }

   return 0;
}
*/

// Construct a spline from a list of point coordinates and normals 
void CBSpline::buildFromPB(PolygonBoundary &pb,int numSplines) {	
	
	// Drop the spline count
	if(numSplines > pb.vertices.size()-2) {
		numSplines = pb.vertices.size() - 2;
	}

	// Create a new spline
	reset();
	spline = new BSpline1D(numSplines+2,2,4);

	// Create a point matrix
	Matrix Q(pb.vertices.size()+1,2);
	for(int i=0;i<pb.vertices.size();i++) {
		Q(i,0) = pb.vertices[i].x.x;
		Q(i,1) = pb.vertices[i].x.y;
	}
	Q(pb.vertices.size(),0) = Q(0,0);
	Q(pb.vertices.size(),1) = Q(0,1);

	// Fit the spline
	// TODO: Fix Me!!! 
    // spline->fitToPoints(Q);
}

void CBSpline::buildFromPB(PolygonBoundary &pb) {
	buildFromPB(pb,pb.vertices.size() / 4);
}

// Get the position along the spline at given t (t ranges from 0 to 1) with same
// point corresponding to t=0 and t=1
void CBSpline::getInterpolation(double t,Vector2D &x) {
	t -= floor(t);
	if(spline) {
		float xy[2];
		int k = spline->kv.getKnotAtParm(t);
		MySMLVec4f W[4];
		spline->kv.basisJet(k,t,W);
		spline->interpolatePoint(k-3,W,0,0,1,xy);
		x.x = xy[0];
		x.y = xy[1];
	}
}
void CBSpline::getInterpolation(double t,Vector2D &x,Vector2D &n) {
	t -= floor(t);
	if(spline) {
		float xy[2],nxy[2];
		int k = spline->kv.getKnotAtParm(t);
		MySMLVec4f W[4];
		spline->kv.basisJet(k,t,W);
		spline->interpolatePoint(k-3,W,0,0,1,xy);
		spline->interpolatePoint(k-3,W,1,0,1,nxy);
		x.x = xy[0];
		x.y = xy[1];
		n.x = -nxy[1];
		n.y = nxy[0];
	}
}
void CBSpline::getInterpolationJet(double t,Vector2D &x,Vector2D &x1,Vector2D &x2) {
	t = fmod(t,tmax());
	if(spline) {
		float xy[2],xy1[2],xy2[2];
		int k = spline->kv.getKnotAtParm(t);
		MySMLVec4f W[4];
		spline->kv.basisJet(k,t,W);
		spline->interpolatePoint(k-3,W,0,0,1,xy);
		spline->interpolatePoint(k-3,W,1,0,1,xy1);
		spline->interpolatePoint(k-3,W,2,0,1,xy2);
		x.x = xy[0];
		x.y = xy[1];
		x1.x = xy1[0];
		x1.y = xy1[1];
		x2.x = xy2[0];
		x2.y = xy2[1];
	}
}

int CBSpline::load(char *fname) {
	return 0;
}

int CBSpline::save(char *fname) {
	return 0;
}




// End namespace
NAMESPACE_PAULY_END
