#ifndef _MREPS_INTERP_
#define _MREPS_INTERP_

#include <mreps2D.h>
#include <optima.h>

/*************************************************************************
 A Class used to interpolate a medial primitive from boundary interpolation
 *************************************************************************/
class MIProblem : public Function {
private:
   MNode *n[2];
	Boundary *b;

	double sLeft,sRight;
	double bLenLeft,bLenRight;
   double t,help;

   // Position and direction of the line that is t-ways between start
   // and end. x(t) = r + t*u
   Vector2D rx,ru;

public:
   MIProblem(MNode *start,MNode *end,double t,double help=1.0);
   MAtom getAtom(const Vector &x);
   double evaluate(const Vector &x);

	void setHelp(double h) {
		help = h;
	}
};

// This method interpolates a selected figure and update its graph
void interpolateFigure(MFigure *figure,int n);

// This method will resample a figure at even intervals and update its graph
void resampleFigureUniform(MFigure *figure);

#endif
