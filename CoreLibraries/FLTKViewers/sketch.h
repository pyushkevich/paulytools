#ifndef _SKETCH_
#define _SKETCH_

#include <Fl/gl.h>
#include <Fl/Fl_Gl_Window.h>

#include <vector>
using namespace std;

#include <optima.h>
#include <phspline.h>

class SketchUI;

struct SketchPoint {
   Vector x;
   double t;
};

class SketchTool {
   Vector getSketchPoint(double t);

public:
   // A list of points that are in the current sketch
   // vector <Vector> sketch;
   // vector <Vector> outline;
	PolygonBoundary sketch;
	CBSpline outline;

   // User interface
   SketchUI *ui;

   void start(double x,double y);
   void end(double x,double y);
   void plot(double x,double y);

   // Current mode of the tool
   bool sketchMode;

};

void getUserSketch(ContinuousBoundary &cb);

#endif
