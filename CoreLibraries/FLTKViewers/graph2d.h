#ifndef _GRAPH_2D_
#define _GRAPH_2D_

#include <Fl/gl.h>
#include <Fl/Fl_Gl_Window.h>

#include <vector>
using namespace std;

#include <optima.h>

class Graph2DUI;

// A graph state.  The graph goes through series of states
// For instance each time you zoom in, you save a state and enable a new state
struct Graph2DState {
   // Domain of the function
   Vector xMin,xMax;

   // Resolution of the graph
   int res;

   // Values of the function
   Matrix F;

   // Max and min values of the function
   double vMax,vMin;
};

class Graph2D {
public:
   // Constructor
   Graph2D(Function &function,Vector min,Vector max,int res);

   // Destructor
   ~Graph2D();

   // Called to compute the new frame
   void compute();

   // Zoom out,in
   void zoomOut();
   void zoomIn();

   // Move to next state
   void next() {
      dassert(state < states.size()-1);
      setActiveState(state+1);
   }

   // Move to prev state
   void back() {
      dassert(state > 0);
      setActiveState(state-1);
   }

   // Show axis
   void showAxis(bool mode);

   // Show grid
   void showGrid(bool mode);

   // Show grid labels
   void showLabels(bool mode);

   // Repaint after colors have changed
   void repaint() {
      paintPixels();
   }

	// Set a callback function for when a mark is made on the graph
	void setMarkCallback(void (*i)(bool,double,double)) {
		markCallback = i;
	}

	// Redraw
	void redraw();

	virtual void drawChild() {};

private:
   // Function to graph
   Function &f;

   // A stack of graph states
   vector<Graph2DState> states;

   // The current state
   int state;
   
   // Are the axis visible on the graph?
   bool axesVisible;
   
   // Is the grid visible on the graph?
   bool gridVisible;

   // Is the grid labeled on the graph?
   bool gridLabeled;

   // Are we dragging right now?
   bool dragMode;

   // The FLTK user interface 
   Graph2DUI *ui;

   // A set of pixels
   struct rgb {
      unsigned char r,g,b;
   } *pixels;

   // Add a state to the set of states
   void createState(Vector min,Vector max,int res,Graph2DState &s);

   // Display current state
   void setActiveState(int s);

   // Paint the pixels
   void paintPixels();

   // Zoom to a box
   void zoom(double x0,double y0,double x1,double y1);

   // Set a mark (call the callback routine if provided...
   void (*markCallback)(bool,double,double);

	void mark(bool visible,double x=0,double y=0); 
	


   friend class Graph2DWindow;
};

#endif
