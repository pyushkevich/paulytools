#ifndef _DISPLAY_WINDOW_
#define _DISPLAY_WINDOW_

#include <Fl/Fl_Gl_Window.h>

#include <matrix.h>
// #include "DSL2DLib/Src/DSL2D.h"
#include <IMAGE2D.h>

class DisplayWindow;
class ToolHandler;

/**
 * A small helper event object to be passed on to tool handlers
 */
struct MouseEvent {
	// Starting mouse position and ending mouse position in space
	Vector sSpace;
	Vector xSpace;

	// Starting mouse position and ending mouse position on the screen
	Vector sScreen;
	Vector xScreen;

	// Target figure and primitive
	int targetFigure;
	int targetPrimitive;

	// Target within the primitive: center,site1,site2
	static const int CENTER,SITE1,SITE2;
	int targetSite;

	// Shift/Ctrl/Alt flags
	bool shift,ctrl,alt;

	// Button tyoe
	int button;

	// Drag length - total distance that the user has dragged the mouse
	double dragLength;

	// Click flag - set to true if dragLength is less than epsilon
	bool click;

	// Default constructor
	MouseEvent() {
		sSpace = Vector(2);
		xSpace = Vector(2);
		sScreen = Vector(2);
		xScreen = Vector(2);
		targetFigure = targetPrimitive = -1;
		shift = ctrl = alt = false;
		dragLength = 0;
		click = false;
	}
};





/**
 * A window in which we will be performing all sorts of GL functions.
 *
 */
class DisplayWindow : public Fl_Gl_Window {
public:
	// GL draw method
	void draw();

	// GL event handling method
	int handle(int);

	// A block of display lists
	int dlImage;//,dlModel,dlSelection;

	// Position and size of viewing box
	Vector winPosition;
	double winSize;
	
	// Are we currently dragging?
	bool dragging;

   // Are we in primitive edit mode
   bool pemMode;
	
	// Current tool
	ToolHandler *tool;

	// The selection box
	bool haveSelection;
	Vector selTopLeft,selBotRight;

	// Convert between space coordinates and screen coordinates
	Vector screenToSpace(double x,double y);	

	// Constructor
	DisplayWindow(int X, int Y, int W, int H, const char *L);

	// Set a new tool (give me a pointer to the tool, it will be destroyed later)
	void setTool(ToolHandler *tool);

	// Build a display list for the image (call whenever image changes)
	void updateImageDL();

	// Build a display list for the model (call whenever model changes)
	//void updateModelDL();

	// Build a display list for the selection (call whenever selection changes)
	//void updateSelectionDL();


	friend class ToolHandler;
};


/**
 * A window for displaying the scale space
 */
class ScaleSpaceWindow : public Fl_Gl_Window {
private:
   // Scale amount
   double scale;

   // Pixel zoom amount
   int zoom;

   // Image to display
   unsigned char *image;

   // Size of it
   int iw,ih;

public:
	// GL draw method
	void draw();

	// GL event handling method
	int handle(int);

	// Constructor
   ScaleSpaceWindow(int X, int Y, int W, int H, const char *L=NULL) : Fl_Gl_Window(X,Y,W,H,L) {
      image = new unsigned char[W*H];
      iw = W-W%4;
      ih = H-H%4;
      zoom = 1;
      setScale(0.0);
   }

   ~ScaleSpaceWindow() {
      delete image;
   }

   // Set the scale
   void setScale(double newScale);

   // Set the zoom (to speed up rendering)
   void setZoom(int newZoom);
};


/**
 * A window for displaying the intensity profile
 */
class IProfileWindow : public Fl_Gl_Window {
private:
   // Boundary Primitive that will be shown in the profile window
   double x, y, n, sigma;

   // Whether or not we are displaying anything
   bool display;

   // Number of samples
   const int samples;

public:
	// GL draw method
	void draw();

	// GL event handling method
   int handle(int) {
      return 0;
   }

	// Constructor
   IProfileWindow(int X, int Y, int W, int H, const char *L=NULL) : Fl_Gl_Window(X,Y,W,H,L),samples(40) {
      display = false;
   }

   // Set the primitive 
   void set(double x,double y,double n,double sigma) {
      this->x = x;
      this->y = y;
      this->n = n;
      this->sigma = sigma;
      display=true;
      redraw();
   }

   // Show nothing
   void unset() {
      display = false;
      redraw();
   }
};


#endif
