#ifndef _DSL_WINDOWS_H_
#define _DSL_WINDOWS_H_

/******************************************************************
 * DSLTOOL                                                        *
 ******************************************************************
 * Author:                 Paul Yushkevich
 *
 * Date:                   July 18, 1999
 *
 * Description             Window handlers for DSLTool
 *									
 *	Sources:                
 *									
 * Dependencies:           matrix,dsl2d,fltk,cimage
 ******************************************************************
 * windows.h
 *	----------
 * Classes that handle display and user input for all custom windows
 * in the application
 ******************************************************************/
#include <matrix.h>
#include <mreps2D.h>
#include <Fl/Fl_Gl_Window.H>
#include <Fl/gl.h>
#include "ispace.h"

class ToolHandler;
class DisplayWindow;
class SplineObject;
class RegularSplineSample;

/******************************************************************
 * Main Display Window                                            *
 * -------------------                                            *
 * In this window you see the image, the model and everything     *
 * related.                                                       *
 ******************************************************************/
class DisplayWindow : public Fl_Gl_Window {
public:
	void drawSpline(RegularSplineSample *spline);
	// GL draw method
	void draw();

	// GL event handling method
	int handle(int);

	// A texture number for the image to display
	GLuint glTextureName;
	GLint glDLName;

	// Position and size of viewing box
	Vector2D winPosition;
	double winSize;
	
	// Are we currently dragging?
	bool dragging;

   // Do we need to recompute texture
   bool needTextureUpdate;
	
	// Current tool
	ToolHandler *tool;

	// The selection box
	bool haveSelection;
	Vector2D selTopLeft,selBotRight;

	// Convert between space coordinates and screen coordinates
	Vector2D screenToSpace(double x,double y);	
	Vector2D spaceToScreen(double x,double y);	

	// Constructor
	DisplayWindow(int X, int Y, int W, int H, const char *L);

	// Set a new tool (give me a pointer to the tool, it will be destroyed later)
	void setTool(ToolHandler *tool);

	// Build a display list for the image (call whenever image changes)
	void drawImage();
	void drawGrid();
	void createTexture();
	void updateTexture() {
		needTextureUpdate = true;
	}

	// Build a display list for the model (call whenever model changes)
	//void updateModelDL();

	// Build a display list for the selection (call whenever selection changes)
	//void updateSelectionDL();


	friend class ToolHandler;
};


/******************************************************************
 * Scale Space Window                                             *
 * ------------------                                             *
 * This window displays the image blurred with Gaussians of       *
 * variuous aperture.                                             *
 ******************************************************************/
class ScaleSpaceWindow : public Fl_Gl_Window {
private:
   // Scale amount
   double scale;

   // Pixel zoom amount
   int zoom;

   // Size of it
   int iw,ih;

public:
   // Image to display
   Array2D<GLbyte> image;

	// GL draw method
	void draw();

	// GL event handling method
	int handle(int);

	// Constructor
   ScaleSpaceWindow(int X, int Y, int W, int H, const char *L=NULL) : Fl_Gl_Window(X,Y,W,H,L) {
      iw = W-W%4;
      ih = H-H%4;
      image.resize(iw,ih);
      zoom = 1;
      setScale(0.0);
   }

   ~ScaleSpaceWindow() {
   }

   // Set the scale
   void setScale(double newScale);

   // Set the zoom (to speed up rendering)
   void setZoom(int newZoom);
};


/******************************************************************
 * Intensity Profile Window                                       *
 * ------------------------                                       *
 * This window displays the profile along the normal vector of    *
 * the boundariness function.                                     *
 ******************************************************************/
class IProfileWindow : public Fl_Gl_Window {
private:
   // Boundary Primitive that will be shown in the profile window
   BAtom ba;

   // Whether or not we are displaying anything
   bool display;

   // Number of samples
   int samples;

public:
	// GL draw method
	void draw();

	// GL event handling method
   int handle(int) {
      return 0;
   }

	// Constructor
   IProfileWindow(int X, int Y, int W, int H, const char *L=NULL) : Fl_Gl_Window(X,Y,W,H,L) {
      display = false;
   }

   // Set the primitive 
   void set(const BAtom &bndAtom) {
      ba = bndAtom;
      display=true;
      redraw();
   }

   // Show nothing
   void unset() {
      display = false;
      redraw();
   }
};


/******************************************************************
 * Intensity Gradient Window                                      *
 * -------------------------                                      *
 * This window displays the direction of the image gradient at    *
 * a boundary site.                                                *
 ******************************************************************/
class IGradientWindow : public Fl_Gl_Window {
private:
   // Boundary Primitive that will be shown in the profile window
   BAtom ba;

   // Whether or not we are displaying anything
   bool display;

public:
	// GL draw method
	void draw();

	// GL event handling method
   int handle(int) {
      return 0;
   }

	// Constructor
   IGradientWindow(int X, int Y, int W, int H, const char *L=NULL) : Fl_Gl_Window(X,Y,W,H,L) {
      display = false;
   }

   // Set the primitive 
   void set(const BAtom &bndAtom) {
      ba = bndAtom;
      display=true;
      redraw();
   }

   // Show nothing
   void unset() {
      display = false;
      redraw();
   }

   // Settings for the display window
   void setPixMapType(int type);
};

// Some constants 
const int MAGIMAGE_PIXELS = 0;
const int MAGIMAGE_BLURRED = 1;
const int MAGIMAGE_GRADIENT = 2;
const int MAGIMAGE_DIRDERIV = 3;
const int MAGIMAGE_POWERB = 4;
const int MAGIMAGE_NOIMAGE = 5;






#endif // _DSL_WINDOWS_H_
