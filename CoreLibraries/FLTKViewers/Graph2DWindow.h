#ifndef _GRAPH_2D_WINDOW_
#define _GRAPH_2D_WINDOW_

#include <Fl/gl.h>
#include <Fl/Fl_Gl_Window.h>

#include <matrix.h>
class Graph2DUI;

/******************************************************
 * A GL Window for FLTK that we will use
 ******************************************************/

class Graph2DWindow : public Fl_Gl_Window {    
   void draw();    
   int handle(int);  

   // Dragging currently?
   bool dragging;

   // Drag start coordinates (in viewport coords);
   double dragX0,dragY0;
   double dragX1,dragY1;

   // Special mark position (in viewport coords);
   bool markVisible;
   double markX,markY;

public:    
   // Constructor
   Graph2DWindow(int X, int Y, int W, int H, const char* L=NULL)
      : Fl_Gl_Window(X,Y,W,H,L) 
   { 
      ui = NULL;
      dragging = false;
      markVisible = false;
   }
   

   // Pointer to the parent UI structure
   Graph2DUI *ui;

   friend class Graph2D;
};  


class Graph2DPalette : public Fl_Gl_Window {    
   void draw();    
   int handle(int);  

   // Windowing points for negative and positive scales
   Matrix M;

   // Convert from function to intensity (range -1 to 1)
   double getIntensity(double x);


public:    
   // Constructor
   Graph2DPalette(int X, int Y, int W, int H, const char* L=NULL)
      : Fl_Gl_Window(X,Y,W,H,L) 
   { 
      M.setSize(2,2);
      M(0,0) = -1;
      M(0,1) = 0;
      M(1,0) = 1;
      M(1,1) = 0;

      ui = NULL;
   }

   void reset();
   
   Graph2DUI *ui;
   friend class Graph2D;
};  


#endif
