#ifndef _SKETCH_WINDOW_
#define _SKETCH_WINDOW_

#include <Fl/gl.h>
#include <Fl/Fl_Gl_Window.h>

#include <matrix.h>

class SketchTool;

/******************************************************
 * A GL Window for FLTK that we will use
 ******************************************************/

class SketchWindow : public Fl_Gl_Window {    
public:    
   void draw();    
   int handle(int);  

   // Dragging currently?
   bool dragging;

   // Drag start coordinates (in viewport coords);
   double dragX0,dragY0;
   double dragX1,dragY1;

   // Parent sketch object
   SketchTool *sketch;

   // Constructor
   SketchWindow(int X, int Y, int W, int H, const char* L=NULL)
      : Fl_Gl_Window(X,Y,W,H,L) 
   { 
      sketch = NULL;
      dragging = false;
   }
   
   friend class SketchTool;
};  


#endif
