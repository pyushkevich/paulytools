#ifndef __TracerMainWindow_h_
#define __TracerMainWindow_h_

#include <FL/Fl_Gl_Window.h>
#include "TracerData.h"
#include "Trackball.h"

class TracerMainWindow : public Fl_Gl_Window
{
public:
  // Constructor
  TracerMainWindow(int x, int y, int w, int h, const char *label) 
    : Fl_Gl_Window(x,y,w,h,label) 
    {
    m_Data = NULL;
    m_DisplayListDirty = false;
    m_DisplayList = -1;
    } 

  // Destructor
  virtual ~TracerMainWindow() {};

  // Associate the data with the window
  void SetTracerData(TracerData *data) 
    {
    m_Data = data;
    m_DisplayListDirty = true;
    }

  // Request an update of the display lists
  void OnMeshUpdate() 
    {
    m_DisplayListDirty = true;
    if(shown()) redraw();
    }

  // Draw method
  void draw();

  // Event handler
  int handle(int);

private:

  // Do the actual job of computing the display list
  void ComputeDisplayList();

  // Display list associated with the brain surface
  int m_DisplayList;
  
  // Whether the display lists requires recomputation
  bool m_DisplayListDirty;

  // Pointer to the tracer data
  TracerData *m_Data;

  // 3D trackball
  Trackball m_Trackball;
};


#endif

