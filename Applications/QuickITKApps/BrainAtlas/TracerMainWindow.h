#ifndef __TracerMainWindow_h_
#define __TracerMainWindow_h_

#include "Trackball.h"
#include <FL/Fl_Gl_Window.h>
#include "TracerData.h"

class TracerMainWindow : public Fl_Gl_Window
{
public:
  // Modes and button information
  enum EditMode { TRACKBALL, TRACER };
  enum TrackballMode { NONE, ROTATE, ZOOM, PAN };

  // Constructor
  TracerMainWindow(int x, int y, int w, int h, const char *label) 
    : Fl_Gl_Window(x,y,w,h,label) 
    {
    m_Data = NULL;
    m_DisplayListDirty = false;
    m_DisplayList = -1;
    m_EditMode = TRACKBALL;
    m_TrackballMode = NONE;
    m_GLStateDirty = true;
    m_CurrentPoint = -1;
    } 

  // Destructor
  virtual ~TracerMainWindow() {};

  // Associate the data with the window
  void SetTracerData(TracerData *data) 
    {
    m_Data = data;
    m_DisplayListDirty = true;
    m_CurrentPoint = -1;
    }
    
  // Update the currently displayed curve
  void OnCurrentCurveChange()
    {
    m_CurrentPoint = -1;
    redraw();
    } 

  // Request an update of the display lists
  void OnMeshUpdate() 
    {
    m_DisplayListDirty = true;
    if(shown()) redraw();
    }

  // Set the tracing/tracking mode
  void SetMode(EditMode mode) 
    {
    m_EditMode = mode;
    redraw();
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

  // Whether GL state needs reinitialization
  bool m_GLStateDirty;

  // Pointer to the tracer data
  TracerData *m_Data;

  // 3D trackball
  Trackball m_Trackball;

  // Mode states
  EditMode m_EditMode;
  TrackballMode m_TrackballMode;

  // Point under the cursor
  vtkIdType m_CurrentPoint;

  // Internal GL init method
  void InitializeGL();
  void SetUpModelMatrix();

  // Perform ray intersection to find point under cursor
  void FindPointUnderCursor();
};


#endif

