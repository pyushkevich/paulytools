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

  // Edge display mode (for edges displayed over the surface of 
  // the mesh
  enum EdgeDisplayMode { EDGE_DISPLAY_NONE, 
    EDGE_DISPLAY_PLAIN, EDGE_DISPLAY_LENGTH, EDGE_DISPLAY_DISTANCE };

  // Constructor
  TracerMainWindow(int x, int y, int w, int h, const char *label);
  
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

    // The edge drawing is no longer clean if drawing accumulated distances
    if(m_EdgeDisplayMode == EDGE_DISPLAY_DISTANCE)
      m_EdgeDisplayListDirty = true;

    } 

  // Request an update of the display lists
  void OnMeshUpdate() 
    {
    m_DisplayListDirty = true;
    m_EdgeDisplayListDirty = true;
    if(shown()) redraw();
    }

  // Set the tracing/tracking mode
  void SetMode(EditMode mode) 
    {
    m_EditMode = mode;
    redraw();
    }

  // Set the edge display mode
  void SetEdgeDisplayMode(EdgeDisplayMode mode)
    {
    m_EdgeDisplayMode = mode;
    m_EdgeDisplayListDirty = true;
    redraw();
    }

  // Draw method
  void draw();

  // Event handler
  int handle(int);

private:

  // Do the actual job of computing the display list
  void ComputeDisplayList();
  void ComputeEdgeDisplayList();

  // Methods for choosing edge colors based on values
  void SetGLColorHSV(double xHue, double xSaturation, double xValue);
  void SetGLEdgeColorFromDistance(double xDistance);
  void SetGLEdgeColorFromWeight(double xWeight);

  // Display list associated with the brain surface and with the
  // overlay edges on the brain surface
  int m_DisplayList, m_EdgeDisplayList;

  // Whether the display lists requires recomputation
  bool m_DisplayListDirty, m_EdgeDisplayListDirty;

  // Whether GL state needs reinitialization
  bool m_GLStateDirty;

  // Current mode for edge overlay display
  EdgeDisplayMode m_EdgeDisplayMode;

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

