#ifndef __TracerMainWindow_h_
#define __TracerMainWindow_h_

#include "Trackball.h"
#include <FL/Fl_Gl_Window.h>
#include "TracerData.h"

class TracerMainWindow 
: public Fl_Gl_Window,
  virtual public ITracerDataListener
{
public:
  // Modes and button information
  enum EditMode { TRACKBALL, TRACER };
  enum TrackballMode { NONE, ROTATE, ZOOM, PAN };

  // Edge display mode (for edges displayed over the surface of 
  // the mesh
  enum EdgeDisplayMode 
    { EDGE_DISPLAY_NONE, EDGE_DISPLAY_PLAIN, 
      EDGE_DISPLAY_LENGTH, EDGE_DISPLAY_DISTANCE };

  // Surface display mode
  enum SurfaceDisplayMode 
    { SURFACE_DISPLAY_ALL, SURFACE_DISPLAY_NEIGHBORHOOD };

  // Constructor
  TracerMainWindow(int x, int y, int w, int h, const char *label);
  
  // Destructor
  virtual ~TracerMainWindow();

  // Associate the data with the window
  void SetTracerData(TracerData *data) 
    {
    // Detach from the old data
    if(m_Data) m_Data->RemoveTracerDataListener(this);

    // Get the new data pointer
    m_Data = data;

    // Reset state
    m_DisplayListDirty = true;
    m_CurrentPoint = -1;

    // Attach to new data
    m_Data->AddTracerDataListener(this);
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

  // Set the surface display mode
  void SetSurfaceDisplayMode(SurfaceDisplayMode mode)
    {
    m_SurfaceDisplayMode = mode;
    m_DisplayListDirty = true;
    m_EdgeDisplayListDirty = true;
    }

  // Draw method
  void draw();

  // Event handler
  int handle(int);

  // Callbacks from TracerData event system
  void OnMeshChange(TracerDataEvent *evt);
  void OnFocusCurveChange(TracerDataEvent *evt);
  void OnFocusPointChange(TracerDataEvent *evt);
  void OnFocusCurveDataChange(TracerDataEvent *evt);
  void OnCurveListChange(TracerDataEvent *evt);
  void OnEdgeWeightsUpdate(TracerDataEvent *evt);

private:

  // Display list associated with the brain surface and with the
  // overlay edges on the brain surface
  int m_DisplayList, m_EdgeDisplayList;

  // Whether the display lists requires recomputation
  bool m_DisplayListDirty, m_EdgeDisplayListDirty;

  // Whether GL state needs reinitialization
  bool m_GLStateDirty;

  // Current mode for edge overlay display
  EdgeDisplayMode m_EdgeDisplayMode;

  // Current mode for surface display
  SurfaceDisplayMode m_SurfaceDisplayMode;

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

  // Do the actual job of computing the display list
  void ComputeDisplayList();
  void ComputeEdgeDisplayList();

  // Methods for choosing edge colors based on values
  void SetGLColorHSV(double xHue, double xSaturation, double xValue);
  void SetGLEdgeColorFromDistance(double xDistance);
  void SetGLEdgeColorFromWeight(double xWeight);

};


#endif

