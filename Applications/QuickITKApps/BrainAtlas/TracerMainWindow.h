#ifndef __TracerMainWindow_h_
#define __TracerMainWindow_h_

#include "Trackball.h"
#include <FL/Fl_Gl_Window.h>
#include "TracerData.h"

class TracerMainWindow 
: public Fl_Gl_Window,
  virtual public ITracerDataListener,
  virtual public VTKMeshShortestDistance::ICellChecher
{
public:
  // Modes and button information
  enum EditMode { TRACKBALL, TRACER, MARKER };
  enum TrackballMode { NONE, ROTATE, ZOOM, PAN };

  // Edge display mode (for edges displayed over the surface of 
  // the mesh
  enum EdgeDisplayMode 
    { EDGE_DISPLAY_NONE, EDGE_DISPLAY_PLAIN, 
      EDGE_DISPLAY_LENGTH, EDGE_DISPLAY_DISTANCE };

  // Surface display mode
  enum SurfaceDisplayMode 
    { SURFACE_DISPLAY_ALL, SURFACE_DISPLAY_NEIGHBORHOOD };

  // Vector type
  typedef TracerData::Vec Vec;

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
    m_FullMeshDisplayListDirty = m_NeighborhoodDisplayListDirty = 
      m_EdgeDisplayListDirty = true;
    
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
    m_EdgeDisplayListDirty = true;
    redraw();
    }

  // Set the neighborhood size for surface display
  void SetNeighborhoodSize(double size)
    {
    m_NeighborhoodSize = size;
    if(m_SurfaceDisplayMode == SURFACE_DISPLAY_NEIGHBORHOOD)
      {
      m_NeighborhoodDisplayListDirty = true;
      m_EdgeDisplayListDirty = true;
      redraw();
      }
    }

  // Set mesh centering mode
  void SetCenterMesh(bool mode)
    {
    m_CenterMesh = mode;
    m_Trackball.ResetPan();
    redraw();
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
  void OnMarkerListChange(TracerDataEvent *evt);
  void OnFocusMarkerChange(TracerDataEvent *evt);
  void OnSegmentationChange(TracerDataEvent *evt);

private:

  // Display list associated with the brain surface and with the
  // overlay edges on the brain surface
  int m_FullMeshDisplayList;
  int m_EdgeDisplayList;
  int m_NeighborhoodDisplayList;

  // Whether the display lists requires recomputation
  bool m_FullMeshDisplayListDirty;
  bool m_EdgeDisplayListDirty;
  bool m_NeighborhoodDisplayListDirty;

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
 
  // Candidate cell for a marker
  vtkIdType m_CurrentMarkerCell;

  // Whether the mesh should be centered on current point or not
  bool m_CenterMesh;

  // Neighborhood size for neighborhood display mode
  double m_NeighborhoodSize;
  
  // Internal GL init method
  void InitializeGL();
  void SetUpModelMatrix();

  // Perform ray intersection to find point under cursor
  void ComputeClickRay(Vec &xStart, Vec &xEnd);
  void FindPointUnderCursor();
  void FindCellUnderCursor();

  // Do the actual job of computing the display list
  void ComputeFullMeshDisplayList();
  void ComputeNeighborhoodDisplayList();
  void ComputeEdgeDisplayList();

  // Submethods in the draw() method, broken up for code clarity
  void DrawMesh();
  void DrawCurves();
  void DrawMarkers();

  // Methods for choosing edge colors based on values
  void SetGLColorHSV(double xHue, double xSaturation, double xValue);
  void SetGLEdgeColorFromDistance(double xDistance);
  void SetGLEdgeColorFromWeight(double xWeight);

  // Draw a sphere using open GL
  void GLDrawSphere(double *x, double r);
  void GLDrawMarker(vtkIdType iCell, const Vec &color);
  void GLDrawStrippedPolyData(vtkPolyData *poly);

  // Callback used to accept/reject cells during ray-mesh intersection
  bool CheckCell(vtkIdType iCell);
};


#endif

