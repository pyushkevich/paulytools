#ifndef __TracerData_h_
#define __TracerData_h_

#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkStripper.h"
#include "vtkPolyDataNormals.h"
#include "VTKMeshShortestDistance.h"
#include "vnl/vnl_cross.h"
#include <list>
#include <ctime>

#include "TracerCurves.h"
#include "EventListenerModel.h"

class TracerData;

/** Event used to communicate changes in the tracer data to the 
  user interface model code */
class TracerDataEvent : public EventObject<TracerData>
{
public:
  typedef EventObject<TracerData> Superclass;

  TracerDataEvent(TracerData *source) : Superclass(source) {}
  virtual ~TracerDataEvent() {}
};

/** Listener for tracer data events */
class ITracerDataListener
{
public:
  /** Fired when the mesh is changed (reloaded) */
  virtual void OnMeshChange(TracerDataEvent *evt) = 0;

  /** Fired when the point (path source) is changed */
  virtual void OnFocusCurveChange(TracerDataEvent *evt) = 0;

  /** Fired when the curve under focus is changed */
  virtual void OnFocusPointChange(TracerDataEvent *evt) = 0;

  /** Fired when a change is made to the curve under focus */
  virtual void OnFocusCurveDataChange(TracerDataEvent *evt) = 0;

  /** Fired when the list of curves is changed, including number of
    * curves and the names of the curves */
  virtual void OnCurveListChange(TracerDataEvent *evt) = 0; 

  /** Fired when the edge weights are updated */
  virtual void OnEdgeWeightsUpdate(TracerDataEvent *evt) = 0;

  /** Fired when the list of markers is updated */
  virtual void OnMarkerListChange(TracerDataEvent *evt) = 0;

  /** Fired when the active marker changes */
  virtual void OnFocusMarkerChange(TracerDataEvent *evt) = 0;
};

/**
 * An encapsulation of the data used by the curve tracer. This class provides
 * a listener interface for the GUI to be alerted when curve information or mesh
 * information changes 
 */
class TracerData
{
public:
  // Typedef for the 3d vector
  typedef vnl_vector_fixed<double, 3> Vec;
  typedef vector<Vec> VecArray;
  typedef TracerCurves::IdType IdType;

  // Tracer data constructors
  TracerData();
  ~TracerData();

  // Add listener method (call AddTracerDataListener)
  pyAddListenerMacro(TracerDataListener, ITracerDataListener);
  
  // Remove listener method (call AddTracerDataListener)
  pyRemoveListenerMacro(TracerDataListener, ITracerDataListener);

  // Load the mesh data from disk
  void LoadInputMesh(const char *file);
  
  // Save the curves to disk
  void SaveCurves(const char *file);

  // Load curves from disk
  bool LoadCurves(const char *file);

  // Check if the mesh is loaded
  bool IsMeshLoaded() 
    { return m_Mesh != NULL; }

  // Get the 'source' mesh
  vtkPolyData *GetInternalMesh() 
    { return m_Mesh; }

  // Get the mesh used for display
  vtkPolyData *GetDisplayMesh() 
    { return m_DisplayMesh; }

  /** Get the distance mapper */
  const VTKMeshShortestDistance *GetDistanceMapper()
    { return m_DistanceMapper; }

  // Get the coordinate for a given point index (in internal mesh)
  void GetPointCoordinate(vtkIdType id, Vec &target)
    { m_Mesh->GetPoint(id,target.data_block()); }

  // Get the coordinate and normal for a given point index (in internal mesh)
  void GetPointCoordinateAndNormal(vtkIdType id, Vec &x, Vec &n)
    {
    m_Mesh->GetPoint(id,x.data_block());
    n.set( m_Mesh-> GetPointData()->GetNormals()->GetTuple3(id));
    }

  // Create a new curve (given a name)
  void AddNewCurve(const char *name)
    {
    SetFocusPoint(NO_FOCUS);
    SetFocusCurve(m_Curves.AddCurve(name));
    }
    
  // Set the current curve
  void SetCurrentCurve(IdType iCurve)
    {
    // Check that the curve exists
    assert(m_Curves.IsCurvePresent(iCurve));

    // Set the focus point
    SetFocusPoint( m_Curves.GetCurveControls(iCurve).size()
      ? (int) m_Curves.GetCurveControls(iCurve).back() : NO_FOCUS );

    // Set the focus curve and notify listeners
    SetFocusCurve(iCurve);
    }

  const TracerCurves *GetCurves()
    { return &m_Curves; }

  void SetCurveName(IdType iCurve, const char *name)
    {
    // Curve must exist
    assert(m_Curves.IsCurvePresent(iCurve));

    // Update the curve's name
    m_Curves.SetCurveName(iCurve, name);

    // Fire appropriate event
    TracerDataEvent evt(this);
    BroadcastOnCurveListChange(&evt);
    }

  /** Delete the last point on the current curve */
  void DeleteCurrentPoint()
    {
    // The current point must be the last point in a curve
    assert(m_FocusPoint != NO_FOCUS && m_FocusCurve != NO_FOCUS);
    assert(m_FocusPoint == 
      m_Curves.GetCurveControls(m_FocusCurve).back());

    // Delete the last point in the curve
    m_Curves.DeleteLastControlPointInCurve(m_FocusCurve);

    // Set the focus point if the curve isn't empty
    SetFocusPoint( m_Curves.GetCurveControls(m_FocusCurve).size()
      ? (int) m_Curves.GetCurveControls(m_FocusCurve).back() : NO_FOCUS );

    // Nofity listeners that the current curve's data has changed
    TracerDataEvent evt(this);
    BroadcastOnFocusCurveDataChange(&evt);
    }
    
  void DeleteCurrentCurve() 
    {
    // We need to have a current curve
    assert(m_FocusCurve != NO_FOCUS);

    // Delete the current curve
    m_Curves.DeleteCurve(m_FocusCurve);

    // Don't focus on any other curve
    SetFocusPoint(NO_FOCUS);
    SetFocusCurve(NO_FOCUS);

    // Notify listeners that the list of curves has changed
    TracerDataEvent evt(this);
    BroadcastOnCurveListChange(&evt);
    }

  /** Is a curve currently selected, or is there no current curve? */
  bool IsCurrentCurveValid()
    { return m_FocusCurve != NO_FOCUS; }
   
  IdType GetCurrentCurve()
    { return m_FocusCurve; }

  unsigned int GetNumberOfCurvePoints(IdType iCurve)
    {
    assert(m_Curves.IsCurvePresent(iCurve));
    return m_Curves.GetCurveVertices(iCurve).size();
    }

  // Create a path from the last point on the focused curve (if any)
  // to the specified point. The path includes the specified point, but
  // not the last point on the focused curve
  bool GetPathToPoint(vtkIdType iPoint, TracerCurves::MeshCurve &target)
    {
    // Clear the return array
    target.clear();

    // Do we have a 'focus' point
    if(m_FocusPoint != NO_FOCUS && m_DistanceMapper->IsVertexConnected(iPoint))
      {
      vtkIdType iLast = m_Curves.GetControlPointVertex(m_FocusPoint);
      vtkIdType iCurrent = iPoint;
      while(iCurrent != iLast)
        {
        target.push_front(iCurrent);
        iCurrent = m_DistanceMapper->GetVertexPredecessor(iCurrent);
        if(iCurrent == target.front())
          {
          cout << "No path exists from " << iPoint << " to " << iLast << 
            " (distance = " << m_DistanceMapper->GetVertexDistance(iPoint) << " ) "
            << endl;
          return false;
          }
        }
      }
    else
      {
      target.push_front(iPoint);
      }

    return true;
    }

  // Check if there is a path source
  bool IsPathSourceSet()
    { return m_FocusPoint != NO_FOCUS; }

  // Get the vertex to which the path currently extends
  vtkIdType GetPathSource()
    { return m_Curves.GetControlPointVertex(m_FocusPoint); }

  // Add a point to the current curve
  void AddNewPoint(vtkIdType iPoint)
    {
    assert(m_FocusCurve >= 0);
    
    // Add the point to the curve info
    IdType iNewControl = m_Curves.AddControlPoint(iPoint);

    // If this is not the first point in the curve, add a link
    if(m_FocusPoint != NO_FOCUS)
      {
      // Get the path to the current point
      TracerCurves::MeshCurve lPoints;
      GetPathToPoint(iPoint, lPoints);
    
      // Remove the first point in the list
      lPoints.pop_front();

      // Add the link
      m_Curves.AddLink(m_FocusPoint, iNewControl, m_FocusCurve, lPoints);
      }

    // Save the added point as the focus point
    SetFocusPoint(iNewControl);

    // Nofity of the change in curve data
    TracerDataEvent evt(this);
    BroadcastOnFocusCurveDataChange(&evt);
    }

  /** Use the Euclidean edge weight function */
  void SetEdgeWeightsToEuclideanDistance()
    { UpdateEdgeWeightFunction(new EuclideanDistanceMeshEdgeWeightFunction); }

  /** Use the Pitch edge weight function with pitch weight factor k */
  void SetEdgeWeightsToPitchDistance(double xPitchFactor)
    {
    PitchBasedMeshEdgeWeightFunction *fnNew = 
      new PitchBasedMeshEdgeWeightFunction();
    fnNew->SetPitchFactor(xPitchFactor);
    UpdateEdgeWeightFunction(fnNew);
    }

  /** Use currect edge weight distance to compute distance between a pair
   * of vertices */
  double GetMeshEdgeWeight(vtkIdType x1, vtkIdType x2)
    { return m_DistanceMapper->GetEdgeWeight(x1, x2); }

  /** Add a new marker with a particular name */
  void AddMarker(vtkIdType idCell, const char *name, 
    double r, double g, double b)
    { 
    // Add the marker
    m_Markers[idCell].Name = name; 
    m_Markers[idCell].Color = Vec(r, g, b);

    // Update the active marker
    m_FocusMarker = idCell;

    // Nofity of the change in curve data
    TracerDataEvent evt(this);
    BroadcastOnMarkerListChange(&evt);
    }

  /** Delete the marker */
  void DeleteMarker(vtkIdType idCell)
    { 
    // Remove the cell
    m_Markers.erase(idCell);

    // Change the active marker
    m_FocusMarker = NO_FOCUS;

    // Nofity of the change in curve data
    TracerDataEvent evt(this);
    BroadcastOnMarkerListChange(&evt);
    }

  /** Rename the marker */
  void RenameMarker(vtkIdType idCell, const char *name)
    { 
    m_Markers[idCell].Name = name;

    // Nofity of the change in curve data
    TracerDataEvent evt(this);
    BroadcastOnMarkerListChange(&evt);
    }

  /** Recolor the marker */
  void RecolorMarker(vtkIdType idCell, double r, double g, double b)
    { 
    m_Markers[idCell].Color = Vec(r, g, b);

    // Nofity of the change in curve data
    TracerDataEvent evt(this);
    BroadcastOnMarkerListChange(&evt);
    }

  /** Change the active marker */
  void SetCurrentMarker(vtkIdType idCell) 
    {
    // Set the focus
    m_FocusMarker = idCell;

    // Fire an event
    TracerDataEvent evt(this);
    BroadcastOnFocusMarkerChange(&evt);
    }

  /** Check if there is an active marker */
  bool IsCurrentMarkerValid() const
    { return m_FocusMarker != NO_FOCUS; }

  /** Get the active marker */
  vtkIdType GetCurrentMarker() const
    { return m_FocusMarker; }
  
  /** Get the number of markers */
  void GetMarkerIds(list<vtkIdType> &outList) const
    { 
    MarkerMap::const_iterator it = m_Markers.begin();
    while(it != m_Markers.end())
      outList.push_back(it++->first);
    }

  /** Get the corners of a marker */
  void GetMarkerVertices(vtkIdType id, Vec &x1, Vec &x2, Vec &x3) const
    {
    vtkIdType nPoints, *iPoints;
    m_Mesh->GetCellPoints(id, nPoints, iPoints);
    
    x1.set(m_Mesh->GetPoint(iPoints[0]));
    x2.set(m_Mesh->GetPoint(iPoints[1]));
    x3.set(m_Mesh->GetPoint(iPoints[2]));
    }

  /** Get the corners, the normal and the center of a marker */
  void GetMarkerGeometry(vtkIdType id, 
    Vec &x1, Vec &x2, Vec &x3, Vec &xNormal, Vec &xCenter) const
    {
    // Get the vertices
    GetMarkerVertices(id, x1, x2, x3);

    // Compute the normal using cross-product
    xNormal = vnl_cross_3d(x2-x1,x3-x1).normalize();

    // Compute the center of the triangle
    xCenter = 0.33333333 * (x1 + x2 + x3);
    }

  /** Get the color of a marker */
  Vec GetMarkerColor(vtkIdType id) const
    { 
    assert(m_Markers.find(id) != m_Markers.end());
    return m_Markers.find(id)->second.Color; 
    }
    
  /** Get the name of a marker */
  const char * GetMarkerName(vtkIdType id) const
    { 
    assert(m_Markers.find(id) != m_Markers.end());
    return m_Markers.find(id)->second.Name.c_str(); 
    }
    
private:
  // Shortest distance computer
  VTKMeshShortestDistance *m_DistanceMapper;

  // The edge weight function object
  MeshEdgeWeightFunction *m_EdgeWeightFunction;

  // The internal mesh (used to find shortest paths)
  vtkPolyData *m_Mesh;

  // The displayed mesh (triangle stripped, used for rendering)
  vtkPolyData *m_DisplayMesh;

  // The reader used to get the mesh
  vtkPolyDataReader *m_DataReader;

  // Triangle strip filter
  vtkStripper *m_Stripper;

  // Convert to triangles filter
  vtkTriangleFilter *m_Triangulator;

  // Compute normals filter
  vtkPolyDataNormals *m_NormalsFilter;

  // Cleaner for polygon data
  vtkCleanPolyData *m_CleanFilter;

  // Curves that have been traced (?)
  TracerCurves m_Curves;

  // Constant indicating no focus
  const static int NO_FOCUS;

  // Curve under focus
  int m_FocusCurve;

  // Point under focus
  int m_FocusPoint;

  // Currently selected marker (not very relevant)
  vtkIdType m_FocusMarker;

  // A structure representing a region marker
  struct MarkerData {
    string Name;
    Vec Color;
  };

  // List of markers (a marker is a cell in the mesh)
  typedef map<vtkIdType, MarkerData> MarkerMap;
  typedef MarkerMap::iterator MarkerIterator;
  MarkerMap m_Markers;

  // Listener list (m_TracerDataListener)
  pyListenerArrayMacro(TracerDataListener,ITracerDataListener);

  // Broadcast methods
  pyBroadcastEventMacro(TracerDataListener, ITracerDataListener,
    OnFocusCurveDataChange, TracerDataEvent);
  
  pyBroadcastEventMacro(TracerDataListener, ITracerDataListener,
    OnFocusPointChange, TracerDataEvent);

  pyBroadcastEventMacro(TracerDataListener, ITracerDataListener,
    OnFocusCurveChange, TracerDataEvent);

  pyBroadcastEventMacro(TracerDataListener, ITracerDataListener,
    OnCurveListChange, TracerDataEvent);

  pyBroadcastEventMacro(TracerDataListener, ITracerDataListener,
    OnMeshChange, TracerDataEvent);

  pyBroadcastEventMacro(TracerDataListener, ITracerDataListener,
    OnEdgeWeightsUpdate, TracerDataEvent);

  pyBroadcastEventMacro(TracerDataListener, ITracerDataListener,
    OnMarkerListChange, TracerDataEvent);

  pyBroadcastEventMacro(TracerDataListener, ITracerDataListener,
    OnFocusMarkerChange, TracerDataEvent);

  // Set the focus curve, firing the associated event
  void SetFocusCurve(int inFocusCurve);

  // Set the focus point, firing the associated event
  void SetFocusPoint(int inFocusPoint);

  // Compute distances to a source point
  void ComputeDistances(int inFocusPoint);
  
  /** Set the edge weight function to another mode */
  void UpdateEdgeWeightFunction(MeshEdgeWeightFunction *fnNew);
};

#endif
