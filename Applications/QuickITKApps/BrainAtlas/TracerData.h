#ifndef __TracerData_h_
#define __TracerData_h_

#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkStripper.h"
#include "vtkPolyDataNormals.h"
#include "VTKMeshShortestDistance.h"
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

  // Get the coordinate for a given point index (in edge-mesh)
  void GetPointCoordinate(vtkIdType id, Vec &target)
    { m_DistanceMapper->GetInputMesh()->GetPoint(id,target.data_block()); }

  // Create a new curve (given a name)
  void AddNewCurve(const char *name)
    {
    SetFocusPoint(-1);
    SetFocusCurve(m_Curves.AddCurve(name));
    }
    
  // Set the current curve
  void SetCurrentCurve(IdType iCurve)
    {
    // Check that the curve exists
    assert(m_Curves.IsCurvePresent(iCurve));

    // Set the focus point
    SetFocusPoint( m_Curves.GetCurveControls(iCurve).size()
      ? (int) m_Curves.GetCurveControls(iCurve).back() : -1 );

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
    assert(m_FocusPoint != -1 && m_FocusCurve != -1);
    assert(m_FocusPoint == 
      m_Curves.GetCurveControls(m_FocusCurve).back());

    // Delete the last point in the curve
    m_Curves.DeleteLastControlPointInCurve(m_FocusCurve);

    // Set the focus point if the curve isn't empty
    SetFocusPoint( m_Curves.GetCurveControls(m_FocusCurve).size()
      ? (int) m_Curves.GetCurveControls(m_FocusCurve).back() : -1 );

    // Nofity listeners that the current curve's data has changed
    TracerDataEvent evt(this);
    BroadcastOnFocusCurveDataChange(&evt);
    }
    
  void DeleteCurrentCurve() 
    {
    // We need to have a current curve
    assert(m_FocusCurve != -1);

    // Delete the current curve
    m_Curves.DeleteCurve(m_FocusCurve);

    // Don't focus on any other curve
    SetFocusPoint(-1);
    SetFocusCurve(-1);

    // Notify listeners that the list of curves has changed
    TracerDataEvent evt(this);
    BroadcastOnCurveListChange(&evt);
    }

  /** Is a curve currently selected, or is there no current curve? */
  bool IsCurrentCurveValid()
    { return m_FocusCurve != -1; }
   
  IdType GetCurrentCurve()
    { return m_FocusCurve; }

  unsigned int GetNumberOfCurvePoints(IdType iCurve)
    {
    assert(m_Curves.IsCurvePresent(iCurve));
    return m_Curves.GetCurveVertices(iCurve).size();
    }

  // Given a ray, find a point closest to that ray
  bool PickPoint(Vec xStart, Vec xRay, vtkIdType &iPoint)
    { return m_DistanceMapper->PickPoint(xStart,xRay,iPoint); }

  // Create a path from the last point on the focused curve (if any)
  // to the specified point. The path includes the specified point, but
  // not the last point on the focused curve
  bool GetPathToPoint(vtkIdType iPoint, TracerCurves::MeshCurve &target)
    {
    // Clear the return array
    target.clear();

    // Do we have a 'focus' point
    if(m_FocusPoint != -1 && m_DistanceMapper->IsVertexConnected(iPoint))
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
          return -1;
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
    { return m_FocusPoint != -1; }

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
    if(m_FocusPoint != -1)
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

  // Curve under focus
  int m_FocusCurve;

  // Point under focus
  int m_FocusPoint;

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
