#ifndef __TracerData_h_
#define __TracerData_h_

#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkStripper.h"
#include "VTKMeshShortestDistance.h"
#include <list>

#include "TracerCurves.h"

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

  // Load the mesh data from disk
  void LoadInputMesh(const char *file);
  
  // Save the curves to disk
  void SaveCurves(const char *file);

  // Load curves from disk
  bool LoadCurves(const char *file);

  // Get the 'source' mesh
  vtkPolyData *GetSourceMesh() 
    {
    return m_Mesh;
    }

  // Get the coordinate for a given point index (in edge-mesh)
  void GetPointCoordinate(vtkIdType id, Vec &target)
    {
    m_DistanceMapper->GetEdgeMesh()->GetPoint(id,target.data_block());
    }

  // Create a new curve (given a name)
  void AddNewCurve(const char *name)
    {
    m_FocusCurve = m_Curves.AddCurve(name);
    m_FocusPoint = -1;
    }
    
  // Set the current curve
  void SetCurrentCurve(IdType iCurve)
    {
    // Check that the curve exists
    assert(m_Curves.IsCurvePresent(iCurve));
    m_FocusCurve = iCurve;

    // Set the focus point
    if(m_Curves.GetCurveControls(iCurve).size())
      {
      // Get the point of focus
      m_FocusPoint = m_Curves.GetCurveControls(iCurve).back();

      // Compute shortest distances to that point
      m_DistanceMapper->ComputeDistances(
        m_Curves.GetControlPointVertex(m_FocusPoint));
      }
    else
      {
      m_FocusPoint = -1;
      }
    }

  const TracerCurves *GetCurves()
    { return &m_Curves; }

  void SetCurveName(IdType iCurve, const char *name)
    {
    assert(m_Curves.IsCurvePresent(iCurve));
    m_Curves.SetCurveName(iCurve, name);
    }
    
  void DeleteCurrentCurve() 
    {
/*
    // Erase the current curve
    m_Curves.erase(m_Curves.begin() + m_FocusCurve);
    
    // Update the curve index
    m_FocusCurve = (m_Curves.size() > 0) ? 
      (m_FocusCurve + 1) % m_Curves.size() : -1;
*/
    }
   
  IdType GetCurrentCurve()
    {
    return m_FocusCurve;
    }

  unsigned int GetNumberOfCurvePoints(IdType iCurve)
    {
    assert(m_Curves.IsCurvePresent(iCurve));
    return m_Curves.GetCurveVertices(iCurve).size();
    }

  /* 
  Vec GetCurvePoint(unsigned int iCurve, unsigned int iPoint)
    {
    assert(m_Curves.IsCurvePresent(iCurve));

    const TracerCurves::MeshCurve &points = m_Curves.GetCurveVertices(iCurve);
    
    assert(iPoint < points.size());

    return GetPointCoordinate(points[iPoint]);
    }
  */

  // Given a ray, find a point closest to that ray
  bool PickPoint(Vec xStart, Vec xRay, vtkIdType &iPoint)
    {
    return m_DistanceMapper->PickPoint(xStart,xRay,iPoint);
    }

  // Create a path from the last point on the focused curve (if any)
  // to the specified point. The path includes the specified point, but
  // not the last point on the focused curve
  void GetPathToPoint(vtkIdType iPoint, TracerCurves::MeshCurve &target)
    {
    // Clear the return array
    target.clear();

    // Do we have a 'focus' point
    if(m_FocusPoint != -1)
      {
      vtkIdType iLast = m_Curves.GetControlPointVertex(m_FocusPoint);
      vtkIdType iCurrent = iPoint;
      while(iCurrent != iLast)
        {
        target.push_front(iCurrent);
        iCurrent = m_DistanceMapper->GetVertexPredecessor(iCurrent);
        }
      }
    else
      {
      target.push_front(iPoint);
      }
    }

  // Check if there is a path source
  bool IsPathSourceSet()
    {
    return m_FocusPoint != -1;
    }

  // Get the vertex to which the path currently extends
  vtkIdType GetPathSource()
    {
    return m_Curves.GetControlPointVertex(m_FocusPoint);
    }

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

    // Compute distances from this point
    cout << "computing distances to point " << iPoint << endl;
    m_DistanceMapper->ComputeDistances(iPoint);
    cout << "done" << endl;

    // Save the added point as the focus point
    m_FocusPoint = iNewControl;
    }

private:
  // Shortest distance computer
  VTKMeshShortestDistance *m_DistanceMapper;

  // The source mesh
  vtkPolyData *m_Mesh;

  // The reader used to get the mesh
  vtkPolyDataReader *m_DataReader;

  // Triangle strip filter
  vtkStripper *m_Stripper;

  // Curves that have been traced (?)
  TracerCurves m_Curves;

  // Curve under focus
  int m_FocusCurve;

  // Point under focus
  int m_FocusPoint;
};

#endif
