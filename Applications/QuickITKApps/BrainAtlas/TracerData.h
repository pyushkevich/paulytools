#ifndef __TracerData_h_
#define __TracerData_h_

#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkStripper.h"
#include "VTKMeshShortestDistance.h"
#include <list>


class TracerData
{
public:
  // Typedef for the 3d vector
  typedef vnl_vector_fixed<double, 3> Vec;
  typedef vector<Vec> VecArray;

  // Structure representing a point on the curve
  struct Vertex {
    vtkIdType i;
    Vec x;
  };

  // Structure representing a curve on the surface of the grey matter
  struct Curve {
    string name;
    vector<Vertex> points;
  };

  // Tracer data constructors
  TracerData();
  ~TracerData();

  // Load the mesh data from disk
  void LoadInputMesh(const char *file);

  // Get the 'source' mesh
  vtkPolyData *GetSourceMesh() 
    {
    return m_Mesh;
    }

  // Get the coordinate for a given point index (in edge-mesh)
  Vec GetPointCoordinate(vtkIdType id)
    {
    Vec x;
    m_DistanceMapper->GetEdgeMesh()->GetPoint(id,x.data_block());
    return x;
    }

  // Create a new curve (given a name)
  void AddNewCurve(const char *name)
    {
    Curve cnew;
    cnew.name = name;
    m_Curves.push_back(cnew);
    m_FocusCurve = m_Curves.size() - 1;
    }

  unsigned int GetNumberOfCurves() 
    {
    return m_Curves.size();
    }

  const char *GetCurveName(unsigned int iCurve)
    {
    return m_Curves[iCurve].name.c_str();
    }

  int GetCurrentCurve()
    {
    return m_FocusCurve;
    }

  unsigned int GetNumberOfCurvePoints(unsigned int iCurve)
    {
    assert(iCurve < m_Curves.size());
    return m_Curves[iCurve].points.size();
    }

  Vec GetCurvePoint(unsigned int iCurve, unsigned int iPoint)
    {
    assert(iCurve < m_Curves.size());
    assert(iPoint < m_Curves[iCurve].points.size());
    return m_Curves[iCurve].points[iPoint].x;
    }

  // Given a ray, find a point closest to that ray
  bool PickPoint(Vec xStart, Vec xRay, vtkIdType &iPoint)
    {
    return m_DistanceMapper->PickPoint(xStart,xRay,iPoint);
    }

  // Create a path from the last point on the focused curve (if any)
  // to the specified point. The path includes the specified point, but
  // not the last point on the focused curve
  list<vtkIdType> GetPathToPoint(vtkIdType iPoint)
    {
    list<vtkIdType> lPoints;

    if(m_Curves[m_FocusCurve].points.size())
      {
      cout << "tracing back points" << endl;
      cout << "distance to point is " << m_DistanceMapper->GetVertexDistance(iPoint) << endl;

      vtkIdType iLast = m_Curves[m_FocusCurve].points.back().i;
      vtkIdType iCurrent = iPoint;
      while(iCurrent != iLast)
        {
        lPoints.push_front(iCurrent);
        iCurrent = m_DistanceMapper->GetVertexPredecessor(iCurrent);
        cout << iCurrent << endl;
        }
      }
    else
      {
      lPoints.push_front(iPoint);
      }

    return lPoints;
    }

  // Add a point to the current curve
  void AddNewPoint(vtkIdType iPoint)
    {
    assert(m_FocusCurve >= 0);
    
    // Get the path to the current point
    list<vtkIdType> lPoints = GetPathToPoint(iPoint);

    // Add the points on the list
    list<vtkIdType>::const_iterator it;
    for(it = lPoints.begin();it != lPoints.end(); it++)
      {
      Vertex p;
      p.i = *it; p.x = GetPointCoordinate(p.i);
      m_Curves[m_FocusCurve].points.push_back(p);
      cout << "adding point " << p.i << endl;
      }

    // Compute distances from this point
    cout << "computing distances to point " << iPoint << endl;
    m_DistanceMapper->ComputeDistances(iPoint);
    cout << "done" << endl;
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
  vector<Curve> m_Curves;

  // Curve under focus
  int m_FocusCurve;

  // Point under focus
  int m_FocusPoint;
};

#endif
