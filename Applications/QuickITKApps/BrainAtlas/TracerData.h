#ifndef __TracerData_h_
#define __TracerData_h_

#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkStripper.h"
#include "VTKMeshShortestDistance.h"



class TracerData
{
public:
  // Typedef for the 3d vector
  typedef vnl_vector_fixed<double, 3> Vec;

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
    m_FocusCurve = m_Curves.size();
    }

  // Given a ray, find a point closest to that ray
  vtkIdType PickPoint(Vec xStart, Vec xRay)
    {
    return m_DistanceMapper->PickPoint(xStart,xRay);
    }

  // Add a point to the current curve
  void AddNewPoint(vtkIdType iPoint)
    {
    assert(m_FocusCurve >= 0);
    
    Vertex p;
    p.i = iPoint;
    p.x = GetPointCoordinate(iPoint);

    m_Curves[m_FocusCurve].points.push_back(p);
 
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
