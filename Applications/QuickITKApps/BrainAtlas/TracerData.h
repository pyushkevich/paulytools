#ifndef __TracerData_h_
#define __TracerData_h_

#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "VTKMeshShortestDistance.h"

class TracerData
{
public:
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
 
private:
  // Shortest distance computer
  VTKMeshShortestDistance *m_DistanceMapper;

  // The source mesh
  vtkPolyData *m_Mesh;

  // The reader used to get the mesh
  vtkPolyDataReader *m_DataReader;

  // Curves that have been traced (?)

  
};

#endif
