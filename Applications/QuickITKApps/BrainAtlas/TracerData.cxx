#include "TracerData.h"

TracerData
::TracerData()
{
  // Create the reader
  m_DataReader = vtkPolyDataReader::New();
  m_DistanceMapper = NULL;
  m_Mesh = NULL;
}

TracerData
::~TracerData()
{
  if(m_DistanceMapper)
    delete m_DistanceMapper;
  m_DataReader->Delete();
}

void
TracerData
::LoadInputMesh(const char *file)
{
  // Read the data from disk
  m_DataReader->SetFileName(file);
  m_DataReader->Update();
  m_Mesh = m_DataReader->GetOutput();
  m_Mesh->BuildCells();
  m_Mesh->BuildLinks();

  // Send the data to the distance mapper
  if(m_DistanceMapper)
    delete m_DistanceMapper;
  m_DistanceMapper = new VTKMeshShortestDistance(m_Mesh);
}
