#include "TracerData.h"

TracerData
::TracerData()
{
  // Create the reader
  m_DataReader = vtkPolyDataReader::New();
  m_DistanceMapper = NULL;
  m_Mesh = NULL;
  m_Stripper = vtkStripper::New();

  m_FocusPoint = -1;
  m_FocusCurve = -1;
}

TracerData
::~TracerData()
{
  if(m_DistanceMapper)
    delete m_DistanceMapper;
  m_DataReader->Delete();
  m_Stripper->Delete();
}

void
TracerData
::LoadInputMesh(const char *file)
{
  // Read the data from disk
  m_DataReader->SetFileName(file);
  m_DataReader->Update();

  // Set up the stripper
  m_Stripper->SetInput(m_DataReader->GetOutput());
  m_Stripper->Update();

  // Triangle strip the data
  m_Mesh = m_Stripper->GetOutput();

  // Send the data to the distance mapper
  if(m_DistanceMapper)
    delete m_DistanceMapper;
  m_DistanceMapper = new VTKMeshShortestDistance(m_Mesh);
}
