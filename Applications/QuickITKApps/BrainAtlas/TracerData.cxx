#include "TracerData.h"
#include <fstream>

TracerData
::TracerData()
{
  // Initialize the mesh edge weight function
  m_EdgeWeightFunction = new EuclideanDistanceMeshEdgeWeightFunction();

  // Initialize the shortest distance computer
  m_DistanceMapper = new VTKMeshShortestDistance();
  
  // Create the necessary filters
  m_DataReader = vtkPolyDataReader::New();
  m_Stripper = vtkStripper::New();
  m_Triangulator = vtkTriangleFilter::New();
  m_NormalsFilter = vtkPolyDataNormals::New();
  m_CleanFilter = vtkCleanPolyData::New();
  
  m_DisplayMesh = m_Mesh = NULL;
  m_FocusPoint = -1;
  m_FocusCurve = -1;
}

TracerData
::~TracerData()
{
  delete m_EdgeWeightFunction;
  delete m_DistanceMapper;
  
  m_DataReader->Delete();
  m_Stripper->Delete();
  m_Triangulator->Delete();
  m_NormalsFilter->Delete();
  m_CleanFilter->Delete();
  
}

void
TracerData
::LoadInputMesh(const char *file)
{
  // Read the data from disk
  m_DataReader->SetFileName(file);
  m_DataReader->Update();
  m_DataReader->GetOutput();

  // Clean the input data
  m_CleanFilter->SetInput(m_DataReader->GetOutput());
  m_CleanFilter->SetTolerance(0);
  m_CleanFilter->Update();

  // Convert the input to triangles
  m_Triangulator->PassLinesOff();
  m_Triangulator->PassVertsOff();
  m_Triangulator->SetInput(m_CleanFilter->GetOutput());
  m_Triangulator->Update();

  // Compute all normals 
  m_NormalsFilter->SetInput(m_Triangulator->GetOutput());
  m_NormalsFilter->ConsistencyOn();
  m_NormalsFilter->AutoOrientNormalsOn();
  m_NormalsFilter->NonManifoldTraversalOn();
  m_NormalsFilter->Update();
  
  // Get the output, store it for subsequent use
  m_Mesh = m_NormalsFilter->GetOutput();

  // Set up the stripper
  m_Stripper->SetInput(m_NormalsFilter->GetOutput());
  m_Stripper->Update();
  m_DisplayMesh = m_Stripper->GetOutput();

  // Send the data to the distance mapper
  m_DistanceMapper->SetInputMesh(m_Mesh);
  m_DistanceMapper->ComputeGraph();
}

void
TracerData
::SaveCurves(const char *file)
{
  ofstream fout(file,ios_base::out);
  m_Curves.SaveAsText(fout);
  fout.close();
}

bool
TracerData
::LoadCurves(const char *file)
{
  ifstream fin(file,ios_base::in);
  bool rc = m_Curves.LoadAsText(fin);
  fin.close();

  if(rc)
    {
    m_FocusCurve = -1;
    m_FocusPoint = -1;
    }

  return rc;
}
