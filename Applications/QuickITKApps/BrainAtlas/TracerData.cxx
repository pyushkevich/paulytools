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

void
TracerData
::SaveCurves(const char *file)
{
  FILE *fout = fopen(file,"wt");
  
  fprintf(fout,"# Mesh surface curve annotation file\n");
  fprintf(fout,"# Commands: \n");
  fprintf(fout,"# NUMCURVES number_of_curves\n");
  fprintf(fout,"# CURVE number_of_points \"name\" \n");
  fprintf(fout,"# POINT vertex_id x y z\n");
  
  fprintf(fout,"NUMCURVES %d\n",m_Curves.size()); 
  for(unsigned int i=0;i<m_Curves.size();i++)
    {
    fprintf(fout,"CURVE %d \"%s\"\n",m_Curves[i].points.size(),m_Curves[i].name.c_str());
    for(unsigned int j=0;j<m_Curves[i].points.size();j++)
      {
      fprintf(fout,"POINT %d %lg %lg %lg\n",
        m_Curves[i].points[j].i,m_Curves[i].points[j].x[0],
        m_Curves[i].points[j].x[1],m_Curves[i].points[j].x[2]);
      }
    }
  fclose(fout);
}