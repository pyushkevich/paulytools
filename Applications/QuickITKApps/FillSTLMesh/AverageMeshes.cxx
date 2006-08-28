#include "vtkPolyData.h"
#include "vtkBYUWriter.h"
#include "vtkBYUReader.h"
#include <vtkSTLReader.h>
#include "vtkPoints.h"
#include "vtkProcrustesAlignmentFilter.h"
#include "vtkLandmarkTransform.h"
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <iostream>

using namespace std;

int usage()
{
  cout << "avgmesh - average a list of BYU meshes" << endl;
  cout << "usage: " << endl;
  cout << "  avgmesh mesh1.buy ... meshN.byu output.byu" << endl;
  return -1;
}

vtkPolyData *ReadVTKData(string fn)
{
  vtkPolyData *p1 = NULL;

  // Choose the reader based on extension

  if(fn.rfind(".byu") == fn.length() - 4)
    {
    vtkBYUReader *reader = vtkBYUReader::New();
    reader->SetFileName(fn.c_str());
    reader->Update();
    p1 = reader->GetOutput();
    }
  else if(fn.rfind(".vtk") == fn.length() - 4)
    {
    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    reader->SetFileName(fn.c_str());
    reader->Update();
    p1 = reader->GetOutput();
    }
  else
    {
    cout << "Could not find a reader for " << fn << endl;
    cout << fn.rfind(".byu") << endl;
    cout << fn.rfind(".vtk") << endl;
    return NULL;
    }

  // Convert the model to triangles
  return p1;
}
int main(int argc, char *argv[])
{
  size_t i, j;
  if(argc < 3) return usage();

  // Load a list of meshes
  unsigned int nMeshes = argc - 2;
  vtkPolyData **xMeshes = new vtkPolyData *[nMeshes];

  // Create a procrustes alignment filter
  vtkProcrustesAlignmentFilter *fltAlign = 
    vtkProcrustesAlignmentFilter::New();
  fltAlign->SetNumberOfInputs(nMeshes);
  fltAlign->GetLandmarkTransform()->SetModeToRigidBody();

  // Read the meshes
  for(i = 0; i < nMeshes; i++)
    {
    // Read the polydata
    xMeshes[i] = ReadVTKData(argv[1 + i]);

    // Plug into procrustes
    fltAlign->SetInput(i, xMeshes[i]);
    
    cout << "." << flush;    
    // reader->Delete();
    }
  cout << endl;

  // Compute the alignment
  fltAlign->Update();
  xMeshes[0]->SetPoints(fltAlign->GetMeanPoints());

  // Add all meshes to the first mesh
  /*
  for(unsigned int j = 0; j < xMeshes[0]->GetNumberOfPoints(); j++)
    {
    double *xTarget = xMeshes[0]->GetPoint(j);
    for(unsigned int d = 0; d < 3; d++)
      {
      for(unsigned int k = 1; k < nMeshes; k++)
        {
        xTarget[d] += xMeshes[k]->GetPoint(j)[d];
        }
      xTarget[d] /= nMeshes;
      }
    }
    */

  // Save the average
  if(strstr(argv[argc - 1], ".byu") != NULL)
    {
    vtkBYUWriter *writer = vtkBYUWriter::New();
    writer->SetGeometryFileName(argv[argc - 1]);
    writer->SetInput(xMeshes[0]);
    writer->Update();
    }
  else
    {
    vtkPolyDataWriter *pdw = vtkPolyDataWriter::New();
    pdw->SetFileName(argv[argc - 1]);
    pdw->SetInput(xMeshes[0]);
    pdw->Update();
    }
}
