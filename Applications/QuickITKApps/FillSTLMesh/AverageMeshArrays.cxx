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
  cout << "avgmesharr - average a list of VTK mesh arrays" << endl;
  cout << "usage: " << endl;
  cout << "  avgmesharr mesh1.vtk ... meshN.vtk ArrName target.vtk output.vtk" << endl;
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
  unsigned int nMeshes = argc - 4;
  vtkPolyData **xMeshes = new vtkPolyData *[nMeshes];

  // Read the meshes
  for(i = 0; i < nMeshes; i++)
    {
    // Read the polydata
    xMeshes[i] = ReadVTKData(argv[1 + i]);

    cout << "." << flush;    
    // reader->Delete();
    }
  cout << endl;

  // HACK: add up arrays
  vtkFloatArray *dist = vtkFloatArray::New();
  dist->SetName(argv[argc-3]);
  dist->SetNumberOfComponents(1);
  dist->SetNumberOfTuples(xMeshes[0]->GetNumberOfPoints());
  for(j = 0; j < xMeshes[0]->GetNumberOfPoints(); j++)
    {
    double m = 0.0;
    for(i = 0; i < nMeshes; i++)
      {
      vtkDataArray *arr = xMeshes[i]->GetPointData()->GetArray(argv[argc-3]);
      m += arr->GetTuple1(j);
      }
    m /= nMeshes;
    dist->SetTuple1(j, m);
    }
  vtkPolyData *target = ReadVTKData(argv[argc-2]);
  target->GetPointData()->AddArray(dist);

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

  vtkPolyDataWriter *pdw = vtkPolyDataWriter::New();
  pdw->SetFileName(argv[argc - 1]);
  pdw->SetInput(target);
  pdw->Update();
}
