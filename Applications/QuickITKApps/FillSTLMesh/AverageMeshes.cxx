#include "vtkPolyData.h"
#include "vtkBYUWriter.h"
#include "vtkBYUReader.h"
#include "vtkPoints.h"
#include <iostream>

using namespace std;

int usage()
{
  cout << "avgmesh - average a list of BYU meshes" << endl;
  cout << "usage: " << endl;
  cout << "  avgmesh mesh1.buy ... meshN.byu output.byu" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  if(argc < 3) return usage();

  // Load a list of meshes
  unsigned int nMeshes = argc - 2;
  vtkPolyData **xMeshes = new vtkPolyData *[nMeshes];
  for(unsigned int i = 0; i < nMeshes; i++)
    {
    vtkBYUReader *reader = vtkBYUReader::New();
    reader->SetFileName(argv[1 + i]);
    reader->Update();

    xMeshes[i] = reader->GetOutput();
    cout << "." << flush;
    // reader->Delete();
    }
  cout << endl;

  // Add all meshes to the first mesh
  for(unsigned int j = 0; j < xMeshes[0]->GetNumberOfPoints(); j++)
    {
    float *xTarget = xMeshes[0]->GetPoint(j);
    for(unsigned int d = 0; d < 3; d++)
      {
      for(unsigned int k = 1; k < nMeshes; k++)
        {
        xTarget[d] += xMeshes[k]->GetPoint(j)[d];
        }
      xTarget[d] /= nMeshes;
      }
    }

  // Save the average
  vtkBYUWriter *writer = vtkBYUWriter::New();
  writer->SetGeometryFileName(argv[argc - 1]);
  writer->SetInput(xMeshes[0]);
  writer->Update();
}
