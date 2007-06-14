#include "ReadWriteVTK.h"
#include "vtkPoints.h"
#include "vtkQuadricClustering.h"
#include "vtkCleanPolyData.h"
#include "vtkTriangleFilter.h"
#include <iostream>

using namespace std;

int usage()
{
  cout << "quadcluster - apply clustering to a mesh" << endl;
  cout << "usage: " << endl;
  cout << "  quadcluster input.vtk output.vtk X Y Z" << endl;
  cout << "parameters: " << endl;
  cout << "  X Y Z       Size of the clustering element in mm" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  size_t i, j;
  if(argc < 5) return usage();

  // Get the cluster size
  double cx = atof(argv[3]);
  double cy = atof(argv[4]);
  double cz = atof(argv[5]);

  // Read the input
  vtkPolyData *input = ReadVTKData(argv[1]);

  // Triangulate ?
  vtkTriangleFilter *fltTriangle = vtkTriangleFilter::New();
  fltTriangle->SetInput(input);
  fltTriangle->Update();
  vtkPolyData *tri = fltTriangle->GetOutput();

  // Get the extents of the mesh
  double bounds[6];
  tri->ComputeBounds();
  tri->GetBounds(bounds);

  // Compute the number of divisions
  int nx = (int) ceil((bounds[1] - bounds[0]) / cx);
  int ny = (int) ceil((bounds[3] - bounds[2]) / cy);
  int nz = (int) ceil((bounds[5] - bounds[4]) / cz);
  printf("Using %dx%dx%d divisions\n", nx, ny, nz);

  // Create the clustering filter
  vtkQuadricClustering *flt = vtkQuadricClustering::New();
  flt->SetInput(tri);
  flt->SetNumberOfDivisions(nx,ny,nz);
  flt->SetNumberOfDivisions(50,50,50);
  // flt->SetNumberOfDivisions((int) cx, (int) cy, (int) cz);
  flt->Update();

  vtkCleanPolyData *cln = vtkCleanPolyData::New();
  cln->SetInput(flt->GetOutput());
  cln->Update();

  WriteVTKData(cln->GetOutput(), argv[2]);
}

