#include "ReadWriteVTK.h"
#include "vtkPoints.h"
#include "vtkButterflySubdivisionFilter.h"
#include <iostream>

using namespace std;

int usage()
{
  cout << "submesh - subdivide a mesh" << endl;
  cout << "usage: " << endl;
  cout << "  submesh mesh.buy output.byu levels" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  size_t i, j;
  if(argc < 4) return usage();

  vtkPolyData *input = ReadVTKData(argv[1]);
  
  vtkButterflySubdivisionFilter *fltButterfly = vtkButterflySubdivisionFilter::New();
  fltButterfly->SetInput(input);
  fltButterfly->SetNumberOfSubdivisions(atoi(argv[3]));
  fltButterfly->Update();

  WriteVTKData(fltButterfly->GetOutput(), argv[2]);
}

