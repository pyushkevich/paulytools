#include "ReadWriteVTK.h"
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <iostream>

using namespace std;

int usage()
{
  cout << "mesharrstat - print out statistics about an array in a mesh" << endl;
  cout << "usage: " << endl;
  cout << "  mesharrstat mesh.vtk arrname" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  if(argc < 3) return usage();


  vtkPolyData *xMesh = ReadVTKData(argv[1]);

  

  // HACK: add up arrays
  double m = 0.0;
  vtkDataArray *arr = xMesh->GetPointData()->GetArray(argv[2]);
  for(int j = 0; j < xMesh->GetNumberOfPoints(); j++)
    {
    m += arr->GetTuple1(j);
    }
    m /= xMesh->GetNumberOfPoints();
  
  std::cout << m << std::endl;

}
