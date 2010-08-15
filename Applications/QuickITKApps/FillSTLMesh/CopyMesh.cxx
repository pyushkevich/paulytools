#include "ReadWriteVTK.h"
#include <iostream>

using namespace std;

int usage()
{
  cout << "usage: copymesh input.vtk output.vtk" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  if(argc < 4) return usage();

  vtkPolyData *p = ReadVTKData(argv[1]);
  WriteVTKData(p, argv[2]);
}
