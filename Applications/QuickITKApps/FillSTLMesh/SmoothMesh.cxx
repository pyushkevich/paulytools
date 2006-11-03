#include "ReadWriteVTK.h"
#include <vtkWindowedSincPolyDataFilter.h>
#include <iostream>

using namespace std;

int usage()
{
  cout << "usage: smoothmesh input.vtk output.vtk pass_band" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  if(argc < 4) return usage();

  vtkPolyData *p = ReadVTKData(argv[1]);
  vtkWindowedSincPolyDataFilter *fltSmooth = vtkWindowedSincPolyDataFilter::New();
  fltSmooth->SetInput(p);
  fltSmooth->SetNumberOfIterations(100);
  fltSmooth->SetPassBand(atof(argv[3]));
  fltSmooth->Update();
  WriteVTKData(fltSmooth->GetOutput(), argv[2]);
}
