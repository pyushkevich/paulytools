#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkPolyDataReader.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char *argv[])
  {
  // Usage
  cout << "USAGE: dumpmeshpoints input.vtk output.txt" << endl;

  // Read the mesh
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(argv[1]);
  reader->Update();

  vtkPolyData *pd = reader->GetOutput();
  vtkPoints *points = pd->GetPoints();

  // Create output file
  ofstream fout(argv[2],ios_base::out);
  fout << points->GetNumberOfPoints() << endl << endl;

  // Dump all the points
  for(unsigned int i=0;i<points->GetNumberOfPoints();i++)
    {
    fout << points->GetPoint(i)[0] << " ";
    fout << points->GetPoint(i)[1] << " ";
    fout << points->GetPoint(i)[2] << endl;
    }

  // Close f
  fout.close();
  
  return 1;
  }

