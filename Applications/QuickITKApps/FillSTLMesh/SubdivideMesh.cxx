#include "vtkPolyData.h"
#include "vtkBYUWriter.h"
#include "vtkBYUReader.h"
#include <vtkSTLReader.h>
#include "vtkPoints.h"
#include "vtkButterflySubdivisionFilter.h"
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPointData.h>
#include <iostream>

using namespace std;

int usage()
{
  cout << "submesh - subdivide a mesh" << endl;
  cout << "usage: " << endl;
  cout << "  submesh mesh.buy output.byu levels" << endl;
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
  if(argc < 4) return usage();

  vtkPolyData *input = ReadVTKData(argv[1]);
  
  vtkButterflySubdivisionFilter *fltButterfly = vtkButterflySubdivisionFilter::New();
  fltButterfly->SetInput(input);
  fltButterfly->SetNumberOfSubdivisions(atoi(argv[3]));
  fltButterfly->Update();

  if(strstr(argv[argc - 1], ".byu") != NULL)
    {
    vtkBYUWriter *writer = vtkBYUWriter::New();
    writer->SetGeometryFileName(argv[2]);
    writer->SetInput(fltButterfly->GetOutput());
    writer->Update();
    }
  else
    {
    vtkPolyDataWriter *pdw = vtkPolyDataWriter::New();
    pdw->SetFileName(argv[2]);
    pdw->SetInput(fltButterfly->GetOutput());
    pdw->Update();
    }
}
