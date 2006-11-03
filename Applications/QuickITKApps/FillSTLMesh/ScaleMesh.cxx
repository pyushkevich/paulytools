#include <iostream>
#include <string>
#include <sstream>

#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkLODActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataMapper.h>
#include <vtkTriangleFilter.h>
#include "ReadWriteVTK.h"

#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkCommand.h>
#include <itkImage.h>

#include "DrawTriangles.h"

using namespace std;
using namespace itk;

#ifndef vtkFloatingPointType
#define vtkFloatingPointType vtkFloatingPointType
typedef float vtkFloatingPointType;
#endif

int usage()
{
  cout << "scalemesh - scale a BYU mesh in x, y and z" << endl;
  cout << "usage: " << endl;
  cout << "   scalemesh input.byu output.byu scale_x scale_y scale_z" << endl;
  return -1;
}

vtkPolyData *ReadAndTriangulateVTKData(string fn)
{
  vtkPolyData *p1 = ReadVTKData(fn);
  if(p1 != NULL)
    {
    // Convert the model to triangles
    vtkTriangleFilter *tri = vtkTriangleFilter::New();
    tri->SetInput(p1);
    cout << "Converting to triangles ..." << endl;
    tri->Update();
    vtkPolyData *pd = tri->GetOutput();
    return pd;
    }
  else return NULL;
}

void ProgressCommand(Object *, const EventObject &, void *)
{ 
  cout << "." << flush; 
}

int main(int argc, char **argv)
{
  // Check the parameters
  if(argc != 6) return usage();

  // Get the file names
  string fn1 = argv[1];
  string fn2 = argv[2];
  double sx = atof(argv[3]);
  double sy = atof(argv[4]);
  double sz = atof(argv[5]);

  // Read the appropriate meshes
  vtkPolyData *p1 = ReadAndTriangulateVTKData(fn1);
  vtkPoints *pts = p1->GetPoints();

  // Scale all the points in the mesh
  for(size_t iPoint = 0; iPoint < pts->GetNumberOfPoints(); iPoint++)
    {
    vtkFloatingPointType *x = pts->GetPoint(iPoint);
    x[0] *= sx;
    x[1] *= sy;
    x[2] *= sz;
    pts->SetPoint(iPoint, x);
    }

  WriteVTKData(p1, fn2);
}
