#include <iostream>
#include <string>
#include <sstream>

#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
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
  cout << "flipmesh - flip mesh across the plane a*x + b*y + c*z + d = 0" << endl;
  cout << "usage: " << endl;
  cout << "   scalemesh input.vtk output.vtk a b c w" << endl;
  return -1;
}

int main(int argc, char **argv)
{
  // Check the parameters
  if(argc != 7) return usage();

  // Get the file names
  string fn1 = argv[1];
  string fn2 = argv[2];
  double ax = atof(argv[3]);
  double ay = atof(argv[4]);
  double az = atof(argv[5]);
  double b = atof(argv[6]);

  // Read the input mesh
  vtkPolyData *p1 = ReadVTKData(fn1);
  
  for(size_t iPoint = 0; iPoint < p1->GetNumberOfPoints(); iPoint++)
    {
    // Reflect the point
    double x[3]; p1->GetPoint(iPoint, x);
    double q = - (ax * x[0] + ay * x[1] + az * x[2] + b) 
      / (ax * ax + ay * ay + az * az);
    x[0] += 2 * q * ax;
    x[1] += 2 * q * ay;
    x[2] += 2 * q * az;
    p1->GetPoints()->SetPoint(iPoint, x);

    // Compute the reflection for all vector data
    for(size_t i = 0; i < p1->GetPointData()->GetNumberOfArrays(); i++)
      {
      vtkDataArray *da = p1->GetPointData()->GetArray(i);
      if(da->GetNumberOfComponents() == 3)
        {
        double v[3];
        da->GetTuple(iPoint, v);

        double r = - (ax * v[0] + ay * v[1] + az * v[2]) 
          / (ax * ax + ay * ay + az * az);
        v[0] += 2 * r * ax;
        v[1] += 2 * r * ay;
        v[2] += 2 * r * az;

        // If we are working with the normals, change the sign for display purposes
        if(da == p1->GetPointData()->GetNormals())
          { v[0] *= -1; v[1] *= -1; v[2] *= -1; }

        da->SetTuple(iPoint, v);
        }
      }
    }

  WriteVTKData(p1, fn2);
}
