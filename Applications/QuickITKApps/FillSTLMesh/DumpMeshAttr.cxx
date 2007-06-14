#include <iostream>
#include <string>
#include <sstream>
#include <set>

#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkLODActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkBYUReader.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkBYUWriter.h>
#include <vtkSTLWriter.h>
#include <vtkTriangleFilter.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>

#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkCommand.h>
#include <itkImage.h>

#include "DrawTriangles.h"
#include "ReadWriteVTK.h"

using namespace std;
using namespace itk;

int usage()
{
  cout << "dumpmeshattr - Dumps an attribute array in a VTK mesh" << endl;
  cout << "usage: " << endl;
  cout << "   dumpmeshattr input.vtk array_name " << endl;
  return -1;
}

int main(int argc, char **argv)
{
  // Check the parameters
  if(argc < 3) return usage();
  
  // Read the appropriate mesh
  vtkPolyData *p = ReadVTKData(argv[1]);
  vtkDataArray *arr = p->GetPointData()->GetArray(argv[2]);

  // Write the array values
  for(size_t iPoint = 0; iPoint < arr->GetNumberOfTuples(); iPoint++)
    {
    for(size_t iComp = 0; iComp < arr->GetNumberOfComponents(); iComp++)
      {
      cout << arr->GetComponent(iPoint, iComp);
      if(iComp + 1 < arr->GetNumberOfComponents())
        cout << " ";
      else 
        cout << endl;
      }
    }
}
