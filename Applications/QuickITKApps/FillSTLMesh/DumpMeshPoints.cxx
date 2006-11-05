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
  cout << "dumpmeshpoints - Dumps points in a VTK mesh, q-hull compatible" << endl;
  cout << "usage: " << endl;
  cout << "   dumpmeshpoints [options] <input.byu|input.vtk> " << endl;
  cout << "options: " << endl;
  cout << "   -a NAME          Dump the contents of array name as a column" << endl;
  return -1;
}

int main(int argc, char **argv)
{
  // Check the parameters
  if(argc < 2) return usage();

  // Read the appropriate mesh
  vtkPolyData *p1 = ReadVTKData(argv[argc-1]);

  // Get the list of arrays that need to be dumped along with points
  std::vector<vtkDataArray *> vdat;
  for(size_t i = 1; i < argc - 1; ++i)
    {
    if(0 == strcmp(argv[i], "-a"))
      {
      vtkDataArray *arr = p1->GetPointData()->GetArray(argv[++i]);
      if(arr)
        vdat.push_back(arr);
      else
        cerr << "No array by name " << argv[i] << " found in the mesh" << endl;
      }
    }

  // Get the points from the mesh
  vtkPoints *pts = p1->GetPoints();

  // Scale all the points in the mesh
  cout << pts->GetNumberOfPoints() << endl;
  cout << 3 + vdat.size() << endl;
  for(size_t iPoint = 0; iPoint < pts->GetNumberOfPoints(); iPoint++)
    {
    cout << pts->GetPoint(iPoint)[0] << " " 
      << pts->GetPoint(iPoint)[1] << " " 
      << pts->GetPoint(iPoint)[2];
    for(size_t k = 0; k < vdat.size(); k++)
      cout << " " << vdat[k]->GetTuple1(iPoint);
    cout << endl;
    }
}
