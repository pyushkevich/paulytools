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

int usage()
{
  cout << "dumpmeshpoints - Dumps points in a VTK mesh, q-hull compatible" << endl;
  cout << "usage: " << endl;
  cout << "   dumpmeshpoints <input.byu|input.vtk> " << endl;
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
  vtkTriangleFilter *tri = vtkTriangleFilter::New();
  tri->SetInput(p1);
  cout << "Converting to triangles ..." << endl;
  tri->Update();
  vtkPolyData *pd = tri->GetOutput();

  return pd;
}

vtkPolyData *WriteVTKData(string fn, vtkPolyData *data)
{
  // Choose the writer based on extension
  if(fn.rfind(".byu") == fn.length() - 4)
    {
    vtkBYUWriter *writer = vtkBYUWriter::New();
    writer->SetGeometryFileName(fn.c_str());
    writer->SetInput(data);
    writer->Update();
    }
  else if(fn.rfind(".vtk") == fn.length() - 4)
    {
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetFileName(fn.c_str());
    writer->SetInput(data);
    writer->Update();
    }
  else
    {
    cout << "Could not find a writer for " << fn << endl;
    }
}


void ProgressCommand(Object *, const EventObject &, void *)
{ 
  cout << "." << flush; 
}

int main(int argc, char **argv)
{
  // Check the parameters
  if(argc != 2) return usage();

  // Get the file names
  string fn = argv[argc-1];

  // Read the appropriate meshes
  vtkPolyData *p1 = ReadVTKData(fn);
  vtkPoints *pts = p1->GetPoints();

  // Scale all the points in the mesh
  cout << pts->GetNumberOfPoints() << endl;
  cout << 3 << endl;
  for(size_t iPoint = 0; iPoint < pts->GetNumberOfPoints(); iPoint++)
    cout << pts->GetPoint(iPoint)[0] << " " 
      << pts->GetPoint(iPoint)[1] << " " 
      << pts->GetPoint(iPoint)[2] << endl;
}
