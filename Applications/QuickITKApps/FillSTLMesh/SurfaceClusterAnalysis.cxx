#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "ReadWriteVTK.h"
#include "vtkTriangleFilter.h"
#include "vtkBandedPolyDataContourFilter.h"

int usage()
{
  cout << "This program performs cluster analysis on a VTK mesh" << endl;
  cout << "Usage: " << endl;
  cout << "    meshcluster mesh.vtk ..." << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  if(argc < 2)
    return usage();

  // Read the input mesh and convert to triangles
  vtkPolyData *mesh = ReadVTKData(argv[1]);
  vtkTriangleFilter *fTri = vtkTriangleFilter::New();
  fTri->SetInput(mesh);
  fTri->Update();

  // Set the active scalar array
  fTri->GetOutput()->GetPointData()->SetActiveScalars(argv[3]);

  // Use the banded contour filter
  vtkBandedPolyDataContourFilter *fContour = 
    vtkBandedPolyDataContourFilter::New();
  fContour->SetInput(fTri->GetOutput());
  fContour->SetValue(0, 0.01);
  fContour->Update();

  // Save the output of this thingy
  WriteVTKData(fContour->GetOutput(), argv[2]);

}
