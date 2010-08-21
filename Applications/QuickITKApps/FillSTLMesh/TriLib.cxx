#include "TriLib.h"

#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkTriangle.h>
#include <vtkCellData.h>

#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_det.h>


using namespace std;

typedef vnl_matrix_fixed<double, 4, 4> MatrixType;



void WriteTriangleOutputAsPolyDataMesh(triangulateio &out, const char *fname)
{
  vtkPoints* outpoints = vtkPoints::New();
  vtkFloatArray* pointarray = vtkFloatArray::New();
  vtkFloatArray* cellarray  = vtkFloatArray::New();
  vtkPolyData* outmesh = vtkPolyData::New();

  // Add the points 
  outpoints->SetNumberOfPoints(out.numberofpoints);
  pointarray->SetName ("Point Area");
  pointarray->Allocate(out.numberofpoints);
  cout <<  out.numberofpoints << " points in triangle mesh" << endl;
  cout <<  out.numberoftriangleattributes << " cell arrays in triangle mesh" << endl;
  cout <<  out.numberofpointattributes << " point arrays in triangle mesh" << endl;
  for(int k = 0; k < out.numberofpoints; k++)
    {
    double pt[3];
    pt[0] = out.pointlist[ k*2 ];
    pt[1] = out.pointlist[ k*2 + 1 ];
    pt[2] = 1.0;
    outpoints->SetPoint(k, pt[0], pt[1], pt[2]);
    pointarray->InsertNextValue(0.0);
    }
  outmesh->SetPoints(outpoints);
  

  // Allocate the mesh
  outmesh->Allocate(out.numberoftriangles);
  
  cellarray->SetName ("Cell Area");
  cellarray->Allocate(out.numberoftriangles);
  // Add the cells
  for(int k = 0; k < out.numberoftriangles; k++)
    {
    vtkIdList* idlist = vtkIdList::New();
    unsigned int ids[4];
    ids[0] = out.trianglelist[ k*3 ];
    ids[1] = out.trianglelist[ k*3 + 1 ];
    ids[2] = out.trianglelist[ k*3 + 2 ];
    
    idlist->InsertNextId (ids[0]);
    idlist->InsertNextId (ids[1]);
    idlist->InsertNextId (ids[2]);
    outmesh->InsertNextCell (VTK_TRIANGLE, idlist);

    // Calculate the area
    double p[3][3];
    outmesh->GetPoint(ids[0], p[0]);
    outmesh->GetPoint(ids[1], p[1]);
    outmesh->GetPoint(ids[2], p[2]);
    double area = vtkTriangle::TriangleArea(p[0], p[1], p[2]);
    cellarray->InsertNextValue(area);

    // Distribute the area
    for(size_t ivert = 0; ivert < 3; ivert++)
       pointarray->SetTuple1(ids[ivert], pointarray->GetTuple1(ids[ivert]) + area/3.0);
    idlist->Delete();
    }
  if (outmesh->GetCellData())
    {
    cout << "Adding cell data.." << endl;
    outmesh->GetCellData()->AddArray(cellarray);
    }
    
  if (outmesh->GetPointData())
    {
    cout << "Adding point data.." << endl;
    outmesh->GetPointData()->AddArray(pointarray);
    }

  cout <<  outmesh->GetNumberOfPoints() << " points in output mesh" << endl;
  // Output mesh to files 'barout.node', 'barout.ele' and 'barout.face'.
  /*
  out.save_nodes("vtkout");
  out.save_poly("vtkout");
  out.save_elements("vtkout");
  out.save_faces("vtkout");
  */

  // Write the vtk output
//  WriteVTKData( outmesh, fname);
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetFileName(fname);
    writer->SetInput(outmesh);
    writer->Update();

}
