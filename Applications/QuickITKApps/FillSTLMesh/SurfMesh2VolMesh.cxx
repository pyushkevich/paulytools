#include <iostream>
#include <string>
#include <sstream>
#include <set>

#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <itkImageFileReader.h>
#include <itkLinearInterpolateImageFunction.h>
//#include <itkOrientedRASImage.h>
#include <itkImage.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_det.h>
#include "ReadWriteVTK.h"
#include "tetgen.h"

using namespace std;
using namespace itk;

typedef vnl_matrix_fixed<double, 4, 4> MatrixType;

int usage()
{
  cout << "SurfMesh2VolMesh - Converts a VTK surface mesh to a VTK volume mesh" << endl;
  cout << "usage: " << endl;
  cout << "   SurfMesh2VolMesh meshin.vtk meshout.vtk [tetgenoptions]" << endl;
  cout << "tetgenoptions: default is pq1.414a0.1" << endl;
  return -1;
}

int main(int argc, char **argv)
{
  
  tetgenio in, out;
  tetgenio::facet *f;
  tetgenio::polygon *p;


  // Check the parameters
  if(argc < 3) return usage();

  // tetgen option 
  std::string optstr = std::string("");
  if (argc > 3)
    {
    for(int i = 3; i < argc; i++)
      optstr = optstr + std::string( argv[i] );
    }
  else
    optstr = optstr + std::string("pq1.414a0.1");


  // Read the mesh
  vtkPolyData *mesh = ReadVTKData(argv[1]);

  // Cast mesh into tetgen format
  // All indices start from 0.
  in.firstnumber = 0;

  // Enter the nodes
  in.numberofpoints = mesh->GetNumberOfPoints();
  in.pointlist = new REAL[in.numberofpoints * 3];
  for(size_t k = 0; k < mesh->GetNumberOfPoints(); k++)
    {
    double *pt = mesh->GetPoint(k);
    in.pointlist[ k*3 ] = pt[0];
    in.pointlist[ k*3 + 1 ] = pt[1];
    in.pointlist[ k*3 + 2 ] = pt[2];
    }
  
  // Enter the cells
  in.numberoffacets = mesh->GetNumberOfCells(); 
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];
  for(size_t k = 0; k < mesh->GetNumberOfCells(); k++)
    {
    vtkCell *cell = mesh->GetCell(k);
    if(cell->GetCellType() != VTK_TRIANGLE)
      throw("Wrong cell type");
    vtkIdType a0 = cell->GetPointId(0);
    vtkIdType a1 = cell->GetPointId(1);
    vtkIdType a2 = cell->GetPointId(2);
    f = &in.facetlist[k];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 3;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = a0;
    p->vertexlist[1] = a1;
    p->vertexlist[2] = a2;


    }

  // Output the PLC to files 'barin.node' and 'barin.poly'.
  /*
  in.save_nodes("vtkin");
  in.save_poly("vtkin");
  in.save_elements("vtkin");
  in.save_faces("vtkin");
  */

  //tetrahedralize("pq1.414a0.1", &in, &out);
  tetrahedralize((char *)optstr.c_str(), &in, &out);
  

  vtkPoints* outpoints = vtkPoints::New();
  vtkFloatArray* pointarray = vtkFloatArray::New();
  vtkFloatArray* cellarray  = vtkFloatArray::New();
  vtkUnstructuredGrid* outmesh = vtkUnstructuredGrid::New();

  // Add the points 
  outpoints->SetNumberOfPoints(out.numberofpoints);
  pointarray->SetName ("Point Volume");
  pointarray->Allocate(out.numberofpoints);
  cout <<  out.numberofpoints << " points in tetgen mesh" << endl;
  cout <<  out.numberoftetrahedronattributes << " cell arrays in tetgen mesh" << endl;
  cout <<  out.numberofpointattributes << " point arrays in tetgen mesh" << endl;
  for(size_t k = 0; k < out.numberofpoints; k++)
    {
    double pt[3];
    pt[0] = out.pointlist[ k*3 ];
    pt[1] = out.pointlist[ k*3 + 1 ];
    pt[2] = out.pointlist[ k*3 + 2 ];
    outpoints->SetPoint(k, pt[0], pt[1], pt[2]);
    pointarray->InsertNextValue(0.0);
    }
  outmesh->SetPoints(outpoints);
  
  // Debug
  /*
  for(size_t k = 0; k < outpoints->GetNumberOfPoints(); k++)
    {
    double *pt = outpoints->GetPoint(k);
    cout << " point " << k << " is " << pt << endl;
    }
  */

  // Allocate the mesh
  outmesh->Allocate(out.numberoftetrahedra);
  
  cellarray->SetName ("Cell Volume");
  cellarray->Allocate(out.numberoftetrahedra);
  // Add the cells
  MatrixType mymat;
  mymat.set_identity();
  for(size_t k = 0; k < out.numberoftetrahedra; k++)
    {
    vtkIdList* idlist = vtkIdList::New();
    unsigned int ids[4];
    ids[0] = out.tetrahedronlist[ k*4 ];
    ids[1] = out.tetrahedronlist[ k*4 + 1 ];
    ids[2] = out.tetrahedronlist[ k*4 + 2 ];
    ids[3] = out.tetrahedronlist[ k*4 + 3 ];
    
    idlist->InsertNextId (ids[0]);
    idlist->InsertNextId (ids[1]);
    idlist->InsertNextId (ids[2]);
    idlist->InsertNextId (ids[3]);
    outmesh->InsertNextCell (VTK_TETRA, idlist);
    // Calculate tet volume
    for(size_t ivert = 0; ivert < 4; ivert++)
       {
       double ithvert[3];
       outmesh->GetPoint( ids[ivert], ithvert );
       for (size_t idir = 0; idir < 3; idir++) 
           mymat(ivert,idir) = ithvert[idir] ;
       mymat(ivert, 3) = 1.0;
       
       }
       
    //cout << "mat " << mymat << endl;
    float volume = vnl_det(mymat)/6.0;   
    cellarray->InsertNextValue(abs(volume));
    // Distribute the volume to all vertices equally
    for(size_t ivert = 0; ivert < 4; ivert++)
       pointarray->SetTuple1(ids[ivert], pointarray->GetTuple1(ids[ivert]) + abs(volume)/4.0);
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
  vtkUnstructuredGridWriter* writer = vtkUnstructuredGridWriter::New();
  writer->SetFileName (argv[2]);
  writer->SetInput (outmesh);
  writer->Write();
  writer->Delete();


  return 0;
}
