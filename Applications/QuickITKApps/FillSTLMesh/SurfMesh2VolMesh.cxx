#include <iostream>
#include <string>
#include <sstream>
#include <set>

#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <itkImageFileReader.h>
#include <itkLinearInterpolateImageFunction.h>
//#include <itkOrientedRASImage.h>
#include <itkImage.h>
#include "ReadWriteVTK.h"
#include "tetgen.h"

using namespace std;
using namespace itk;

int usage()
{
  cout << "SurfMesh2VolMesh - Converts a VTK surface mesh to a VTK volume mesh" << endl;
  cout << "usage: " << endl;
  cout << "   SurfMesh2VolMesh [options] meshin.vtk meshout.vtk " << endl;
  cout << "options: " << endl;
  return -1;
}

int main(int argc, char **argv)
{
  
  tetgenio in, out;
  tetgenio::facet *f;
  tetgenio::polygon *p;


  // Check the parameters
  if(argc < 3) return usage();


  // Read the mesh
  vtkPolyData *mesh = ReadVTKData(argv[1]);

  // Cast mesh into tetgen format
  // All indices start from 1.
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
  in.save_nodes("vtkin");
  in.save_poly("vtkin");
  in.save_elements("vtkin");
  in.save_faces("vtkin");

  tetrahedralize("pq1.414a0.1", &in, &out);
  
/*
  // Add the points to the poly data
  vtkPoints *outpoints = vtkPoints::New();
  outpoints->Allocate(out.numberofpoints);
    
  // Allocate the polydata
  vtkPolyData *outmesh = vtkPolyData::New();
  outmesh->Allocate(xModel->GetNumberOfTriangles());
  outmesh->SetPoints(outpoints);
*/

  // Output mesh to files 'barout.node', 'barout.ele' and 'barout.face'.
  out.save_nodes("vtkout");
  in.save_poly("vtkout");
  out.save_elements("vtkout");
  out.save_faces("vtkout");



/*
  // Update the coordinates
  for(size_t k = 0; k < p->GetNumberOfPoints(); k++)
    {
    // Get the point (in RAS coords)
    double *pt = p->GetPoint(k);
    ImageType::PointType itk_point;
    itk_point[0] = pt[0]; itk_point[1] = pt[1]; itk_point[2] = pt[2];
    
    // Map the point to a continuous index
    FuncType::ContinuousIndexType idx;
    wimg[0]->TransformRASPhysicalPointToContinuousIndex(itk_point, idx); 
    cout << "Evaluate at index " << idx[0] << " " << idx[1] << " " << idx[2] << endl;

    // Compute the transformation. We assume the transformation is in ITK
    // physical space. 
    vnl_vector_fixed<double, 4> w;
    w[0] = warp[0]->EvaluateAtContinuousIndex(idx);
    w[1] = warp[1]->EvaluateAtContinuousIndex(idx);
    w[2] = warp[2]->EvaluateAtContinuousIndex(idx);
    w[3] = 0.0;

    // Map the transformation to RAS
    // vnl_vector_fixed<double, 4> w_ras = 
    //  wimg[0]->GetSpacingOriginPhysicalSpaceToRASPhysicalSpaceMatrix() * w;

    // vnl_vector_fixed<double, 4> w_ras = w;
      
    // Assume transformation is in spacing units
    itk::FixedArray<double, 3> aw, awras;
    aw[0] = w[0];
    aw[1] = w[1];
    aw[2] = w[2];
    wimg[0]->TransformLocalVectorToPhysicalVector(aw, awras);

    p->GetPoints()->SetPoint(k, pt[0] + awras[0], pt[1] + awras[1], pt[2] + awras[2]);
    // p->GetPoints()->SetPoint(k, pt[0] + w_ras[0], pt[1] + w_ras[1], pt[2] + w_ras[2]);
    }
    
  // Write the mesh
  WriteVTKData(p, argv[3]);
*/
  return 0;
}
