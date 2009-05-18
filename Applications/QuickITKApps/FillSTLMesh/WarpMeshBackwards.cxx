#include <iostream>
#include <string>
#include <sstream>
#include <set>

#include <vnl/vnl_inverse.h>

#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPoints.h>
#include <itkImageFileReader.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageToVectorImageFilter.h>
#include "ReadWriteVTK.h"

using namespace std;
using namespace itk;

int usage()
{
  cout << "warpmesh - Applies a warp field (Brian's format) to a VTK mesh" << endl;
  cout << "usage: " << endl;
  cout << "   warpmesh [options] mesh.vtk out.vtk warp_images" << endl;
  cout << "options: " << endl;
  cout << "   -w spec    : Warp coordinate specification. The warp field gives " << endl;
  cout << "                 a displacement in some coordinate space. It can be " << endl;
  cout << "                 voxel space (ijk), scaled+offset voxel space (ijkos) " << endl;
  cout << "                 lps or ras. Also, 'ants' for ANTS warps " << endl;
  cout << "   -m spec    : Mesh coordinate specification" << endl;
  return -1;
}

enum Coord { IJK, IJKOS, LPS, RAS, ANTS };

bool scan_coord(const char *text, Coord &val)
{
  if(!strcmp(text, "ijk")) val = IJK;
  else if(!strcmp(text, "ijkos")) val = IJKOS;
  else if(!strcmp(text, "lps")) val = LPS;
  else if(!strcmp(text, "ras")) val = RAS;
  else if(!strcmp(text, "ants")) val = ANTS;
  else return false;
  return true;
}

/** 
 * This static function constructs a NIFTI matrix from the ITK direction
 * cosines matrix and Spacing and Origin vectors
 */
vnl_matrix_fixed<double,4,4> ConstructNiftiSform(
  vnl_matrix<double> m_dir, 
  vnl_vector<double> v_origin,
  vnl_vector<double> v_spacing)
{
  // Set the NIFTI/RAS transform
  vnl_matrix<double> m_ras_matrix;
  vnl_diag_matrix<double> m_scale, m_lps_to_ras;
  vnl_vector<double> v_ras_offset;

  // Compute the matrix
  m_scale.set(v_spacing);
  m_lps_to_ras.set(vnl_vector<double>(3, 1.0));
  m_lps_to_ras[0] = -1;
  m_lps_to_ras[1] = -1;
  m_ras_matrix = m_lps_to_ras * m_dir * m_scale;

  // Compute the vector
  v_ras_offset = m_lps_to_ras * v_origin;

  // Create the larger matrix
  vnl_vector<double> vcol(4, 1.0);
  vcol.update(v_ras_offset);

  vnl_matrix_fixed<double,4,4> m_sform;
  m_sform.set_identity();
  m_sform.update(m_ras_matrix);
  m_sform.set_column(3, vcol);
  return m_sform;
}

vnl_matrix_fixed<double,4,4> ConstructVTKtoNiftiTransform(
  vnl_matrix<double> m_dir, 
  vnl_vector<double> v_origin,
  vnl_vector<double> v_spacing)
{
  vnl_matrix_fixed<double,4,4> vox2nii = ConstructNiftiSform(m_dir, v_origin, v_spacing);
  vnl_matrix_fixed<double,4,4> vtk2vox; 
  vtk2vox.set_identity();
  for(size_t i = 0; i < 3; i++)
    {
    vtk2vox(i,i) = 1.0 / v_spacing[i];
    vtk2vox(i,3) = - v_origin[i] / v_spacing[i];
    }
  return vox2nii * vtk2vox;
}

struct WarpMeshParam
{
  string fnMeshIn;
  string fnMeshOut;
  Coord mesh_coord, warp_coord;
  string fnWarp[3];
  WarpMeshParam()
    {
    fnMeshIn = "";
    fnMeshOut = "";
    for(size_t d = 0; d < 3; d++)
      fnWarp[d] = "";
    mesh_coord = RAS;
    warp_coord = RAS;
    }
};

template <class TMeshType> 
TMeshType * ReadMesh(const char *fname)
{ return NULL; }

template <>
vtkUnstructuredGrid *ReadMesh<>(const char *fname)
{
  vtkUnstructuredGridReader *reader = vtkUnstructuredGridReader::New();
  reader->SetFileName(fname);
  reader->Update();
  return reader->GetOutput();
}

template <>
vtkPolyData *ReadMesh<>(const char *fname)
{
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(fname);
  reader->Update();
  return reader->GetOutput();
}


template <class TMeshType> 
void WriteMesh(TMeshType *mesh, const char *fname)
{ }

template <>
void WriteMesh<>(vtkUnstructuredGrid *mesh, const char *fname)
{
  vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
  writer->SetFileName(fname);
  writer->SetInput(mesh);
  writer->Update();
}

template <>
void WriteMesh<>(vtkPolyData *mesh, const char *fname)
{
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(fname);
  writer->SetInput(mesh);
  writer->Update();
}

/**
 * The actual method is templated over the VTK data type (unstructured/polydata)
 */
template <class TMeshType>
int WarpMesh(WarpMeshParam &parm)
{
  // Read warp field
  typedef itk::Image<double, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::LinearInterpolateImageFunction<ImageType> FuncType;

  // Read each of the images in turn
  ImageType::Pointer warp[3];
  ReaderType::Pointer reader[3];
  FuncType::Pointer func[3];
  
  for(size_t d = 0; d < 3; d++)
    {
    reader[d] = ReaderType::New();
    reader[d]->SetFileName(parm.fnWarp[d].c_str());
    reader[d]->Update();
    warp[d] = reader[d]->GetOutput();
    func[d] = FuncType::New();
    func[d]->SetInputImage(warp[d]);
    }

  // Read the mesh
  TMeshType *mesh = ReadMesh<TMeshType>(parm.fnMeshIn.c_str());

  // Set up the transforms
  vnl_matrix_fixed<double, 4, 4> ijk2ras = ConstructNiftiSform(
    warp[0]->GetDirection().GetVnlMatrix(),
    warp[0]->GetOrigin().GetVnlVector(),
    warp[0]->GetSpacing().GetVnlVector());

  vnl_matrix_fixed<double, 4, 4> vtk2ras = ConstructVTKtoNiftiTransform(
    warp[0]->GetDirection().GetVnlMatrix(),
    warp[0]->GetOrigin().GetVnlVector(),
    warp[0]->GetSpacing().GetVnlVector());

  vnl_matrix_fixed<double, 4, 4> lps2ras;
  lps2ras.set_identity();
  lps2ras(0,0) = -1; lps2ras(1,1) = -1;

  vnl_matrix_fixed<double, 4, 4> ras2ijk = vnl_inverse(ijk2ras);
  vnl_matrix_fixed<double, 4, 4> ras2vtk = vnl_inverse(vtk2ras);
  vnl_matrix_fixed<double, 4, 4> ras2lps = vnl_inverse(lps2ras);

  cout << "RAS transform " << endl;
  cout << ijk2ras << endl;

  // Update the coordinates
  for(int k = 0; k < mesh->GetNumberOfPoints(); k++)
    {
    // Get the point (in whatever format that it's stored)
    vnl_vector_fixed<double, 4> x_mesh, x_ras, x_ijk, v_warp, v_ras;
    x_mesh[0] = mesh->GetPoint(k)[0]; x_mesh[1] = mesh->GetPoint(k)[1]; x_mesh[2] = mesh->GetPoint(k)[2];
    x_mesh[3] = 1.0;

    // Map the point into RAS coordinates
    if(parm.mesh_coord == RAS)
      x_ras = x_mesh;
    else if(parm.mesh_coord == LPS)
      x_ras = lps2ras * x_mesh;
    else if(parm.mesh_coord == IJKOS)
      x_ras = vtk2ras * x_mesh;
    else 
      x_ras = ijk2ras * x_mesh;

    // Map the point to IJK coordinates (continuous index)
    x_ijk = ras2ijk * x_ras;
    FuncType::ContinuousIndexType idx;
    idx[0] = x_ijk[0]; idx[1] = x_ijk[1]; idx[2] = x_ijk[2];

    // Interpolate the warp at the point
    // cout << "Evaluate at index " << idx[0] << " " << idx[1] << " " << idx[2] << endl;
    for(size_t d = 0; d < 3; d++)
      v_warp[d] = func[d]->EvaluateAtContinuousIndex(idx);
    v_warp[3] = 0.0;

    // Compute the displacement in RAS coordinates
    if(parm.warp_coord == RAS)
      v_ras = v_warp;
    else if(parm.warp_coord == LPS)
      v_ras = lps2ras * v_warp;
    else if(parm.warp_coord == IJKOS)
      v_ras = vtk2ras * v_warp;
    else if(parm.warp_coord == ANTS)
      {
      // vector is multiplied by the direction matrix ??? Really???
      v_ras = lps2ras * v_warp;
      }
    else 
      v_ras = ijk2ras * v_warp;

    // Add displacement
    vnl_vector_fixed<double, 4> y_ras, y_mesh;
    y_ras = x_ras + v_ras;

    // Map new coordinate to desired system
    if(parm.mesh_coord == RAS)
      y_mesh = y_ras;
    else if(parm.mesh_coord == LPS)
      y_mesh = ras2lps * y_ras;
    else if(parm.mesh_coord == IJKOS)
      y_mesh = ras2vtk * y_ras;
    else 
      y_mesh = ras2ijk * y_ras;

    // Map the transformation to RAS
    // vnl_vector_fixed<double, 4> w_ras = 
    //  wimg[0]->GetSpacingOriginPhysicalSpaceToRASPhysicalSpaceMatrix() * w;

    // vnl_vector_fixed<double, 4> w_ras = w;
      
    /*
     * OLD CODE THT WORKED WITH ANTS
    // Assume transformation is in spacing units
    itk::FixedArray<double, 3> aw, awras;
    aw[0] = w[0];
    aw[1] = w[1];
    aw[2] = w[2];
    wimg[0]->TransformLocalVectorToPhysicalVector(aw, awras);
    p->GetPoints()->SetPoint(k, pt[0] + awras[0], pt[1] + awras[1], pt[2] + awras[2]);
    */

    mesh->GetPoints()->SetPoint(k, y_mesh[0], y_mesh[1], y_mesh[2]);
    }
    
  // Write the mesh
  WriteMesh<TMeshType>(mesh, parm.fnMeshOut.c_str());
  return 0;
}


  

int main(int argc, char **argv)
{
  // Check the parameters
  if(argc < 4) return usage();

  // Parse the optional parameters
  int ch;
  WarpMeshParam parm;
  while ((ch = getopt(argc, argv, "m:w:")) != -1)
    {
    switch(ch)
      {
    case 'm': if(!scan_coord(optarg, parm.mesh_coord) || parm.mesh_coord==ANTS)
                return usage(); 
              break;
    case 'w': if(!scan_coord(optarg, parm.warp_coord))
                return usage(); 
              break;
    default: return usage();
      }
    }

  // Parse the filenames
  if(optind + 5 != argc) return usage();
  parm.fnMeshIn = argv[optind++];
  parm.fnMeshOut = argv[optind++];
  for(size_t d = 0; d < 3; d++)
    parm.fnWarp[d] = argv[optind++];

  // Check the data type of the input file
  vtkDataReader *reader = vtkDataReader::New();
  reader->SetFileName(parm.fnMeshIn.c_str());
  reader->OpenVTKFile();
  reader->ReadHeader();

  // Is this a polydata?
  if(reader->IsFileUnstructuredGrid())
    {
    reader->Delete();
    return WarpMesh<vtkUnstructuredGrid>(parm);
    }
  else if(reader->IsFilePolyData())
    {
    reader->Delete();
    return WarpMesh<vtkPolyData>(parm);
    }
  else
    {
    reader->Delete();
    cerr << "Unsupported VTK data type in input file" << endl;
    return -1;
    }
}
