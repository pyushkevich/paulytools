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
#include <vtkFloatArray.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <itkImageFileReader.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageToVectorImageFilter.h>
#include "ReadWriteVTK.h"

using namespace std;
using namespace itk;

int usage()
{
  cout << "tetsample - samples an image at points on a VTK mesh" << endl;
  cout << "usage: " << endl;
  cout << "   warpmesh [options] mesh.vtk out.vtk image scalarname" << endl;
  cout << "options: " << endl;
  cout << "   -i spec    :  Image coordinate specification. It can be " << endl;
  cout << "                 voxel space (ijk), scaled+offset voxel space (ijkos) " << endl;
  cout << "                 lps or ras. " << endl;
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

struct TetSampleParam
{
  string fnMeshIn;
  string fnMeshOut;
  Coord mesh_coord, warp_coord;
  string fnImage;
  string scalarName;
  TetSampleParam()
    {
    fnMeshIn = "";
    fnMeshOut = "";
    fnImage = "";
    scalarName = "";
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

template <class X> X mymode(X *data, int size)
{
  register int t, w;
  X md, oldmd;
  int count, oldcount;

  oldmd = 0;
  oldcount = 0;
  for(t=0; t<size; t++) {
    md = data[t];
    count = 1;
    for(w = t+1; w < size; w++) 
      if(md==data[w]) count++;
    if(count > oldcount) {
      oldmd = md;
      oldcount = count;
    }
  }
  return oldmd;
}

/**
 * The actual method is templated over the VTK data type (unstructured/polydata)
 */
template <class TMeshType>
int TetSample(TetSampleParam &parm)
{
  // Read warp field
  typedef itk::Image<double, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType> FuncType;

  // Read each of the images in turn
  ImageType::Pointer sampim;
  ReaderType::Pointer reader;
  FuncType::Pointer func;
  
  reader = ReaderType::New();
  reader->SetFileName(parm.fnImage.c_str());
  reader->Update();
  sampim = reader->GetOutput();
  func = FuncType::New();
  func->SetInputImage(sampim);

  // Read the mesh
  TMeshType *mesh = ReadMesh<TMeshType>(parm.fnMeshIn.c_str());

  // Set up the transforms
  vnl_matrix_fixed<double, 4, 4> ijk2ras = ConstructNiftiSform(
    sampim->GetDirection().GetVnlMatrix(),
    sampim->GetOrigin().GetVnlVector(),
    sampim->GetSpacing().GetVnlVector());

  vnl_matrix_fixed<double, 4, 4> vtk2ras = ConstructVTKtoNiftiTransform(
    sampim->GetDirection().GetVnlMatrix(),
    sampim->GetOrigin().GetVnlVector(),
    sampim->GetSpacing().GetVnlVector());

  vnl_matrix_fixed<double, 4, 4> lps2ras;
  lps2ras.set_identity();
  lps2ras(0,0) = -1; lps2ras(1,1) = -1;

  vnl_matrix_fixed<double, 4, 4> ras2ijk = vnl_inverse(ijk2ras);
  vnl_matrix_fixed<double, 4, 4> ras2vtk = vnl_inverse(vtk2ras);
  vnl_matrix_fixed<double, 4, 4> ras2lps = vnl_inverse(lps2ras);

  cout << "RAS transform " << endl;
  cout << ijk2ras << endl;

  // Create the volume array (cell-wise) for future jacobian computation
  vtkFloatArray *sampledScalarCell = vtkFloatArray::New();
  sampledScalarCell->SetName(parm.scalarName.c_str());
  sampledScalarCell->Allocate(mesh->GetNumberOfCells());
  
  vtkFloatArray *sampledScalarPoint = vtkFloatArray::New();
  sampledScalarPoint->SetName(parm.scalarName.c_str());
  sampledScalarPoint->Allocate(mesh->GetNumberOfPoints());

  // Update the coordinates
  for(int k = 0; k < mesh->GetNumberOfPoints(); k++)
    {
    // Get the point (in whatever format that it's stored)
    vnl_vector_fixed<double, 4> x_mesh, x_ras, x_ijk, v_ras;
    float v_image;
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

    // Interpolate the image at the point
    v_image = func->EvaluateAtContinuousIndex(idx);

    sampledScalarPoint->InsertNextValue( static_cast<float>(v_image));
    // cout << "Evaluate at index " << idx[0] << " " << idx[1] << " " << idx[2] << " " << v_image << " " <<  sampledScalarPoint->GetTuple1(k) << endl;

    }
    
  // Assign cell labels as well TODO do it in a sensible way like mode of the point values
  for(int i = 0; i < mesh->GetNumberOfCells(); i++)
    {
      vtkCell *cell =  mesh->GetCell(i);
      float *arr ;
      arr = (float *)malloc(sizeof(float)*cell->GetNumberOfPoints());
      for (int j = 0; j < cell->GetNumberOfPoints(); j++)
          arr[j] = sampledScalarPoint->GetTuple1(cell->GetPointId(j));
      sampledScalarCell->InsertNextValue(mymode<float>(arr, cell->GetNumberOfPoints()));


    }

  // Add both arrays
  mesh->GetCellData()->AddArray(sampledScalarCell);
  mesh->GetPointData()->AddArray(sampledScalarPoint);
  mesh->GetCellData()->SetActiveScalars(parm.scalarName.c_str());
  mesh->GetPointData()->SetActiveScalars(parm.scalarName.c_str());
 

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
  TetSampleParam parm;
  while ((ch = getopt(argc, argv, "m:i:")) != -1)
    {
    switch(ch)
      {
    case 'm': if(!scan_coord(optarg, parm.mesh_coord) || parm.mesh_coord==ANTS)
                return usage(); 
              break;
    case 'i': if(!scan_coord(optarg, parm.warp_coord))
                return usage(); 
              break;
    default: return usage();
      }
    }

  // Parse the filenames
  if(optind + 4 != argc) return usage();
  parm.fnMeshIn = argv[optind++];
  parm.fnMeshOut = argv[optind++];
  parm.fnImage = argv[optind++];
  parm.scalarName = argv[optind];

  // Check the data type of the input file
  vtkDataReader *reader = vtkDataReader::New();
  reader->SetFileName(parm.fnMeshIn.c_str());
  reader->OpenVTKFile();
  reader->ReadHeader();

  // Is this a polydata?
  if(reader->IsFileUnstructuredGrid())
    {
    reader->Delete();
    return TetSample<vtkUnstructuredGrid>(parm);
    }
  else if(reader->IsFilePolyData())
    {
    reader->Delete();
    return TetSample<vtkPolyData>(parm);
    }
  else
    {
    reader->Delete();
    cerr << "Unsupported VTK data type in input file" << endl;
    return -1;
    }
}
