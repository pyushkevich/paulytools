#include "MedialModelIO.h"
#include "vtkPolyData.h"
#include "vtkOBJReader.h"
#include "vtkBYUReader.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"

#include <itksys/SystemTools.hxx>

using namespace std;

GenericMedialModel *
MedialModelIO
::ReadModel(const char *file)
{
  // Change to the directory of the registry file
  std::string dirWork = itksys::SystemTools::GetCurrentWorkingDirectory();
  std::string dirFile = itksys::SystemTools::GetFilenamePath(file);

  // Open the file as a registry
  Registry R(file);

  // Read the type of the model
  std::string type = R["Grid.Type"][""];

  // If the grid type is Cartesian, create a new model
  if(type == "Cartesian") 
    {
    CartesianMedialModel *model = CartesianMedialModelIO::ReadModel(R);
    return model;
    }
  else if(type == "LoopSubdivision")
    {
    itksys::SystemTools::ChangeDirectory(dirFile.c_str());
    SubdivisionMedialModel *model = SubdivisionMedialModelIO::ReadModel(R);
    itksys::SystemTools::ChangeDirectory(dirWork.c_str());
    return model;
    }
  else throw ModelIOException("Unknown or missing Grid.Type in the model file");
}

void
MedialModelIO
::WriteModel(GenericMedialModel *model, const char *file)
{
  // Try casting to a cartesian model
  CartesianMedialModel *cart = dynamic_cast<CartesianMedialModel *>(model);
  if(cart != NULL)
    {
    CartesianMedialModelIO::WriteModel(cart, file);
    return;
    }

  // Try casting to a subdivision model
  SubdivisionMedialModel *subd = dynamic_cast<SubdivisionMedialModel *>(model);
  if(subd != NULL)
    {
    SubdivisionMedialModelIO::WriteModel(subd, file);
    return;
    }

  // Something is wrong!
  throw ModelIOException("Unknown model type in MedialModelIO::WriteModel");
}

CartesianMedialModel *
CartesianMedialModelIO::ReadModel(Registry &R)
{
  // Target pointer
  CartesianMedialModel *cmm;

  // Read the grid size specification
  size_t m = R["Grid.Size.U"][0];
  size_t n = R["Grid.Size.V"][0];

  // Make sure the grid specification exists
  if(m == 0 || n == 0)
    throw ModelIOException("Grid specification is missing in Cartesian model");

  // Check if the grid spacing is available
  if(R.Folder("Grid.Spacing.U").GetArraySize() == m && 
    R.Folder("Grid.Spacing.V").GetArraySize() == n)
    {
    // Read the grid spacing
    CartesianMedialModel::Vec uGrid(m, 0.0), vGrid(n, 0.0);
    R.Folder("Grid.Spacing.U").GetArray(uGrid.data_block(), 0.0);
    R.Folder("Grid.Spacing.V").GetArray(vGrid.data_block(), 0.0);

    // Use the grid-based initialization
    cmm = new CartesianMedialModel(uGrid, vGrid);
    }
  else
    {
    cmm = new CartesianMedialModel(m, n);
    }

  // Now, load the surface specification
  if(R["SurfaceModel"][""] != Registry::StringType("Fourier"))
    throw ModelIOException(
      "SurfaceModel is not Fourier when reading CartesianMedialModel");

  // Read the dimensions of the surface model
  size_t mc = R["Fourier.Size.U"][0];
  size_t nc = R["Fourier.Size.V"][0];
  if(mc == 0 || nc == 0 )
    throw ModelIOException(
      "Number of components missing in Fourier surface model specification");

  // Allocate the Fourier surface model
  FourierSurface *surface = new FourierSurface(mc, nc);

  // Stick the surface into the medial model
  cmm->AdoptMedialSurface(surface);

  // Pass the registry to the model to do the rest
  cmm->ReadFromRegistry(R);
  cmm->ComputeAtoms();

  // Return the model
  return cmm;
}

void CartesianMedialModelIO
::WriteModel(CartesianMedialModel *model, const char *file)
{
  // This code really needs to be unified!
  Registry R;
  model->WriteToRegistry(R);
  R.WriteToFile(file);
}

vtkPolyData *
SubdivisionMedialModelIO
::ReadMesh(const string &file, const string &type)
{
  try 
    {
    vtkPolyData *poly = vtkPolyData::New();

    if(type == "OBJ") 
      {
      vtkOBJReader *reader = vtkOBJReader::New();
      reader->SetOutput(poly);
      reader->SetFileName(file.c_str());
      reader->Update();
      reader->Delete();
      }
    else if (type == "BYU") 
      {
      vtkBYUReader *reader = vtkBYUReader::New();
      reader->SetOutput(poly);
      reader->SetFileName(file.c_str());
      reader->Update();
      reader->Delete();
      }
    else if (type == "VTK") 
      {
      vtkPolyDataReader *reader = vtkPolyDataReader::New();
      reader->SetOutput(poly);
      reader->SetFileName(file.c_str());
      reader->Update();
      reader->Delete();
      }
    else
      throw ModelIOException("Unknown mesh type in Subdivision Model");

    return poly;
    }
  catch(...)
    {
    throw ModelIOException("Error loading VTK mesh in Subdivision Model");
    }
} 

SubdivisionMedialModel *
SubdivisionMedialModelIO
::ReadModel(Registry &R)
{
  size_t i;

  // Read the information about the coefficient-level mesh
  string fnMesh = R["Grid.Model.Coefficient.FileName"][""];
  string sMeshType = R["Grid.Model.Coefficient.FileType"][""];

  // Read the polydata containing the coefficient map
  vtkPolyData *poly = ReadMesh(fnMesh, sMeshType);

  // Read the number of subdivisions to apply to the data
  size_t nSubs = R["Grid.Model.Atom.SubdivisionLevel"][0];
  if(nSubs == 0)
    throw ModelIOException("Missing SubdivisionLevel in model file");

  // Generate a mesh level from the model
  SubdivisionSurface::MeshLevel mesh;
  SubdivisionSurface::ImportLevelFromVTK(poly, mesh);

  // Now we must read the coefficient information. It is normally contained
  // in the poly itself, although rho might not be available
  SMLVec3d *xCtl = new SMLVec3d[mesh.nVertices];
  for(i = 0; i < mesh.nVertices; i++)
    xCtl[i].copy_in(poly->GetPoint(i));

  // Get the default rho (constant)
  double xDefaultRho = R["Grid.Model.Coefficient.ConstantRho"][-0.25];

  // Generate the rho array. We initialize to a constant, and then use the rho
  // array in the polydata. This may not exist.
  vnl_vector<double> rho(mesh.nVertices, xDefaultRho);
  if(poly->GetPointData()->GetArray("Rho") != NULL) 
    {
    vtkDataArray *daRho = poly->GetPointData()->GetScalars("Rho");
    for(i = 0; i < mesh.nVertices; i++)
      rho[i] = daRho->GetTuple1(i);
    }

  // Delete the poly data
  poly->Delete();

  // Create the medial model
  SubdivisionMedialModel *smm = new SubdivisionMedialModel();
  smm->SetMesh(mesh, xCtl, rho.data_block(), nSubs, 0);
  smm->ComputeAtoms();

  return smm;
}

void
SubdivisionMedialModelIO
::WriteModel(SubdivisionMedialModel *model, const char *file)
{
  // First of all, we have to determine a filename where to save the mesh. It
  // should be in the same directory and name as the file, but with a different
  // extension. So basically, we need to strip the extension on the file
  std::string fnreg(file);
  std::string fnbase = itksys::SystemTools::GetFilenameWithoutLastExtension(fnreg);
  std::string fnmesh = fnbase + ".vtk";

  // Create a registry to save the data
  Registry R;

  // Save the model file information
  R["Grid.Type"] << "LoopSubdivision";
  R["Grid.Model.Coefficient.FileName"] << fnmesh;
  R["Grid.Model.Coefficient.FileType"] << "VTK";

  // Save the subdivision level info
  R["Grid.Model.Atom.SubdivisionLevel"] << model->GetSubdivisionLevel();

  // Save the registry
  R.WriteToFile(file);

  // Save the mesh level as a VTK mesh
  vtkPolyData *poly = vtkPolyData::New();

  // Allocate the array for the vertices
  vtkPoints *points = vtkPoints::New();
  points->Allocate(model->mlCoefficient.nVertices);
  poly->SetPoints(points);

  // Allocate the array for rho
  vtkDoubleArray *xRho = vtkDoubleArray::New();
  xRho->SetName("Rho");
  xRho->SetNumberOfComponents(1);
  xRho->SetNumberOfTuples(model->mlCoefficient.nVertices);
  poly->GetPointData()->AddArray(xRho);

  // Set the values of X and Rho
  for(size_t i = 0; i < model->mlCoefficient.nVertices; i++)
    {
    // Set the point's coordinates
    points->InsertNextPoint(model->xCoefficients.data_block() + i * 4);

    // Set rho
    xRho->SetTuple1(i, model->xCoefficients[i * 4 + 3]);
    }

  // Populate the cells
  SubdivisionSurface::ExportLevelToVTK(model->mlCoefficient, poly);

  // Save the model
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(poly);
  writer->SetFileName(fnmesh.c_str());
  writer->SetFileTypeToBinary();
  writer->Update();

  // Clean up
  writer->Delete();
  xRho->Delete();
  points->Delete();
  poly->Delete();
}
