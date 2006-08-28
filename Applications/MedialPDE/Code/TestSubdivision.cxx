#include "vnl/vnl_random.h"
#include "MedialAtom.h"
#include "MedialModelIO.h"
#include "vtkOBJReader.h"
#include "vtkBYUWriter.h"
#include "vtkPolyData.h"
#include "SubdivisionSurface.h"
#include "SubdivisionMedialModel.h"
#include "SubdivisionSurfaceMedialIterationContext.h"
#include "ScriptInterface.h"

#include <string>
#include <iostream>

void ExportMedialMeshToVTK(GenericMedialModel *model, const char *file);
void ExportBoundaryMeshToVTK(GenericMedialModel *model, const char *file);

using namespace std;

/**
 * Test subdivision surface functionality
 */
int TestSubdivisionSurface(const char *objMesh)
{
  // Read the input mesh
  vtkOBJReader *reader = vtkOBJReader::New();
  reader->SetFileName(objMesh);
  reader->Update();
  vtkPolyData *poly = reader->GetOutput();

  // Create a subdivision surface
  SubdivisionSurface::MeshLevel mesh;
  SubdivisionSurface::ImportLevelFromVTK(poly, mesh);

  // Describe the mesh
  cout << "Input mesh has " << mesh.triangles.size() << " triangles";
  cout << " and " << mesh.nVertices << " vertices" << endl;

  // Save the input surface for comparison
  vtkBYUWriter *writer = vtkBYUWriter::New();
  writer->SetInput(poly);
  writer->SetGeometryFileName("byu_input.byu");
  writer->Update();

  // Check the mesh after loading
  if(!SubdivisionSurface::CheckMeshLevel(mesh))
    return 1;

  // Subdivide the mesh once
  SubdivisionSurface::MeshLevel meshsub;
  SubdivisionSurface::Subdivide(&mesh, &meshsub);

  cout << "Subdivided mesh has " << meshsub.triangles.size() << " triangles";
  cout << " and " << meshsub.nVertices << " vertices" << endl;

  // Check the subdivided mesh
  if(!SubdivisionSurface::CheckMeshLevel(meshsub))
    return 1;

  // Compute the subdivision surface
  vtkPolyData *polysub = vtkPolyData::New();
  SubdivisionSurface::ApplySubdivision(poly, polysub, meshsub);

  // Save the subdivision surface
  writer->SetInput(polysub);
  writer->SetGeometryFileName("byu_subdivide.byu");
  writer->Update();

  // Subdivide the mesh one more time
  SubdivisionSurface::MeshLevel meshresub;
  SubdivisionSurface::Subdivide(&meshsub, &meshresub);

  cout << "Subdivided mesh has " << meshresub.triangles.size() << " triangles";
  cout << " and " << meshresub.nVertices << " vertices" << endl;

  // Check the subdivided mesh
  if(!SubdivisionSurface::CheckMeshLevel(meshresub))
    return 1;

  // Save the subdivision surface
  SubdivisionSurface::ApplySubdivision(poly, polysub, meshresub);
  writer->SetInput(polysub);
  writer->SetGeometryFileName("byu_resubdivide.byu");
  writer->Update();



  return 0;
}

/**
 * Test sparse matrix multiplication
 */
int TestSparseMatrix()
{
  // Create a sparse matrix with random elements
  size_t i, j, k, N1 = 5, N2 = 4, N3 = 6;
  typedef vnl_sparse_matrix<int> Mutable;
  typedef ImmutableSparseMatrix<int> Immutable;

  // Create two mutable matrices
  Mutable A(N1, N2), B(N2, N3), C(N1, N3);

  // Initialize them with some random values (-10 to 10)
  vnl_random rnd;
  for(i = 0; i < N1; i++) for(j = 0; j < N2; j++)
    if(rnd.lrand32(0, 4) == 0)
      A(i,j) = rnd.lrand32(0,18) - 9;

  for(i = 0; i < N2; i++) for(j = 0; j < N3; j++)
    if(rnd.lrand32(0, 4) == 0)
      B(i,j) = rnd.lrand32(0,18) - 9;
  
  // Compute their product using VNL
  A.mult(B, C);

  // Compute the product using immutable matrices
  Immutable AI, BI, CI, DI;
  AI.SetFromVNL(A); 
  BI.SetFromVNL(B);
  CI.SetFromVNL(C);
  Immutable::Multiply(DI, AI, BI);

  // Print the results
  cout << "A is " << AI << endl;
  cout << "B is " << BI << endl;
  cout << "C is " << CI << endl;
  cout << "D is " << DI << endl;

  // Compare CI and DI
  return CI == DI ? 0 : 1;
}

int TestSubdivisionPDE(const char *objMesh)
{
  // Read the input mesh
  vtkOBJReader *reader = vtkOBJReader::New();
  reader->SetFileName(objMesh);
  reader->Update();
  vtkPolyData *poly = reader->GetOutput();

  // Create a subdivision surface
  SubdivisionSurface::MeshLevel mesh;
  SubdivisionSurface::ImportLevelFromVTK(poly, mesh);

  // Get the polygon data
  SMLVec3d *xCtl = new SMLVec3d[mesh.nVertices];
  double *rhoCtl = new double[mesh.nVertices];
  for(size_t i = 0; i < mesh.nVertices; i++)
    {
    xCtl[i][0] = poly->GetPoint(i)[0];
    xCtl[i][1] = poly->GetPoint(i)[1];
    xCtl[i][2] = poly->GetPoint(i)[2];
    rhoCtl[i] = -0.25;
    }

  // Create the medial model
  SubdivisionMedialModel model;
  
  // Pass the mesh to the medial model with specified number of levels
  model.SetMesh(mesh, xCtl, rhoCtl, 2, 0); 
  model.ComputeAtoms();

  // Temp: test derivatives
  int iTest = model.GetSolver()->TestPartialDerivatives();

  // Save the result
  ExportMedialMeshToVTK(&model, "mesh_medial.vtk");
  ExportBoundaryMeshToVTK(&model, "mesh_boundary.vtk");

  // Save the model itself
  MedialModelIO::WriteModel(&model, "subdivision_model.cmrep");

  return iTest;
}

int TestModelSubdivision(const char *file)
{
  medialpde::SubdivisionMPDE mp1(file), mp2(file);
  mp2.SubdivideMeshes(1,0);

  GenericMedialModel *mm1 = mp1.GetMedialModel();
  GenericMedialModel *mm2 = mp2.GetMedialModel();

  double maxerr = 0.0;
  for(size_t i = 0; i < mm1->GetNumberOfAtoms(); i++)
    {
    double del = (mm1->GetAtomArray()[i].X - mm2->GetAtomArray()[i].X).magnitude();
    if(del > maxerr) maxerr = del;
    }

  return (maxerr > 1.0e-10);
}


int usage()
{
  cout << "testsub: MedialPDE Test Module" << endl;
  cout << "  usage: testpde TEST_ID [parameters] " << endl;
  cout << "  tests: " << endl;
  cout << "    SUBSURF1 XX.obj            Test subdivision with OBJ mesh" << endl;
  cout << "    SUBSURF2 XX.obj            Test subdivision PDE" << endl;
  cout << "    SPARSEMAT                  Test sparse matrix routines" << endl;
  cout << "    MODELSUB                   Test model subdivistion" << endl;
  cout << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  // Different tests that can be executed
  if(argc == 1) return usage();

  // Choose a test depending on the parameters
  if(0 == strcmp(argv[1], "SUBSURF1") && argc > 2)
    return TestSubdivisionSurface(argv[2]);
  else if(0 == strcmp(argv[1], "SUBSURF2") && argc > 2)
    return TestSubdivisionPDE(argv[2]);
  else if(0 == strcmp(argv[1], "SPARSEMAT"))
    return TestSparseMatrix();
  else if(0 == strcmp(argv[1], "MODELSUB"))
    return TestModelSubdivision(argv[2]);
  else 
    return usage();
}

