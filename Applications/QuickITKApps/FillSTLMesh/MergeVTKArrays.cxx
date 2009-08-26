#include <iostream>
#include <string>
#include <vector>
#include <vtksys/RegularExpression.hxx>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>

using namespace std;
using namespace vtksys;

int usage()
{
  cout <<
    "vtkmergearr: merge arrays across many VTK meshes\n"
    "usage:\n"
    "  vtkmergearr ref.vtk out.vtk list_file
    "parameters:\n"
    "  ref.vtk  Mesh that will be used as reference. Must match all inputs\n"
    "  out.vtk  Output filename\n"
    "list_file format:\n"
    "  filename.vtk source_array target_array [source_array target_array] ... \n"
    "array specification format:\n"
    "  
      
  return -1;
}


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

struct MergeParms
{
  const char *target, *source, *output;
  const char *target_arr, *source_arr;
}

template <class TMeshType>
int MergeArrays


template <class TMeshType>
int MyMain(int argc, char *argv[])
{
  if(argc < 6)
    return usage();

  const char *fn_ref = argv[1];
  const char *fn_out = argv[2];
  const char *nm_array = argv[3];
  const char *regex = argv[4];

  size_t ifarg = 5;

  // Set up the regular expression search
  RegularExpression re(regex);

  // Array of vtk objects
  std::vector<TMeshType *> mesh;

  // Read the reference object
  TMeshType *mref = ReadMesh<TMeshType>(fn_ref);

  // Search all filenames
  for(size_t i = ifarg; i < argc; i++)
    {
    // Read the file as a mesh
    TMeshType *m = ReadMesh<TMeshType>(argv[i]);
    if(!m || m->GetNumberOfPoints() == 0)
      cerr << "File " << argv[i] << " does not contain a valid mesh\n";
   
    // Process the regular expression
    re.find(argv[i]);
    printf("Matched %s in %s\n", re.match(1).c_str(), argv[i]);

    // Get the array from the file
    if(m->GetPointData()->HasArray(nm_array))
      {
      // Create the corresponding array in the reference mesh
      vtkDataArray *arr = m->GetPointData()->GetArray(nm_array);
      arr->SetName(re.match(1).c_str());
      mref->GetPointData()->AddArray(arr);
      }
    if(m->GetCellData()->HasArray(nm_array))
      {
      // Create the corresponding array in the reference mesh
      vtkDataArray *arr = m->GetCellData()->GetArray(nm_array);
      arr->SetName(re.match(1).c_str());
      mref->GetCellData()->AddArray(arr);
      }

    // Clear memory
    m->Delete();
    }

  // Save the reference mesh
  WriteMesh<TMeshType>(mref, fn_out);

  return 0;
}

int main(int argc, char *argv[])
{
  return MyMain<vtkUnstructuredGrid>(argc, argv);
}

