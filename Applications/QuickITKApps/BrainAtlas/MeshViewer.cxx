/**
 * This is a trivial VTK mesh viewer
 */
#include <iostream>

#include "vtkPolyDataReader.h"
#include "vtkPolyDataMapper.h"
#include "vtkOBJReader.h"
#include "vtkBYUReader.h"
#include "vtkLODActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkStripper.h"


using namespace std;


int main(int argc, char *argv[])
{
  if(argc < 2)
    {
    cout << "usage: [options] MeshViewer mesh.vtk" << endl;
    cout << "options:" << endl;
    cout << "   -a NAME         Specify the name of a scalar array to overlay " << endl;
    cout << "   -s X.XX X.XX    Specify the scalar range to show the array " << endl;
    return -1;
    }

  // Parameter scanning
  const char *fnInput = argv[argc - 1];
  const char *sScalarArray = NULL;
  double xScalarMin, xScalarMax;
  bool flagSetScalars = false;
  for(int iArg = 1; iArg < argc-1; iArg++)
    {
    if(!strcmp(argv[iArg],"-a"))
      {
      sScalarArray = argv[++iArg];
      }
    else if(!strcmp(argv[iArg],"-s"))
      {
      xScalarMin = atof(argv[++iArg]);
      xScalarMax = atof(argv[++iArg]);
      flagSetScalars = true;
      }
    }
      
  // Load the mesh
  vtkPolyData *poly = NULL;
  if(0 != strstr(fnInput, ".byu"))
    {
    vtkBYUReader *byu = vtkBYUReader::New();
    byu->SetFileName(fnInput);
    byu->Update();
    poly = byu->GetOutput();
    }
  else if(0 != strstr(fnInput, ".obj"))
    {
    vtkOBJReader *obj = vtkOBJReader::New();
    obj->SetFileName(fnInput);
    obj->Update();
    poly = obj->GetOutput();
    }
  else if(0 != strstr(fnInput, ".vtk"))
    {
    vtkPolyDataReader *fltReader = vtkPolyDataReader::New();
    fltReader->SetFileName(fnInput);
    fltReader->Update();
    poly = fltReader->GetOutput();
    }
  else
    {
    cerr << "Unknown format!" << endl; return -1;
    }

  // Triangulate the mesh
  vtkStripper *fltStripper = vtkStripper::New();
  fltStripper->SetInput(poly);
  fltStripper->Update();
  cout << "Converted to triangle strips " << endl;

  // Check what additional information the mesh has
  vtkPolyData *mesh = fltStripper->GetOutput();

  // If the scalar arrays are specified ...
  if(sScalarArray) 
    {
    // Check how many arrays there are
    unsigned int nArrays = mesh->GetPointData()->GetNumberOfArrays();
    cout << "Mesh has " << nArrays << " arrays " << endl;
    
    // List arrays 
    for(unsigned int i = 0; i < nArrays; i++)
      {
      cout << " Array " << i << "\tnamed " << mesh->GetPointData()->GetArrayName(i) << "\thas "
        << mesh->GetPointData()->GetArray(i)->GetNumberOfTuples() << " tuples " << endl;
      }

    // Assign the selected array to scalars
    vtkDataArray *array = mesh->GetPointData()->GetArray(sScalarArray);
    if(array)
      mesh->GetPointData()->SetScalars(array);
    else
      cout << "Your array could not be found!" << endl;
    }
  
  // Draw the mesh
  vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
  mapper->SetInput(mesh);

  if(flagSetScalars)
    mapper->SetScalarRange(xScalarMin,xScalarMax);

  vtkLODActor *actor = vtkLODActor::New();
  actor->SetMapper(mapper);

  vtkRenderer *ren = vtkRenderer::New();
  ren->AddActor(actor);
  ren->SetBackground(0.1,0.2,0.4);

  vtkRenderWindow *renWin = vtkRenderWindow::New();
  renWin->AddRenderer(ren);
  renWin->SetSize(500,500);

  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);
  iren->Initialize();

  renWin->Render();
  iren->Start();

  iren->Delete();
  renWin->Delete();
  ren->Delete();
  actor->Delete();
  mapper->Delete();


}

