/**
 * This is a trivial VTK mesh viewer
 */
#include <iostream>

#include "vtkPolyDataReader.h"
#include "vtkPolyDataMapper.h"
#include "vtkLODActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

using namespace std;


int main(int argc, char *argv[])
{
  cout << "usage: MeshViewer mesh.vtk" << endl;

  // Load the mesh
  vtkPolyDataReader *fltReader = vtkPolyDataReader::New();
  fltReader->SetFileName(argv[1]);
  fltReader->Update();
  
  // Draw the mesh
  vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
  mapper->SetInput(fltReader->GetOutput());

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

