#include "TracerUserInterfaceLogic.h"
#include "FL/Fl_File_Chooser.H"

void
TracerUserInterfaceLogic
::OnMenuLoadMesh()
{
  char *file = fl_file_chooser("Select a VTK mesh", "*.vtk", NULL, 1);
  if(file)
    {
    m_Data->LoadInputMesh(file);
    m_WinTrace->OnMeshUpdate();
    }
}

void 
TracerUserInterfaceLogic
::OnMenuQuit()
{


}

void 
TracerUserInterfaceLogic
::OnButtonModeTrackball()
{
  m_WinTrace->SetMode(TracerMainWindow::TRACKBALL);
}

void 
TracerUserInterfaceLogic
::OnButtonModeTracer()
{
  m_WinTrace->SetMode(TracerMainWindow::TRACER);
}

void 
TracerUserInterfaceLogic
::OnButtonDeleteCurve()
{
  
}

void 
TracerUserInterfaceLogic
::OnButtonStartCurve()
{
  
}

void 
TracerUserInterfaceLogic
::OnSelectCurve()
{
  
}

