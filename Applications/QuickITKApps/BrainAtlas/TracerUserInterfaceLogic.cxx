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
