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
::OnButtonEditCurve()
{
  
}

void 
TracerUserInterfaceLogic
::OnButtonStartCurve()
{
  // Prompt the user for the name of the curve
  const char *name = fl_input(
    "Please enter the name of the new curve: ");
  if(name)
    {
    // Add the curve
    m_Data->AddNewCurve(name);

    // Rebuild the curve list
    RebuildCurveList();

    // Enable curve editing
    m_BtnModeTracer->activate();
    m_BtnDeleteCurve->activate();
    m_BtnEditCurve->activate();
    m_ChcCurve->activate();
    }
}

void
TracerUserInterfaceLogic
::RebuildCurveList()
{
  // Clear the choice of curves
  m_ChcCurve->clear();

  // Get all the curves
  for(unsigned int iCurve=0;iCurve<m_Data->GetNumberOfCurves();iCurve++)
    m_ChcCurve->add(m_Data->GetCurveName(iCurve));

  // Set the current item
  m_ChcCurve->value(m_Data->GetCurrentCurve());
}

void 
TracerUserInterfaceLogic
::OnSelectCurve()
{
  
}

