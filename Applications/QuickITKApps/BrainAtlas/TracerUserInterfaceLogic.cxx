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
::OnMenuSaveCurves()
{
  // Save the curves to a file
  char *file = fl_file_chooser("Select a text file to save the curves","*.txt","curves.txt",1);
  if(file)
    {
    m_Data->SaveCurves(file);
    }
}

void 
TracerUserInterfaceLogic
::OnMenuLoadCurves()
{
  // Save the curves to a file
  char *file = fl_file_chooser("Select a text file containing the curves","*.txt","curves.txt",1);
  if(file)
    {
    // m_Data->LoadCurves(file);
    m_WinTrace->SetTracerData(m_Data);
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
  if(fl_ask("Do you really want to delete\ncurve %s ?",
    m_Data->GetCurveName(m_Data->GetCurrentCurve())))
    {
    // Remove the currently selected curve
    m_Data->DeleteCurrentCurve();
    
    // Let the listener know of the curve update
    m_WinTrace->OnCurrentCurveChange();

    // Potentially, clear the list of names
    if(m_Data->GetNumberOfCurves() == 0)
      {
      // Enter trackball mode
      m_BtnModeTrackball->value(1);
      this->OnButtonModeTrackball();

      // Adjust controls
      DeactivateCurveEditControls();
      }
    else
      {
      // Rebuild the choice list
      RebuildCurveList();
      }    
    }
}

void 
TracerUserInterfaceLogic
::OnButtonEditCurve()
{
  // Update the name of the curve
  const char *name = 
    fl_input("Please enter the new name for the curve: ",
      m_Data->GetCurveName(m_Data->GetCurrentCurve()));
  
  // Change the name
  m_Data->SetCurveName(m_Data->GetCurrentCurve(),name);
  
  // Update drop box
  RebuildCurveList();
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
    ActivateCurveEditControls();
    }
}

void
TracerUserInterfaceLogic
::ActivateCurveEditControls()
{
  // Enable curve editing
  m_BtnModeTracer->activate();
  m_BtnDeleteCurve->activate();
  m_BtnEditCurve->activate();
  m_ChcCurve->activate();
  m_MenuSaveCurves->activate();
}

void
TracerUserInterfaceLogic
::DeactivateCurveEditControls()
{
  // Enable curve editing
  m_BtnModeTracer->deactivate();
  m_BtnDeleteCurve->deactivate();
  m_BtnEditCurve->deactivate();
  m_ChcCurve->deactivate();
  m_MenuSaveCurves->deactivate();
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
  // Get the currently selected curve
  int i = m_ChcCurve->value();
  
  // Select the curve in the data
  m_Data->SetCurrentCurve(i);
  
  // Pass the current curve information to the GL window
  m_WinTrace->OnCurrentCurveChange();
}

