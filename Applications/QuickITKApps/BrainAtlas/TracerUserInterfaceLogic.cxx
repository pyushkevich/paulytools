#include "TracerUserInterfaceLogic.h"
#include "FL/Fl_File_Chooser.H"

bool
TracerUserInterfaceLogic
::PromptForSaveOrCancel(const char *message)
{
  // If nothing needs saving, return true
  if(m_Data->GetCurves()->GetNumberOfCurves() == 0 || m_CurvesDirty == false)
    return true;

  // Check what the user wants
  int rc = fl_choice(message,
    "Cancel", "Save Curves","Discard Curves");

  // If 0, the user hit cancel, don't allow the requester to continue
  if(rc == 0) 
    return false;

  // If 2, the user hit don't save, allow to continue 
  else if(rc == 2)
    return true;

  // The user wants to save. Try that
  OnMenuSaveCurves();

  // Check that the curves were actually saved. If not, an error occurred, and
  // we can't let the requesting operation to continue
  return (m_CurvesDirty == false);
}

void
TracerUserInterfaceLogic
::OnMenuLoadMesh()
{
  // Check if we can continue with this destructive operation
  if(!PromptForSaveOrCancel(
    "After loading the mesh, all curves will be lost.\n"
    "Do you want to save the curves first?")) return;
  
  // Get the mesh file name
  char *file = fl_file_chooser("Select a VTK mesh", "*.vtk", NULL, 1);
  if(file)
    m_Data->LoadInputMesh(file);
}

void 
TracerUserInterfaceLogic
::OnMenuSaveCurves()
{
  // If there is a filename, save into it
  if(m_CurvesFile.size())
    {
    m_Data->SaveCurves(m_CurvesFile.c_str());
    m_CurvesDirty = false;
    }
  else
    {
    // Call the save-as code
    OnMenuSaveCurvesAs();
    }
}


void 
TracerUserInterfaceLogic
::OnMenuSaveCurvesAs()
{
  // Prompt for a file name using curves.txt as default
  char *file = fl_file_chooser(
    "Select a text file to save the curves","*.txt",
    m_CurvesFile.size() ? "curves.txt" : m_CurvesFile.c_str(), 1);
  
  // Do the actual save
  if(file) 
    {
    m_CurvesFile = file;
    OnMenuSaveCurves();
    }
}

void 
TracerUserInterfaceLogic
::OnMenuLoadCurves()
{
  // Check if we can continue with this destructive operation
  if(!PromptForSaveOrCancel(
    "Do you want to save the current curves before loading new ones?")) return;
  
  // Load curves from a file
  char *file = fl_file_chooser(
    "Select a text file containing the curves","*.txt",m_CurvesFile.c_str(),1);

  // Update the user interface
  if(file)
    {
    // Try loading
    if(m_Data->LoadCurves(file))
      { 
      m_CurvesDirty = false;
      m_CurvesFile = file;
      }
    }
}

void 
TracerUserInterfaceLogic
::OnMenuQuit()
{
  // Check if we can continue with this destructive operation
  if(!PromptForSaveOrCancel(
    "Do you want to save the curves before exiting?")) return;

  // Quit
  m_WinTrace->hide();
  m_WinMain->hide();
}

void 
TracerUserInterfaceLogic
::OnButtonModeTrackball()
{
  m_BtnModeTrackball->value(1);
  m_BtnModeTracer->value(0);
  m_MenuModeTrackball->setonly();
  
  m_WinTrace->SetMode(TracerMainWindow::TRACKBALL);
}

void 
TracerUserInterfaceLogic
::OnButtonModeTracer()
{
  m_BtnModeTrackball->value(0);
  m_BtnModeTracer->value(1);
  m_MenuModeTracer->setonly();
  
  m_WinTrace->SetMode(TracerMainWindow::TRACER);
}

void
TracerUserInterfaceLogic
::OnButtonDeleteLastPoint()
{
  // Delete the last point
  m_Data->DeleteCurrentPoint();
}

void 
TracerUserInterfaceLogic
::OnButtonDeleteCurve()
{
  // Get the curves info
  const TracerCurves *curves = m_Data->GetCurves();
  
  if(fl_ask("Do you really want to delete\ncurve %s ?",
    curves->GetCurveName(m_Data->GetCurrentCurve())))
    {
    // Remove the currently selected curve
    m_Data->DeleteCurrentCurve();
    }
}

void 
TracerUserInterfaceLogic
::OnButtonEditCurve()
{
  // Update the name of the curve
  const char *name = 
    fl_input("Please enter the new name for the curve: ",
      m_Data->GetCurves()->GetCurveName(m_Data->GetCurrentCurve()));
  
  // Change the name
  m_Data->SetCurveName(m_Data->GetCurrentCurve(),name);
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
    }
}

void
TracerUserInterfaceLogic
::OnCurveStateChange()
{
}

void
TracerUserInterfaceLogic
::ActivateCurveEditControls()
{
  // Enable curve editing
  m_MenuSaveCurves->activate();
  m_MenuSaveCurvesAs->activate();
  m_MenuModeTracer->activate();
  m_MenuDeleteCurve->activate();
  m_MenuEditCurve->activate();
  
  m_BtnModeTracer->activate();
  m_BtnDeleteCurve->activate();
  m_BtnEditCurve->activate();
}

void
TracerUserInterfaceLogic
::DeactivateCurveEditControls()
{
  // Disable curve editing
  m_BtnModeTracer->deactivate();
  m_BtnDeleteCurve->deactivate();
  m_BtnEditCurve->deactivate();
  
  m_MenuSaveCurves->deactivate();
  m_MenuSaveCurvesAs->deactivate();
  m_MenuModeTracer->deactivate();
  m_MenuDeleteCurve->deactivate();
  m_MenuEditCurve->deactivate();
}

void
TracerUserInterfaceLogic
::RebuildCurveList()
{
  // Clear the choice of curves
  m_ChcCurve->clear();

  // If there are no curves, disable the drop down
  if(m_Data->GetCurves()->GetNumberOfCurves() > 0)
    {
    // Get the alphabetical curve list
    typedef TracerCurves::StringIdPair StringIdPair;
    typedef list<StringIdPair> StringIdList;
    StringIdList lCurveNames;
    m_Data->GetCurves()->GetAlphabeticCurveList(lCurveNames);

    // Keep track of which item is current
    int iCurrent = -1;

    // Add the curves to the choice
    StringIdList::iterator it=lCurveNames.begin();
    while(it != lCurveNames.end())
      {
      int iPos = 
        m_ChcCurve->add(it->first.c_str(), 0, NULL, (void *) it->second, 0);

      if(m_Data->GetCurrentCurve() == it->second)
        iCurrent = iPos; 

      ++it;
      }

    // Set the current item
    m_ChcCurve->value(iCurrent);
    m_ChcCurve->activate();
    }
  else
    {
    m_ChcCurve->deactivate();
    }
}

void 
TracerUserInterfaceLogic
::OnSelectCurve()
{
  // Get the currently selected curve
  int i = m_ChcCurve->value();

  // Make sure this is a valid item
  if(i < 0 || i >= m_ChcCurve->size()) 
    return;
  
  // Get the user data for the selected curve
  TracerCurves::IdType id = 
    (TracerCurves::IdType) m_ChcCurve->menu()[i].user_data_;
  
  // Select the curve in the data
  m_Data->SetCurrentCurve(id);
  
  // Pass the current curve information to the GL window
  OnCurveStateChange();
}

void
TracerUserInterfaceLogic
::OnInputSulcalFactor(double value)
{
  if(value == 0)
    m_Data->SetEdgeWeightsToEuclideanDistance();
  else
    m_Data->SetEdgeWeightsToPitchDistance(value);
}

void 
TracerUserInterfaceLogic
::OnSelectEdgeColoring(int value)
{
  if(value == 0)
    m_WinTrace->SetEdgeDisplayMode(TracerMainWindow::EDGE_DISPLAY_DISTANCE);
  else if(value == 1)
    m_WinTrace->SetEdgeDisplayMode(TracerMainWindow::EDGE_DISPLAY_LENGTH);
  else if(value == 2)
    m_WinTrace->SetEdgeDisplayMode(TracerMainWindow::EDGE_DISPLAY_PLAIN);
}

void 
TracerUserInterfaceLogic
::OnCheckDisplayEdges(int value)
{
  if(value)
    {
    // Enable the edge coloring drop box
    m_ChcEdgeColoring->activate();

    // Enable edge display in the window based on current
    // edge coloring
    OnSelectEdgeColoring(m_ChcEdgeColoring->value());
    }
  else
    {
    // Disable edge coloring choice
    m_ChcEdgeColoring->deactivate();

    // Disable edge display
    m_WinTrace->SetEdgeDisplayMode(TracerMainWindow::EDGE_DISPLAY_NONE);
    }
}

void 
TracerUserInterfaceLogic
::OnMeshChange(TracerDataEvent *evt)
{
  // If the mesh exists, draw the mesh
  if( m_Data->IsMeshLoaded() )
    {
    m_BtnStartCurve->activate();
    m_MenuStartCurve->activate();
    m_MenuLoadCurves->activate();
    }
  else
    {
    m_BtnStartCurve->deactivate();
    m_MenuStartCurve->deactivate();
    m_MenuLoadCurves->deactivate();
    }
}

void 
TracerUserInterfaceLogic
::OnFocusCurveChange(TracerDataEvent *evt)
{
  // If there is no current curve, disable controls
  if(!m_Data->IsCurrentCurveValid())
    {
    // Enter trackball mode
    this->OnButtonModeTrackball();

    // Adjust controls
    DeactivateCurveEditControls();
    }
  else
    {
    // Activate the controls
    ActivateCurveEditControls();
    }

  // Change the current item in the drop-down
  RebuildCurveList();
}

void 
TracerUserInterfaceLogic
::OnFocusPointChange(TracerDataEvent *evt)
{
  // Check if the delete last point is valid
  if(m_Data->IsPathSourceSet())
    {
    m_BtnDeleteLastPoint->activate();
    m_MenuDeleteLastPoint->activate();
    }
  else
    {
    m_BtnDeleteLastPoint->deactivate();
    m_MenuDeleteLastPoint->deactivate();
    }
}

void 
TracerUserInterfaceLogic
::OnFocusCurveDataChange(TracerDataEvent *evt)
{
  m_CurvesDirty = true;
}

void
TracerUserInterfaceLogic
::OnCurveListChange(TracerDataEvent *evt)
{
  // The list of curves has changed, so we need to rebuild
  // the curve drop down box. We don't worry about the buttons
  // because that's handled by the Focus and Data events
  RebuildCurveList();

  // If there is more than one curve, enable the saving and loading
  if(m_Data->GetCurves()->GetNumberOfCurves())
    {
    m_MenuSaveCurvesAs->activate();
    m_MenuSaveCurves->activate();
    }
  else
    {
    m_MenuSaveCurvesAs->deactivate();
    m_MenuSaveCurves->deactivate();
    }

  m_CurvesDirty = true;
}


