#include "TracerUserInterfaceLogic.h"
#include "FL/Fl_File_Chooser.H"
#include "FL/Fl_Color_Chooser.H"

// Constructor
TracerUserInterfaceLogic 
::TracerUserInterfaceLogic() 
: TracerUserInterface() 
{
  m_Data = NULL;
  m_CurvesFile = "";
  m_CurvesDirty = false;

  // Set up the dependences between activation flags
  m_Activation.SetFlagImplies(AF_TRACER_DATA, AF_MESH);
  m_Activation.SetFlagImplies(AF_CURVES_EXIST, AF_TRACER_DATA);
  m_Activation.SetFlagImplies(AF_MARKERS_EXIST, AF_TRACER_DATA);
  m_Activation.SetFlagImplies(AF_ACTIVE_CURVE, AF_CURVES_EXIST);
  m_Activation.SetFlagImplies(AF_ACTIVE_CONTROL, AF_ACTIVE_CURVE);
  m_Activation.SetFlagImplies(AF_ACTIVE_MARKER, AF_MARKERS_EXIST);
  m_Activation.SetFlagImplies(AF_TRACER_MODE, AF_ACTIVE_CURVE);
  m_Activation.SetFlagImplies(AF_MARKER_MODE, AF_MESH);
}

void
TracerUserInterfaceLogic
::MakeWindow()
{
  // Call the parent method that creates all the controls
  TracerUserInterface::MakeWindow();

  // Set up the state-based control activation
  m_Activation.AddMenuItem(m_MenuLoadCurves, AF_MESH); 
  m_Activation.AddMenuItem(m_MenuSaveCurves, AF_TRACER_DATA);
  m_Activation.AddMenuItem(m_MenuSaveCurvesAs, AF_TRACER_DATA ); 
  m_Activation.AddMenuItem(m_MenuImportCurves, AF_MESH); 
  m_Activation.AddMenuItem(m_MenuExportCurves, AF_TRACER_DATA);
  m_Activation.AddMenuItem(m_MenuModeTracer, AF_ACTIVE_CURVE); 
  m_Activation.AddMenuItem(m_MenuModeMarker, AF_MESH);
  m_Activation.AddMenuItem(m_MenuStartCurve, AF_MESH); 
  m_Activation.AddMenuItem(m_MenuEditCurve, AF_ACTIVE_CURVE); 
  m_Activation.AddMenuItem(m_MenuDeleteCurve,    AF_ACTIVE_CURVE); 
  m_Activation.AddMenuItem(m_MenuDeleteLastPoint, AF_ACTIVE_CONTROL); 
  m_Activation.AddMenuItem(m_MenuRenameMarker, AF_ACTIVE_MARKER);
  m_Activation.AddMenuItem(m_MenuDeleteMarker, AF_ACTIVE_MARKER);
  m_Activation.AddMenuItem(m_MenuRecolorMarker, AF_ACTIVE_MARKER);
  m_Activation.AddMenuItem(m_MenuComputeRegions, AF_MARKERS_EXIST);
  
  m_Activation.AddWidget(m_ChcCurve, AF_TRACER_DATA); 
  m_Activation.AddWidget(m_BtnStartCurve, AF_MESH); 
  m_Activation.AddWidget(m_BtnDeleteCurve, AF_ACTIVE_CURVE); 
  m_Activation.AddWidget(m_BtnEditCurve, AF_ACTIVE_CURVE); 
  m_Activation.AddWidget(m_BtnDeleteLastPoint, AF_ACTIVE_CONTROL); 
  m_Activation.AddWidget(m_BtnModeTracer, AF_ACTIVE_CURVE); 
  m_Activation.AddWidget(m_BtnModeMarker, AF_MESH);
  m_Activation.AddWidget(m_BtnRenameMarker, AF_ACTIVE_MARKER);
  m_Activation.AddWidget(m_BtnDeleteMarker, AF_ACTIVE_MARKER);
  m_Activation.AddWidget(m_BtnRecolorMarker, AF_ACTIVE_MARKER);
  m_Activation.AddWidget(m_BrsMarkers, AF_MARKERS_EXIST);
  m_Activation.AddWidget(m_BtnComputeRegions, AF_MARKERS_EXIST);
}

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
::OnMenuExportCurves()
{
  // Prompt for a file name using curves.txt as default
  char *file = fl_file_chooser(
    "Select a text file to export the curves","*.txt",
    "curves_export.txt", 1);
  
  // Do the actual save
  if(file) 
    m_Data->ExportCurves(file);
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
::OnMenuImportCurves()
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
    if(m_Data->ImportCurves(file))
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
  m_BtnModeMarker->value(0);
  m_MenuModeTrackball->setonly();
  
  m_WinTrace->SetMode(TracerMainWindow::TRACKBALL);
}

void 
TracerUserInterfaceLogic
::OnButtonModeTracer()
{
  m_BtnModeTrackball->value(0);
  m_BtnModeTracer->value(1);
  m_BtnModeMarker->value(0);
  m_MenuModeTracer->setonly();
  
  m_WinTrace->SetMode(TracerMainWindow::TRACER);
}

void
TracerUserInterfaceLogic
::OnButtonModeMarker()
{
  m_BtnModeTrackball->value(0);
  m_BtnModeTracer->value(0);
  m_BtnModeMarker->value(1);
  m_MenuModeMarker->setonly();
  
  m_WinTrace->SetMode(TracerMainWindow::MARKER);
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
::OnCheckCenterMesh(int value)
{
  // Tell the window to center mesh when possible
  m_WinTrace->SetCenterMesh(value != 0);
}
  
void 
TracerUserInterfaceLogic
::OnCheckDisplayNeighborhood(int value)
{
  if(value)
    {
    m_WinTrace->SetSurfaceDisplayMode(
      TracerMainWindow::SURFACE_DISPLAY_NEIGHBORHOOD);
    m_WinTrace->SetNeighborhoodSize(m_InNeighborhoodSize->value());
    m_InNeighborhoodSize->activate();
    }
  else
    {
    m_WinTrace->SetSurfaceDisplayMode(
      TracerMainWindow::SURFACE_DISPLAY_ALL);
    m_InNeighborhoodSize->deactivate();
    }
}
  
void 
TracerUserInterfaceLogic
::OnInputNeighborhoodRadius(double value)
{
  m_WinTrace->SetNeighborhoodSize(value);
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
  if(value == 2)
    m_WinTrace->SetEdgeDisplayMode(TracerMainWindow::EDGE_DISPLAY_DISTANCE);
  else if(value == 1)
    m_WinTrace->SetEdgeDisplayMode(TracerMainWindow::EDGE_DISPLAY_LENGTH);
  else if(value == 0)
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
  // Activate the controls that are tied to the mesh being present
  m_Activation.UpdateFlag( AF_MESH, m_Data->IsMeshLoaded() );
}

void 
TracerUserInterfaceLogic
::OnFocusCurveChange(TracerDataEvent *evt)
{
  // Update the flags based on current state
  m_Activation.UpdateFlag( AF_ACTIVE_CURVE, m_Data->IsCurrentCurveValid());
  
  // If there is no current curve, disable controls
  if(!m_Data->IsCurrentCurveValid())
    this->OnButtonModeTrackball();

  // Change the current item in the drop-down
  RebuildCurveList();
}

void 
TracerUserInterfaceLogic
::OnFocusPointChange(TracerDataEvent *evt)
{
  // Check if the delete last point is valid
  m_Activation.UpdateFlag(AF_ACTIVE_CONTROL, m_Data->IsPathSourceSet());
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
  // Update the active and inactive controls
  m_Activation.UpdateFlag( AF_CURVES_EXIST, 
    m_Data->GetCurves()->GetNumberOfCurves() > 0);
  
  // The list of curves has changed, so we need to rebuild
  // the curve drop down box. We don't worry about the buttons
  // because that's handled by the Focus and Data events
  RebuildCurveList();

  // Set the dirty flag
  m_CurvesDirty = true;
}

void
TracerUserInterfaceLogic
::OnButtonRenameMarker()
{
  // The marker id
  vtkIdType id = m_Data->GetCurrentMarker();
  
  // Show the color chooser
  const char *name = 
    fl_input("Please enter the name of the marker", m_Data->GetMarkerName(id));

  if(name) m_Data->RenameMarker(id, name);
}

void
TracerUserInterfaceLogic
::OnButtonRecolorMarker()
{
  // The marker id
  vtkIdType id = m_Data->GetCurrentMarker();
  
  // Get the current color of the marker
  TracerData::Vec x = m_Data->GetMarkerColor(id);

  // Show the color chooser
  if(fl_color_chooser("Please select the marker color", x[0],x[1],x[2]))
    m_Data->RecolorMarker(id, x[0],x[1],x[2]);
}

void
TracerUserInterfaceLogic
::OnButtonDeleteMarker()
{
  // Must have an active marker
  assert(m_Data->IsCurrentMarkerValid());
  
  // Delete the active marker
  m_Data->DeleteMarker(m_Data->GetCurrentMarker());
}

void
TracerUserInterfaceLogic
::OnButtonComputeRegions()
{
  m_Data->ComputeMarkerSegmentation();
}

void
TracerUserInterfaceLogic
::OnSelectMarker(int value)
{
  // Get the id associated with the marker
  vtkIdType id = (vtkIdType) m_BrsMarkers->data(value);

  // Set the id as active id
  m_Data->SetCurrentMarker(id);
}

void 
TracerUserInterfaceLogic::
OnMarkerListChange(TracerDataEvent *evt)
{
  // Rebuild the list of markers
  m_BrsMarkers->clear();

  // Get a list of marker ID's
  list<vtkIdType> lMarkers;
  m_Data->GetMarkerIds(lMarkers);

  // Populate the list of markers
  list<vtkIdType>::iterator it = lMarkers.begin();
  while(it != lMarkers.end())
    {
    m_BrsMarkers->add( m_Data->GetMarkerName(*it), (void *)(*it) );
    ++it;
    }
  
  // Enable/disable controls
  m_Activation.UpdateFlag(AF_MARKERS_EXIST, m_Data->GetNumberOfMarkers() > 0);

  // Select the active item (use method below)
  OnFocusMarkerChange(evt);
}

void 
TracerUserInterfaceLogic::
OnFocusMarkerChange(TracerDataEvent *evt)
{
  // Deselect all items in the browser
  m_BrsMarkers->deselect();
  
  // Update the flags
  m_Activation.UpdateFlag(AF_ACTIVE_MARKER, m_Data->IsCurrentMarkerValid());

  // Select the active item
  if(m_Data->IsCurrentMarkerValid())
    {
    
    // Select the active item
    for(unsigned int i = 1; i <= m_BrsMarkers->size(); i++)
      {
      if(((vtkIdType)m_BrsMarkers->data(i)) == m_Data->GetCurrentMarker())
        {
        m_BrsMarkers->select(i);
        break;
        }
      }
    }
}
