#ifndef __TracerUserInterfaceLogic_h_
#define __TracerUserInterfaceLogic_h_

#include "TracerUserInterface.h"
#include "TracerData.h"
#include "TracerMainWindow.h"

class TracerUserInterfaceLogic
: public TracerUserInterface, 
  virtual public ITracerDataListener
{
public:
  // Constructor
  TracerUserInterfaceLogic() : TracerUserInterface() 
    {
    m_Data = NULL;
    m_CurvesFile = "";
    m_CurvesDirty = false;
    }

  // Destructor
  virtual ~TracerUserInterfaceLogic() 
    { if(m_Data) m_Data->RemoveTracerDataListener(this); }
  
  // Callbacks from the user interface
  void OnMenuLoadMesh();
  void OnMenuLoadCurves();
  void OnMenuSaveCurves();
  void OnMenuSaveCurvesAs();
  void OnMenuQuit();
  void OnButtonModeTrackball();
  void OnButtonModeTracer();
  void OnButtonDeleteCurve();
  void OnButtonDeleteLastPoint();
  void OnButtonEditCurve();
  void OnButtonStartCurve();
  void OnInputSulcalFactor(double value);
  void OnSelectCurve();
  void OnSelectEdgeColoring(int value);
  void OnCheckDisplayEdges(int value);

  // Callbacks from TracerData event system
  void OnMeshChange(TracerDataEvent *evt);
  void OnFocusCurveChange(TracerDataEvent *evt);
  void OnFocusPointChange(TracerDataEvent *evt);
  void OnFocusCurveDataChange(TracerDataEvent *evt);
  void OnCurveListChange(TracerDataEvent *evt);
  void OnEdgeWeightsUpdate(TracerDataEvent *evt) {};

  void ShowWindows() 
    {
    m_WinMain->show();
    m_WinTrace->show();
    }

  /** Set the tracer data */
  void SetTracerData(TracerData *data) 
    {
    // Remove ourselves from the event list of the old data
    if(m_Data) m_Data->RemoveTracerDataListener(this);
    
    // Set the data
    m_Data = data;
    
    // Pass the data to the tracing GL window
    m_WinTrace->SetTracerData(data);

    // Set the event listeners
    m_Data->AddTracerDataListener(this);
    }
  
private:
  // Tracer data
  TracerData *m_Data;

  // Whether the curves are dirty (need saving) or clean
  bool m_CurvesDirty;

  // The name of the last saved curve file
  string m_CurvesFile;

  void RebuildCurveList();
  void DeactivateCurveEditControls();
  void ActivateCurveEditControls();
  void OnCurveStateChange();

  /** Check if it's ok to proceed with a destructive operation or
    * if the user wants to save the curves first */
  bool PromptForSaveOrCancel(const char *message);
};

#endif
