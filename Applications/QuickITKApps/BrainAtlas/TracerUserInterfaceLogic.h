#ifndef __TracerUserInterfaceLogic_h_
#define __TracerUserInterfaceLogic_h_

#include "FLWidgetActivationManager.h"
#include "TracerUserInterface.h"
#include "TracerData.h"
#include "TracerMainWindow.h"

class TracerUserInterfaceLogic
: public TracerUserInterface, 
  virtual public ITracerDataListener
{
public:
  // Constructor
  TracerUserInterfaceLogic();

  // Destructor
  virtual ~TracerUserInterfaceLogic() 
    { if(m_Data) m_Data->RemoveTracerDataListener(this); }

  /** The MakeWindow method, in order to perform actions after the 
    * controls in the window have been initialized */
  virtual void MakeWindow();
  
  // Callbacks from the user interface
  void OnMenuLoadMesh();
  void OnMenuLoadCurves();
  void OnMenuSaveCurves();
  void OnMenuSaveCurvesAs();
  void OnMenuQuit();

  void OnButtonModeTrackball();
  void OnButtonModeTracer();
  void OnButtonModeMarker();

  void OnButtonDeleteCurve();
  void OnButtonDeleteLastPoint();
  void OnButtonEditCurve();
  void OnButtonStartCurve();
  void OnSelectCurve();

  void OnButtonRenameMarker();
  void OnButtonRecolorMarker();
  void OnButtonDeleteMarker();
  void OnButtonComputeRegions();
  void OnSelectMarker(int value);
 
  void OnCheckDisplayEdges(int value);
  void OnSelectEdgeColoring(int value);
  void OnCheckCenterMesh(int value);
  void OnCheckDisplayNeighborhood(int value);
  void OnInputNeighborhoodRadius(double value);
  void OnInputSulcalFactor(double value);

  // Callbacks from TracerData event system
  void OnMeshChange(TracerDataEvent *evt);
  void OnFocusCurveChange(TracerDataEvent *evt);
  void OnFocusPointChange(TracerDataEvent *evt);
  void OnFocusCurveDataChange(TracerDataEvent *evt);
  void OnCurveListChange(TracerDataEvent *evt);
  void OnEdgeWeightsUpdate(TracerDataEvent *evt) {};
  void OnMarkerListChange(TracerDataEvent *evt);
  void OnFocusMarkerChange(TracerDataEvent *evt);
  void OnSegmentationChange(TracerDataEvent *evt) {};

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
  
  /** State flags used in conjunction with the activation manager */
  enum ActivationFlag {
    AF_NULL, AF_MESH, AF_TRACER_DATA, AF_CURVES_EXIST, AF_MARKERS_EXIST, 
    AF_ACTIVE_CURVE, AF_ACTIVE_CONTROL, AF_ACTIVE_MARKER, 
    AF_TRACER_MODE, AF_MARKER_MODE };

  /** Control activation manager used to enable and disable groups
    of controls based on a set of flags */
  FLWidgetActivationManager<ActivationFlag> m_Activation;

  /** Check if it's ok to proceed with a destructive operation or
    * if the user wants to save the curves first */
  bool PromptForSaveOrCancel(const char *message);
};

#endif
