#ifndef __TracerUserInterfaceLogic_h_
#define __TracerUserInterfaceLogic_h_

#include "TracerUserInterface.h"
#include "TracerMainWindow.h"

class TracerUserInterfaceLogic : public TracerUserInterface 
{
public:
  void OnMenuLoadMesh();
  void OnMenuLoadCurves();
  void OnMenuSaveCurves();
  void OnMenuQuit();
  void OnButtonModeTrackball();
  void OnButtonModeTracer();
  void OnButtonDeleteCurve();
  void OnButtonEditCurve();
  void OnButtonStartCurve();
  void OnInputSulcalFactor(double value);
  void OnSelectCurve();

  void ShowWindows() 
    {
    m_WinMain->show();
    m_WinTrace->show();
    }

  /** Set the tracer data */
  void SetTracerData(TracerData *data) 
    {
    m_WinTrace->SetTracerData(data);
    m_Data = data;
    }
  
private:
  // Tracer data
  TracerData *m_Data;

  void RebuildCurveList();
  void DeactivateCurveEditControls();
  void ActivateCurveEditControls();
};

#endif
