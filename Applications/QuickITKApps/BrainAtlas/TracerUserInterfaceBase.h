#ifndef __TracerUserInterfaceBase_h_
#define __TracerUserInterfaceBase_h_

class TracerUserInterfaceBase 
{
public:
  virtual void OnMenuLoadMesh() = 0;
  virtual void OnMenuLoadCurves() = 0;
  virtual void OnMenuSaveCurves() = 0;
  virtual void OnMenuQuit() = 0;
  virtual void OnButtonModeTrackball() = 0;
  virtual void OnButtonModeTracer() = 0;
  virtual void OnButtonDeleteCurve() = 0;
  virtual void OnButtonEditCurve() = 0;
  virtual void OnButtonStartCurve() = 0;
  virtual void OnInputSulcalFactor(double value) = 0;
  virtual void OnSelectCurve() = 0;
  virtual void OnSelectEdgeColoring(int value) = 0;
  virtual void OnCheckDisplayEdges(int value) = 0;
};

#endif
