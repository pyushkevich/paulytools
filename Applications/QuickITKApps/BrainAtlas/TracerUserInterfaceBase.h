#ifndef __TracerUserInterfaceBase_h_
#define __TracerUserInterfaceBase_h_

class TracerUserInterfaceBase 
{
public:
  virtual void OnMenuLoadMesh() = 0;
  virtual void OnMenuQuit() = 0;
};

#endif
