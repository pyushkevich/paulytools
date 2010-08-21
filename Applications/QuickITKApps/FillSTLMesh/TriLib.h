#ifndef __TriLib_h_
#define __TriLib_h_
#define TRILIBRARY
#define ANSI_DECLARATORS

#include "triangle.h"
// Write triangulateio object to VTK
void WriteTriangleOutputAsPolyDataMesh(triangulateio &out, const char *fname);
  

#endif
