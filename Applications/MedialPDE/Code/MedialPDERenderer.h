#ifndef __MedialPDERenderer_h_
#define __MedialPDERenderer_h_

#include "glengine.h"

class MedialPDESolver;

void glDrawWireframeElements(unsigned short width, unsigned short height);
void glDrawQuadStripElements(unsigned short width,unsigned short height);

class PDESplineRenderer : public GLDisplayListRenderer
{
public:
  // Constructor
  PDESplineRenderer(MedialPDESolver *solver);

  // Build the display list
  virtual void build();

private:
  // INternal methods
  void DrawInternalPoints( size_t nCuts );

  MedialPDESolver *solver;
  GLMaterial *matMedial, *matBoundary;
};

#endif
