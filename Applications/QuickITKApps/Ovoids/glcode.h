/************************************************************************
 * COMP 257 Final Project
 * Differential Geometry of Implicit Surfaces
 * Author:		Paul Yushkevich
 * Module:
 * Last Update: Dec 13, 1998
 *
 * Description:
 *
 *
 *************************************************************************/

#ifndef _GLCODE_H_
#define _GLCODE_H_

#include "implicit.h"
#include "blobmodl.h"
#include <FL/Fl_Window.h>

int isurfCallback(int v1,int v2,int v3,IMP_VERTICES vv);
int isurfCallbackColored(int v1,int v2,int v3,IMP_VERTICES vv);

void initGLDisplay();

extern Fl_Window *win3d;

// A list of displayed ISurfaces.  Add surfaces you need drawn onto this list
extern PList *iSurfList;

// The current point
extern PointData *currentPoint;

#endif
