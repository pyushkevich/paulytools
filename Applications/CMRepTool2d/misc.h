#ifndef _DSL_MISC_H_
#define _DSL_MISC_H_

#include "matrix.h"
#include <mreps2D.h>
#include "glutil.h"
#include "scripting.h"
#include <Fl/gl.h>

struct MouseEvent;
class ImageSpace;

// Miscellaneous functions

// Sine and cosine in degreed
// For PI conversion
inline double dsin(double x) {
	return sin(M_PI*x/180.0);
}
inline double dcos(double x) {
	return cos(M_PI*x/180.0);
}

double computeRadius(const MAtom &atom);
void computeSelectionBox();
bool locateBAtom(MouseEvent &e,BRef &outBRef);

// We will have a script created for UI commands
extern UserInterfaceScript uiScript;

extern const int TARGET_CENTER;
extern const int TARGET_NOSE;
extern const int TARGET_LEFT;
extern const int TARGET_RIGHT;
extern const double DEG_TO_RAD;

// Super helpful little cheat
void drawSingleAtom(const MAtom &a);
void draw(MNode *);
void drawPrimitiveEditMode(MNode *);
void draw(MGraph *graph,bool fast,bool pem);
void drawBranchMode(MGraph *graph);

// Return a single selected figure and deselect everything else
MFigure *getSingleSelectedFigure();

#endif
