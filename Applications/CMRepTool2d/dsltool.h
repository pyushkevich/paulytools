#ifndef _DSLTOOL_H_
#define _DSLTOOL_H_

#include "mreps2D.h"
#include "ispace.h"
#include "registry.h"
#include "Fl/gl.h"
#include "Fl/Fl_Widget.H"
#include "undo.h"
#include "likehood.h"

/***********************************************************
 * Global variables for the DSLTool application				  *
 ***********************************************************/
class ImageSpace;


// The image that is being displayed
// extern cimage baseImage;
// extern Array2D<GLbyte> baseBytes;

// Scale space built on top of the image
extern ImageSpace *imageSpace;

// Number of selected primitives in current model
extern int nSelected;

// The model that we will be working with
extern MGraph model;

// Settings of all sorts
extern Registry regOptions;

// Undo buffer, undo fuctions
extern UndoBuffer<MGraph> *undoBuffer;

// Medialness computer
extern IMatchComputer iMatchComputer;

void pushModel();
void updateUndoMenus();

// External functions
void uiSetZoom(double factor);
void uiNewModel(int,int,bool,double);
void uiSetCmd(Fl_Widget *widget);
void uiModelThin(double amount);

// Current tool mode enumeration and mode
enum DisplayModeEnum {
	DISPLAY_MREP=0,DISPLAY_SPLINE,DISPLAY_PEM,
		DISPLAY_BRANCH, DISPLAY_BEND};

extern DisplayModeEnum eDisplayMode;

// Enumeration of tool types
enum ToolEnum {
	TOOL_ZOOM = 0, TOOL_PAN, TOOL_SELECT, 
	TOOL_TRANSLATE, TOOL_SCALE, TOOL_ROTATE, 
	TOOL_STRETCH, TOOL_BEND, TOOL_SPLINE, 
	TOOL_PEM, TOOL_BRANCH
};

// Execute one of the tools
void uiRunTool(ToolEnum toolID);

extern ToolEnum eToolId;


// Called whenever model or selection change
void onModelOrSelectionChange();

// This method adds all kinds of additiocnal data to the model
void updatePrimitiveEditor(MGraph *model);
void updatePrimitiveEditor(MNode *node);

#endif
