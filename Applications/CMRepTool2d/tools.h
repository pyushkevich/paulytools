#ifndef _TOOLS_H_
#define _TOOLS_H_

/******************************************************************
* DSLTOOL                                                        *
******************************************************************
* Author:                 Paul Yushkevich
*
* Date:                   July 18, 1999
*
* Description             Tool handlers for DSLTool
*									
*	Sources:                
*									
* Dependencies:           matrix,dsl2d,fltk
******************************************************************
* tools.h
*	---------
* This module contains the handlers for all the tools in DSLTool.
******************************************************************/
#include <matrix.h>
#include <mreps2D.h>
// #include "SplineMrep.h"

class Fl_Return_Button;

/******************************************************************
* Support Structures                                             *
******************************************************************/
// This structure represents a point on the path drawn by the user.
struct BendPathPoint {
	// Position
	Vector2D x;

	// Slope
	double t;

	// Distance from beginning
	double s;
};

// A small helper event object to be passed on to tool handlers 
struct MouseEvent {
	// Starting mouse position and ending mouse position in space
	Vector2D sSpace;
	Vector2D xSpace;

	// Starting mouse position and ending mouse position on the screen
	Vector2D sScreen;
	Vector2D xScreen;

	// Same, but unnormalized by h,w
	Vector2D sRawScreen;
	Vector2D xRawScreen;

	// Target figure and primitive
	BRef target;

	// Shift/Ctrl/Alt flags
	bool shift,ctrl,alt;

	// Button tyoe
	int button;

	// Drag length - total distance that the user has dragged the mouse
	double dragLength;

	// Click flag - set to true if dragLength is less than epsilon
	bool click;

	// Default constructor
	MouseEvent() {
		shift = ctrl = alt = false;
		dragLength = 0;
		click = false;
	}
};

// Classes defined elsewhere
class DisplayWindow;

/******************************************************************
* Generic Tool Handler                                           *
* --------------------                                           *
* This class is the parent for all other tool handlers           *
******************************************************************/
class ToolHandler {
protected:
	DisplayWindow &dw;
public:
	// Respond to user's mouse actions
	virtual void press(MouseEvent &ev) {};
	virtual void release(MouseEvent &ev) {};
	virtual void drag(MouseEvent &ev) {};

	// Called by window's draw routine - paint tool-specific stuff
	virtual void display() {};

	// Called when tool is assigned to / unassigned from the window
	virtual void init() {};
	virtual void close() {};

	// Get the name and context help for the tool
	virtual const char *getLabel() = 0;
	virtual const char *getContextHelp() = 0;

	// Constructor that takes a window as a parameter
	ToolHandler(DisplayWindow &win) : dw(win) {
	}

	// Virtual destructor is needed because we have virtual members
	virtual ~ToolHandler() {};
};

/******************************************************************
* Zoom and Pan Tool Handler                                      *
******************************************************************/
class ZoomPanToolHandler : public ToolHandler {
	// Window size/pos when zoom started
	double initWinSize;
	Vector2D initWinPosition;

	bool zoomMode;

public:
	void press(MouseEvent &ev);
	void release(MouseEvent &ev);
	void drag(MouseEvent &ev);
	void close();

	// Get the label name for this tool
	const char *getLabel() {
		return zoomMode ? "Zoom Tool" : "Panning Tool";
	}

	// Get the context help for the tool
	const char *getContextHelp() {
		return zoomMode 
			? "Left-drag to zoom in/out on a point on the image.  Right-drag to pan image" 
			: "Left-drag to pan image.  Right-drag to zoom in/out on a point on the image.";
	}

	ZoomPanToolHandler(DisplayWindow &win,bool zoomMode);
};


/******************************************************************
* Selection Tool Handler                                         *
******************************************************************/
class SelectionToolHandler : public ToolHandler {
	MouseEvent lastEvent;
	void processSelectionClick(MouseEvent &e);
	void processSelectionDrag(MouseEvent &e);

public:
	void init();
	void press(MouseEvent &ev);
	void release(MouseEvent &ev);
	void drag(MouseEvent &ev);
	void display();

	// Get the label name for this tool
	char *getLabel() {
		return "Selection Tool";
	}

	// Get the context help for the tool
	char *getContextHelp() {
		return "Click on a primitive once to select.  Click twice to select figure.  Drag mouse to select primitives "
			"inside a box.  Use SHIFT to add to selection, CTRL to toggle selection.";
	}

	SelectionToolHandler(DisplayWindow &win) : ToolHandler(win) {
	}
};


/******************************************************************
* Translation Tool Handler                                       *
******************************************************************/
class TranslateToolHandler : public ToolHandler {
	// The model that we will be working on
	MGraph workModel;
	Vector2D lastDX;

public:
	void init();
	void press(MouseEvent &ev);
	void release(MouseEvent &ev);
	void drag(MouseEvent &ev);
	void display();

	TranslateToolHandler(DisplayWindow &win);

	// Get the label name for this tool
	char *getLabel() {
		return "Translation Tool";
	}

	// Get the context help for the tool
	char *getContextHelp() {
		return "Drag selection to move it around.";
	}
};


/******************************************************************
* Scaling and Rotation Tool Handler                              *
******************************************************************/
class ScaleRotateToolHandler : public ToolHandler {
	bool scaleMode;
	MouseEvent lastEvent;
	Vector2D center;

	double angle,scale,lastAngle,lastScale;

	// The model that we will be working on
	MGraph workModel;

	void computeAngleAndScale(MouseEvent &e);
public:
	void init();
	void press(MouseEvent &ev);
	void release(MouseEvent &ev);
	void drag(MouseEvent &ev);
	void display();

	ScaleRotateToolHandler(DisplayWindow &win,bool scaleMode);

	// Get the label name for this tool
	const char *getLabel() {
		return scaleMode ? "Scaling Tool" : "Rotation Tool";
	}

	// Get the context help for the tool
	const char *getContextHelp() {
		return scaleMode 
			? "Drag the mouse button to scale selection.  Click anywhere in the window to position the center of"
			" scaling.  Hold SHIFT key to allow simultaneous rotation and scaling."
			: "Drag the mouse button to rotate selection.  Click anywhere in the window to position the center of"
			" rotation.  Hold SHIFT key to allow simultaneous rotation and scaling.";
	}
};


/******************************************************************
* Stretch/Fatten Tool Handler                                    *
******************************************************************/
class StretchToolHandler : public ToolHandler {
public:
	void init();
	void release(MouseEvent &ev);
	void drag(MouseEvent &ev);
	void display();

	StretchToolHandler(DisplayWindow &win) : ToolHandler(win) {
	}

	// Get the label name for this tool
	char *getLabel() {
		return "Stretch Tool";
	}

	// Get the context help for the tool
	char *getContextHelp() {
		return "Use the window to stretch the model.";
	}
};


/******************************************************************
* Primitive Edit Mode Tool                                       *
* ------------------------                                       *
* In this mode the user edits one primitive at a time.           *
* Any primitive can be edited, regardless of selection.          *
******************************************************************/
class PEMToolHandler : public ToolHandler {
private:
	MGraph workModel;
	BRef target;
	MNode *node;

	bool alt,ctrl;

public:
	void press(MouseEvent &ev);
	void release(MouseEvent &ev);
	void drag(MouseEvent &ev);
	void display();
	void init();
	void close();

	PEMToolHandler(DisplayWindow &win);

	// Get the label name for this tool
	char *getLabel() {
		return "Atom Editor Tool";
	}

	// Get the context help for the tool
	char *getContextHelp() {
		return 
			"Drag red center of a primitive to move it.  Drag green circle to rotate primitive.  "
			"Drag blue circle to change radius and object angle.  Hold down CTRL to keep the left (blue) boundary "
			"primitive in place.  Hold ALT to hold the right (purple) boundary primitive in place";
	}

};


/*****************************************************************
 Path Recorder - Shared between tools
 ***************************************************************/
class PathRecorder {
protected:
	vector<BendPathPoint> path;
	void record(MouseEvent &ev);
	void drawPath();
public:
};

/******************************************************************
* Bend Tool Handler                                              *
******************************************************************/
class BendToolHandler : public ToolHandler, public PathRecorder {
private:
	// vector<BendPathPoint> path;
	MFigure *figure;
	// void record(MouseEvent &ev);

public:
	void init();
	void press(MouseEvent &ev);
	void release(MouseEvent &ev);
	void drag(MouseEvent &ev);
	void display();

	BendToolHandler(DisplayWindow &win) : ToolHandler(win) {
		this->figure = NULL;
	}

	// Get the label name for this tool
	char *getLabel() {
		return "Figure Bending Tool";
	}

	// Get the context help for the tool
	char *getContextHelp() {
		return 
			"Drag the mouse to paint a new medial track.  Model will be fitted to it.";
	}
};


/******************************************************************
* M-Spline Tool Handler                                          *
******************************************************************/
class SplineObject;

class SplineToolHandler : public ToolHandler, public PathRecorder {
private:
	// The curve that has been picked
	int pickedCurve,pickedNode;

	// The current mode of operation
	enum Modes {NOTHING=0,CHANGE_X,CHANGE_R,DRAW_SUB,DRAW_CURVE,JOIN_CURVES} mode;

	// Find point under cursor
	void findPoint(MouseEvent &ev);

	// Starting r value
	float rStart;

	friend void uiCreateSplineCallback(Fl_Return_Button*, void*) ;
	void addSplineCurve();

public:
	void init();
	void press(MouseEvent &ev);
	void release(MouseEvent &ev);
	void drag(MouseEvent &ev);
	void display();
	

	SplineToolHandler(DisplayWindow &win);
	~SplineToolHandler();

	// Get the label name for this tool
	char *getLabel() {
		return "Spline Editing Tool";
	}

	// Get the context help for the tool
	char *getContextHelp() {
		return 
			"Move the control points around with the left mouse button.  Change the radius at control "
			"points with the right mouse button";
	}
};


/******************************************************************
* Branch Link Tool Handler                                       *
******************************************************************/
class BranchToolHandler : public ToolHandler {
private:
	BRef brHead,brTail;
	Vector2D dragPoint;
	bool flipped;

public:
	void init();
	void close();
	void press(MouseEvent &ev);
	void release(MouseEvent &ev);
	void drag(MouseEvent &ev);
	void display();

	BranchToolHandler(DisplayWindow &win) : ToolHandler(win) {
	}

	// Get the label name for this tool
	const char *getLabel() {
		return "Branching Tool";
	}

	// Get the context help for the tool
	const char *getContextHelp() {
		return 
			"Click on the boundary atom from which to draw the link.  "
			"Drag the mouse to another atom to form a link, or off the "
			"screen to delete the existing link.  To make an indentation type "
			"link, hold down shift!";
	}

};


#endif //_TOOLS_H_
