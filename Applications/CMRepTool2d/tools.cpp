#include "misc.h"
#include "tools.h"
#include "dslwin.h"
#include "dsltool.h"
#include "undo.h"
#include <Fl/gl.h>
#include <gl/glu.h>

#include "ui.h"
#include "SplineMrep.h"

#include <minmax.h>

extern SplineObject *spline;

void onModelOrSelectionChange();

/*************************************************************
Zoom tool
**************************************************************/
ZoomPanToolHandler::ZoomPanToolHandler(DisplayWindow &win,bool zoom) : ToolHandler(win),zoomMode(zoom) 
{
}

void ZoomPanToolHandler::press(MouseEvent &ev) {
	// Store the original zoom factor and the original pan amount 
	initWinSize = dw.winSize;
	initWinPosition = dw.winPosition;
}

void ZoomPanToolHandler::drag(MouseEvent &e) {
	// Displacement
	Vector2D dx = e.xScreen - e.sScreen;

	// Performing a zoom-in operation
	if((zoomMode && e.button==1) || (!zoomMode && e.button!=1)) {
		// Adjust the window size
		dw.winSize = initWinSize * pow(2.0,4.0*dx.y);

		// Adjust the window position
		dw.winPosition = e.sSpace - (e.sSpace - initWinPosition)*(dw.winSize/initWinSize);
	}

	// Performing a panning operation
	else {
		// Adjust the window position
		dw.winPosition = initWinPosition-dx*dw.winSize;
	}

	dw.redraw();
}

void ZoomPanToolHandler::release(MouseEvent &e) {
}

void ZoomPanToolHandler::close() {
}


/*************************************************************
Selection tool
***********************************************************/
void SelectionToolHandler::press(MouseEvent &e) {
	// Record the event
	lastEvent = e;

	// Crear the selection in the figure
	if(!e.shift && !e.ctrl)
		dw.haveSelection = false;
}



void SelectionToolHandler::drag(MouseEvent &e) {
	// Record the event
	lastEvent = e;

	// Force redisplay
	dw.redraw();
}

void SelectionToolHandler::processSelectionClick(MouseEvent &e) {	

	// Find who we clicked on
	MNode *target = NULL;
	for(MNode *node=model.first();node!=NULL;node=model.next()) {
		double r = computeRadius(*node);			
		if(e.sSpace.distanceTo(node->x()) < r) {
			target = node;
			break;
		}
	}

	// If we clicked off the model, get out, deselect it if we didn't hole shift or
	// control
	if(target == NULL) {
		if(!e.shift && !e.ctrl) 
			model.deselect();
		return;
	}

	MFigure *figure = node->figure();
	MShape *shape = figure->shape();

	// If the node we clicked on has been selected already, select the whole figure
	// If the figure is selected, select the shape.  If the shape is selected, select the
	// whole graph.  If the graph is selected, select just the node
	if(!e.ctrl) {
		if(target->selected()) {
			if(figure->nSelected() == figure->size()) {
				if(shape->nSelected() == shape->size()) {
					if(model.nSelected() == model.nodes().size()) {
						// Deselect all, and select just our node
						if(!e.shift)
							model.deselect();
						target->select();
					}
					else {
						// Select the entire model
						model.select();
					}
				}
				else {
					if(!e.shift)
						model.deselect();
					shape->select();
				}
			}
			else {
				if(!e.shift)
					model.deselect();
				figure->select();
			}
		}
		else {
			// Deselect all
			if(!e.shift)
				model.deselect();
			target->select();
		}
	}
	else if(e.ctrl) {
		target->toggle();
	}
}

// #undef min
// #undef max

void SelectionToolHandler::processSelectionDrag(MouseEvent &e) {

	// Adjust the selection box
	Vector2D lt = Vector2D::minimal(e.sSpace,e.xSpace);
	Vector2D br = Vector2D::maximal(e.sSpace,e.xSpace);

	// If we are in plain mode, deselect all
	if(!e.shift && !e.ctrl) {
		model.deselect();
	}

	// Get a group of the guys in this selection box
	MNodeGroup ng = model.ggInBox(lt,br);

	// Toggle if in control mode, select if in regular mode
	if(e.ctrl)
		ng.toggle();
	else
		ng.select();
}

void SelectionToolHandler::release(MouseEvent &e) {
	// Record the event
	lastEvent = e;

	// Select all the primitives inside the box
	if(e.click)
		processSelectionClick(e);
	else
		processSelectionDrag(e);

	// Update the display list
	onModelOrSelectionChange();

	// Force redisplay
	dw.redraw();
}

void SelectionToolHandler::display() {
	if(dw.dragging) {
		glEnable(GL_LINE_STIPPLE);
		glLineStipple(3,0x5555);

		glColor3d(1.0,1.0,0.0);
		glBegin(GL_LINE_LOOP);
		glVertex2d(lastEvent.sSpace.x,lastEvent.sSpace.y);
		glVertex2d(lastEvent.sSpace.x,lastEvent.xSpace.y);
		glVertex2d(lastEvent.xSpace.x,lastEvent.xSpace.y);
		glVertex2d(lastEvent.xSpace.x,lastEvent.sSpace.y);
		glEnd();

		glDisable(GL_LINE_STIPPLE);
	}
}

/*************************************************************
Rotate tool
***********************************************************/
ScaleRotateToolHandler::ScaleRotateToolHandler(DisplayWindow &win,bool scaleMode) : ToolHandler(win) {
	this->scaleMode = scaleMode;
	center = model.selection()->getCOG();
	scale = 1.0;
	angle = 0.0;
}

void ScaleRotateToolHandler::computeAngleAndScale(MouseEvent &e) {
	Vector2D dx1 = e.sSpace-center;
	Vector2D dx2 = e.xSpace-center;

	double a = dx2.getSlopeDeg() - dx1.getSlopeDeg();
	double s = dx2.twoNorm() / dx1.twoNorm();

	angle = (!scaleMode || e.shift) ? a : 0.0;
	scale = (scaleMode || e.shift) ? s : 1.0;
}

void ScaleRotateToolHandler::press(MouseEvent &e) {
	lastEvent = e;

	// Create a copy of the model on which we will be operating
	workModel = model;

	lastScale = 1.0;
	lastAngle = 0.0;
}

void ScaleRotateToolHandler::drag(MouseEvent &e) {

	lastEvent = e;
	computeAngleAndScale(e);

	// Rotate the working model
	workModel.selection()->scaleBy(scale/lastScale,center);
	workModel.selection()->rotateBy(angle-lastAngle,center);

	// Save the old angles
	lastScale = scale;
	lastAngle = angle;

	// Update the editors ???
	updatePrimitiveEditor(&workModel);
	dw.redraw();
}

void ScaleRotateToolHandler::release(MouseEvent &e) {
	lastEvent = e;
	if(e.click) {
		center = e.xSpace;
	}
	else {
		computeAngleAndScale(e);

		// Save model on undo stack before making changes
		pushModel();

		model.selection()->scaleBy(scale,center);
		model.selection()->rotateBy(angle,center);

		// Model is updated
		onModelOrSelectionChange();
	}

	dw.redraw();
}

void ScaleRotateToolHandler::display() {
	// Draw a target at the center
	glColor4d(0,0,0,0.5);
	drawCircle(center,0.015*dw.winSize);
	glColor4d(1.0,0.3,0.3,0.5);
	drawCircle(center,0.013*dw.winSize);
	glColor4d(1,1,1,0.5);
	drawCircle(center,0.009*dw.winSize);
	glColor4d(1.0,0.3,0.3,0.5);
	drawCircle(center,0.006*dw.winSize);

	// When dragging, draw the selected object
	if(dw.dragging) {
		// Compute the rotate-into point
		Vector2D dx = (lastEvent.sSpace - center)*scale;
		Transform2D rtn = Transform2D::rotation(-M_PI*angle/180);
		Vector2D tip = center + rtn * dx;

		// Quickly display the working model
		draw(&workModel,true,false);

		glEnable(GL_LINE_STIPPLE);
		glLineStipple(3,0x5555);
		glColor4d(0.0,1.0,0.0,0.5);
		glBegin(GL_LINE_STRIP);
		glVertex2d(lastEvent.sSpace);
		glVertex2d(center);
		glVertex2d(tip);
		glEnd();
		glDisable(GL_LINE_STIPPLE);
		drawOutlineCircle(lastEvent.sSpace,0.01*dw.winSize,Gl_Color::rgb(0,1,0,0.5));
		drawOutlineCircle(tip,0.01*dw.winSize,Gl_Color::rgb(0,1,0));
	}
}

/*************************************************************
Translate tool
***********************************************************/
TranslateToolHandler::TranslateToolHandler(DisplayWindow &win) : ToolHandler(win) {   
}

void TranslateToolHandler::press(MouseEvent &e) {
	// Create a copy of the model on which we will be operating
	workModel = model;

	// Last displacement (need that so we don't reinitialize the work model each drag moment
	lastDX = Vector2D(0,0,0);
}

void TranslateToolHandler::drag(MouseEvent &e) {
	Vector2D dx = e.xSpace - e.sSpace;   

	workModel.selection()->translateBy(dx-lastDX);
	lastDX = dx;

	dw.redraw();
	updatePrimitiveEditor(&workModel);
}

void TranslateToolHandler::release(MouseEvent &e) {
	Vector2D dx = e.xSpace - e.sSpace;

	// Save model on undo stack before making changes
	pushModel();

	model.selection()->translateBy(dx);

	onModelOrSelectionChange();
	dw.redraw();
}

void TranslateToolHandler::display() {
	if(dw.dragging) {
		// Quickly display the working model
		draw(&workModel,true,false);
	}
}


/*************************************************************
Stretch tool
***********************************************************/
void StretchToolHandler::drag(MouseEvent &e) {
}

void StretchToolHandler::release(MouseEvent &e) {
}

void StretchToolHandler::display() {
}

/*************************************************************
Primitive Edit Mode tool
***********************************************************/
PEMToolHandler::PEMToolHandler(DisplayWindow &win) : ToolHandler(win) {
}

void PEMToolHandler::close() {
	onModelOrSelectionChange();
	dw.redraw();
}

void PEMToolHandler::init() {

	model.deselect();

	// Mode change
	eDisplayMode = DISPLAY_PEM;

	dw.haveSelection = false;

	onModelOrSelectionChange();
	dw.redraw();
}

void PEMToolHandler::press(MouseEvent &ev) {
	// Find what it was we clicked on
	locateBAtom(ev,target);

	// If we clicked on something good, select it
	model.deselect();
	if(target.node != NULL) {
		target.node->select();

		// Create a working copy of the model
		workModel = model;

		// Get a reference to the matching node
		node = workModel.selection()->node(0);
	}

	// Save the event
	alt = ev.alt;
	ctrl = ev.ctrl;
}

void PEMToolHandler::release(MouseEvent &ev) {
	// Do nothing if nothing is selected
	if(target.node==NULL)
		return;

	// Save for undo purposes
	pushModel();

	// Get data from work node to target node
	target.node->setData(node->getData(true),true);

	// Select the node 
	model.deselect();
	target.node->select();

	// Update selection box and model display
	onModelOrSelectionChange();

	// Redraw
	dw.redraw();
}

void PEMToolHandler::drag(MouseEvent &ev) {
	// Do nothing if nothing is selected
	if(target.node==NULL)
		return;

	node->setData(target.node->getData(true),true);

	// Shorthand
	int bidx = target.bIndex;

	// Displacement
	Vector2D m = node->x();
	Vector2D br = node->x(MAtom::RIGHT);
	Vector2D bl = node->x(MAtom::LEFT);
	Vector2D dx = ev.xSpace - ev.sSpace;

	// Apply change to the atom
	if(bidx == MAtom::MIDDLE) {
		// Atom is displaced
		node->translateBy(dx);

		// Primitive is displaced in such a way that the position of the boundary primitive remains the same,
		// but the radius and object angle change
		if(ctrl) {
			Vector2D blink = br - (m + dx);
			node->r(blink.twoNorm());
			node->fa(-blink.getSlope() - node->oa());
		}
		else if(alt) {
			Vector2D blink = bl-(m+dx);
			node->r(blink.twoNorm());
			node->fa(-blink.getSlope() + node->oa());
		}
	}

	else if(bidx == MAtom::HEAD) {
		// Primitive is rotated
		Vector2D dx1 = ev.sSpace - m;
		Vector2D dx2 = ev.xSpace - m;
		double a = dx1.getSlopeDeg() - dx2.getSlopeDeg();
		node->rotateBy(a);

		// Rotate keeping a primitive in place
		if(ctrl) {
			node->changeAngleBy(a);
		}
		else if(alt) {
			node->changeAngleBy(-a);
		}      
	}

	else if(bidx == MAtom::LEFT || bidx == MAtom::RIGHT) {
		if(ctrl || alt) {
			// The vertices of a triangle are the static point, the moving point and the new center
			Vector2D xMoving = ev.xSpace;
			Vector2D xStatic = node->x(bidx == MAtom::LEFT ? MAtom::RIGHT : MAtom::LEFT);

			// The normal vector
			Vector2D nStatic = xStatic-m;
			nStatic.normalize();

			// The midpoint of the segment
			Vector2D xMidpoint = (xStatic+xMoving) / 2;

			// Compute the center -<d,d>/2<d,n>
			Vector2D d = xStatic - xMoving;
			double r = -0.5 * d.dotProduct(d) / d.dotProduct(nStatic);
			Vector2D xCenter = xStatic + nStatic*r;

			// OK, now we have a triangle.  Make it into a primitive
			double aAxial = - (xMidpoint-xCenter).getSlope();
			double radius = fabs(r);

			// Compute object angle.
			double chordLength = d.twoNorm();
			double aObject = asin(chordLength / (2*radius));

			// It is possible that we need to flip the axial angle towards the next primitive
			if(cos(aAxial-node->fa()) < 0) {
				aAxial += M_PI;
				aObject = M_PI-aObject;
			}

			node->x(xCenter);
			node->r(radius);
			node->fa(aAxial);
			node->oa(aObject);
		}
		else {
			// Primitive is rotated
			Vector2D dx1 = ev.sSpace - m;
			Vector2D dx2 = ev.xSpace - m;
			double a = dx2.getSlopeDeg() - dx1.getSlopeDeg();
			double s = dx2.twoNorm() / dx1.twoNorm();

			node->changeAngleBy((bidx == MAtom::LEFT) ? -a : a);
			node->scaleBy(s);
		}
	}

	// Update the editor
	updatePrimitiveEditor(node);

	// Redraw
	dw.redraw();
}

void PEMToolHandler::display() {
	if(dw.dragging) {
		// Quickly display the working model
		draw(&workModel,true,true);

		// Draw a circle for boundary primitive in place
		if(alt) {
			glColor3d(1,1,1);
			drawCircle(node->x(MAtom::LEFT),0.4*computeRadius(*node));
		}
		if(ctrl) {
			glColor3d(1,1,1);
			drawCircle(node->x(MAtom::RIGHT),0.4*computeRadius(*node));
		}
	}   
}

/***********************************************************
Spline Tool
***********************************************************/
void fitSplineToPoints(SplineCurve &C,vector<BendPathPoint> &Q) {
	int k;

	// M is the index of the last point Q	
	int m = (int)Q.size()-1;

	// N is the index of the last control point
	int n = C.size()-1;

	// D is the length of the curve Q
	double d = 0;
	for(k=1;k<=m;k++) 
		d += Q[k].x.distanceTo(Q[k-1].x);
	
	// Construct an array of u-values for the Q points
	vector<float> uk;
	uk.push_back(0.0f);
	for(k=1;k<m;k++) {
		uk.push_back(uk.back() + Q[k].x.distanceTo(Q[k-1].x) / d);
	}
	uk.push_back(1.0f);

	// Find knot indices of each uk
	vector<int> uki;
	C.kv.getKnotIndexSequence(uk,uki);
	
	// Constuct a 'band' matrix of N_i,p(u_j)
	Matrix BN(m+1,n+1);
	for(int r=0;r<=m;r++) {
		MySMLVec4f W[4];
		C.kv.basisJet(uki[r],uk[r],W);
		int c0 = uki[r]-3;
		int c1 = min(c0+3,n);
		for(int c=c0;c<=c1;c++) {			
			BN(r,c) = W[0].data()[c-c0];
		}
	}

	// Construct matrices N,R
	Matrix R(n-1,2);
	Matrix N(m-1,n-1);

	// R is computed as a sum
	for(int c=1;c<=n-1;c++) {
		for(int r=1;r<=m-1;r++) {
			Vector2D Rk = Q[r].x-Q[0].x*BN(r,0)-Q[m].x*BN(r,n);
			R(c-1,0) += Rk.x * BN(r,c);
			R(c-1,1) += Rk.y * BN(r,c);
		}
	}

	// N is a chunk of BN
	BN.extractMatrix(1,1,N);
	Matrix NTN = N.t() * N;

	// Solve for P
	NTN.ipSolveLinearSystem(R);

	// P containts the results
	C.updateControlX(0,Q[0].x);
	for(int i=1;i<=n-1;i++) {
		C.updateControlX(i,Vector2D(R(i-1,0),R(i-1,1)));
	}
	C.updateControlX(n,Q[m].x);
}

/**
 * Called when user closes the create spline dialog
 */
SplineToolHandler *uiCallbackSplineTool = NULL;
void uiCreateSplineCallback(Fl_Return_Button*, void*) {
	winNewSpline->hide();
	uiCallbackSplineTool->addSplineCurve();
}



void SplineToolHandler::findPoint(MouseEvent &ev) {
	pickedCurve = -1;
	pickedNode = -1;

	// Check which point was pressed
	for(unsigned int c=0;c<spline->curves.size();c++) {
		SplineCurve *curve = spline->curves[c];
		for(int i=0;i<curve->size();i++) {
			MySMLVec3f P = curve->getControlVec(i);
						
			Vector2D x = dw.spaceToScreen(P.x,P.y);
			if(x.distanceTo(ev.xRawScreen) < 6) {
				pickedCurve = c;
				pickedNode = i;
			}
		}
	}
}

void SplineToolHandler::press(MouseEvent &ev) {
	
	findPoint(ev);

	if(ev.shift) {
		// The user is drawing a branch
		mode = (pickedCurve >= 0) ? mode = DRAW_SUB : DRAW_CURVE;
		path.clear();
		record(ev);
	}
	else if(ev.ctrl && pickedCurve >= 0) {
		spline->removeCurve(pickedCurve);
		pickedCurve = -1;
		onModelOrSelectionChange();
		dw.redraw();
		mode = NOTHING;
	}
	else if(ev.alt && pickedCurve >= 0) {
		// Create a branch tool
		mode = JOIN_CURVES;
	}
	else if(pickedCurve >= 0) {
		rStart = spline->curves[pickedCurve]->getControlVec(pickedNode).z;
		mode = (ev.button==1) ? CHANGE_X : CHANGE_R;
	}
	else {
		mode = NOTHING;
	}
}

void SplineToolHandler::addSplineCurve() {
	int np = inNewSplineControls->value()-3;
	double rho = inNewSplineRho->value();
	double r = inNewSplineRadius->value();

	// Create a new curve
	SplineCurve *curve = new SplineCurve(np+3,rho,r);

	// Compute the curve by fitting to points
	fitSplineToPoints(*curve,path);

	// Add branch or curve
	if(mode == DRAW_CURVE) {
		spline->addCurve(curve);
	}
	else {
		spline->addBranch(spline->curves[pickedCurve],pickedNode,curve,0);
	}
	
	onModelOrSelectionChange();
	dw.redraw();
	mode = NOTHING;
}

void SplineToolHandler::release(MouseEvent &ev) {
	if(mode == DRAW_CURVE || mode == DRAW_SUB) {
		// Save the last point
		record(ev);

		// Prompt for the number of patches to use
		uiCallbackSplineTool = this;
		winNewSpline->show();
		return;
	}
	else if(mode == JOIN_CURVES) {
		int pc = pickedCurve;
		int pn = pickedNode;
		findPoint(ev);

		if(pickedCurve >= 0 && pickedNode >= 0) {
			spline->addBranch(spline->curves[pickedCurve],pickedNode,spline->curves[pc],pn);
			onModelOrSelectionChange();
		}
	}
	dw.redraw();
	mode = NOTHING;
}

void SplineToolHandler::drag(MouseEvent &ev) {
	if(mode == CHANGE_X) {
		spline->curves[pickedCurve]->updateControlX(pickedNode,ev.xSpace);
		onModelOrSelectionChange();
		dw.redraw();
		
	}
	else if(mode == CHANGE_R) {
		float dist = ev.xSpace.y - ev.sSpace.y;
		spline->curves[pickedCurve]->updateControlR(pickedNode,max(rStart + dist,0.0f));
		onModelOrSelectionChange();
		dw.redraw();
	}
	else if(mode == DRAW_CURVE || mode == DRAW_SUB) {
		record(ev);
		dw.redraw();
	}	
	
}


void SplineToolHandler::display() {
	if(mode == DRAW_CURVE || mode == DRAW_SUB) {
		drawPath();	
	}
	else if(mode == CHANGE_R) {		
		MySMLVec3f C = spline->curves[pickedCurve]->getControlVec(pickedNode);
		glColor3d(1,1,0);

		glPushAttrib(GL_LINE_BIT | GL_COLOR_BUFFER_BIT);
		glEnable(GL_BLEND);
		glEnable(GL_LINE_SMOOTH);

		glPushMatrix();
		glTranslated(C.x,C.y,0);

		glBegin(GL_LINE_LOOP);
		for(double a=0;a<M_PI*2;a+=M_PI / 90) {
			glVertex2d(C.z * cos(a),C.z * sin(a));
		}
		glEnd();

		glPopMatrix();

		glPopAttrib();

		// drawCircle(Vector2D(C.x,C.y),C.z);
	}
}

SplineToolHandler::SplineToolHandler(DisplayWindow &win) : ToolHandler(win) {

}

SplineToolHandler::~SplineToolHandler() {

}

/*************************************************************
Path Recorder
***********************************************************/
void PathRecorder::record(MouseEvent &e) {
	// Record the point's position
	if(path.size() == 0 || path.back().x.distanceTo(e.xSpace) > 0.000001) {
		BendPathPoint p;
		p.x = e.xSpace;
		path.push_back(p);
	}
}

void PathRecorder::drawPath() {
	glPushAttrib(GL_LINE_BIT);
	
	glLineWidth(3.0);
	glColor3d(1,0,0);

	glBegin(GL_LINE_STRIP);
	for(unsigned int i=0;i<path.size();i++) {
		glVertex2d(path[i].x);
	}
	glEnd(); 

	glPopAttrib();
}

/*************************************************************
Bend tool
***********************************************************/
/*
void BendToolHandler::record(MouseEvent &e) {
// Record the point's position
BendPathPoint p;
p.x = e.xSpace;
path.push_back(p);
}*/

void BendToolHandler::press(MouseEvent &e) {
	// Clear the path
	path.clear();

	// Get a figure
	figure = getSingleSelectedFigure();

	// Record the position
	record(e);
}

void BendToolHandler::drag(MouseEvent &e) {
	// Only applicable if a figure exists
	if(!figure)
		return;

	// Record the position
	record(e);

	// Request redisplay
	dw.redraw();
}

void BendToolHandler::release(MouseEvent &e) {
	// Only applicable if a figure exists
	if(!figure)
		return;

	// Record the position
	record(e);

	// Need a table for the weights of these things
	const double wTable[] = {0,20,15,6,1};

	// Compute distance along the path
	unsigned int i;
	for(i=0;i<path.size();i++) {
		// Set the distance parameter
		if(i==0)
			path[i].s = 0;
		else
			path[i].s = path[i-1].s + path[i].x.distanceTo(path[i-1].x);
	}

	// OK. Now the path is computed.  Lets see the total length of the path
	if(path.size() == 0 || path[path.size()-1].s < dw.winSize / 100.0) {
		// The path is too short
		path.clear();
		dw.redraw();
		return;
	}

	// Save the model
	pushModel();

	// Compute normal at each point in the path
	for(i=0;i<path.size();i++) {
		double tSum=0,wSum=0;

		// Compute the tagent
		for(int j=-4;j<=4;j++) {
			// Normal is not affected by the point itself
			if(j==0 || i+j<0 || i+j >= path.size()) 
				continue;

			// The contributing point
			int k = i+j;

			// The tangent vector estimated by this pair of points
			double t = (k<i) ? -(path[i].x-path[k].x).getSlopeDeg() : -(path[k].x-path[i].x).getSlopeDeg();

			// Assign a weight to this component
			double w = wTable[abs(j)];
			tSum+=t*w;wSum+=w;
		}

		path[i].t = tSum/wSum;
	}

	// First find the total length of the figure
	double fLen = 0,dist = 0,pLen = path[path.size()-1].s;
	Vector2D lastX;
	int pp = 0;

	for(int p=0;p<figure->size();p++) {
		if(p > 0)
			fLen += figure->node(p)->x().distanceTo(lastX);
		lastX = figure->node(p)->x();
	} 

	// Now, for each atom find the closest point on the path
	for(p=0;p<figure->size();p++) {
		MNode *node = figure->node(p);
		Vector2D x = node->x();

		// This is very simple for first and last primitives
		if(p==0) {
			node->x(path[0].x);
			node->faDeg(path[0].t);
		}
		else if(p==figure->size()-1) {
			node->x(path[path.size()-1].x);
			node->faDeg(path[path.size()-1].t);
		}
		else {
			// Calculate cumulative distance along model in units of length of the path
			dist += x.distanceTo(lastX) * pLen / fLen;

			// Find the first path point beyond that distance
			while(path[pp].s < dist)
				pp++;

			// If pp is zero, the model is screwed up and we will skip that primitive
			if(!pp)
				continue;

			// Position of the 'other' path point we'll interpolate with
			double ppLast = pp-1;
			double w1 = (dist-path[ppLast].s)/(path[pp].s-path[ppLast].s);
			double w2 = (path[pp].s-dist)/(path[pp].s-path[ppLast].s);


			// OK, now interpolate the position and the normal
			node->x(path[pp].x*w1+path[ppLast].x*w2);
			node->faDeg(path[pp].t*w1+path[ppLast].t*w2);
		}

		lastX = x;
	}

	// Clear the path - we don't need it any more
	path.clear();

	// Model has changed - call the method
	onModelOrSelectionChange();
	dw.redraw();
}

void BendToolHandler::display() {
	// OK, here we just draw the path
	if(dw.dragging) {
		drawPath();
	}
}

/*************************************************************
Branch tool
***********************************************************/
// Mouse is pressed at the link at which we are going to want to begin.  
void BranchToolHandler::press(MouseEvent &ev) {
	// Find out which atom it is
	locateBAtom(ev,brTail);         
	brHead.node = NULL;

	// Middles don't count
	if(brTail.bIndex == MAtom::MIDDLE)
		brTail.node = NULL;

	flipped = ev.shift;
}

// As we drag, we need to redisplay
void BranchToolHandler::drag(MouseEvent &ev) {
	if(brTail.node) {
		// Find an atom under the mouse
		locateBAtom(ev,brHead);

		// Middles don't count
		if(brHead.bIndex == MAtom::MIDDLE)
			brHead.node = NULL;

		// Record the position, snap to center
		dragPoint = (brHead.node) ? brHead.node->x(brHead.bIndex) : ev.xSpace;

		// Repaint
		dw.redraw();
	}
}

// As we realease, do the magic
void BranchToolHandler::release(MouseEvent &ev) {
	if(brTail.node) {
		// Find an atom under the mouse
		locateBAtom(ev,brHead);

		// Middles don't count
		if(brHead.bIndex == MAtom::MIDDLE)
			brHead.node = NULL;

		// A click (drag from node to itself removes the arrow as well
		if(brHead.node == brTail.node)
			brHead.node = NULL;


		if(brHead.node) {
			// Backup model
			pushModel();

			// Create a branch
			model.beginTopologyEdit();
			model.addBranch(brTail,brHead,flipped);
			model.endTopologyEdit();
		}

		else {
			// Remove the node if it exists
			model.beginTopologyEdit();
			model.removeBranch(brTail);
			model.endTopologyEdit();
		}
	}

	brTail.node = NULL;
	brHead.node = NULL;

	// Repaint
	dw.redraw();
}

// Display method
void BranchToolHandler::display() {   
	// Draw all the existing branches
	for(MLinkItr it = model.links().begin();it != model.links().end();it++) {
		MLink *link = *it;
		if(link->type() == MLink::BRANCH) {
			BranchLink *bl = (BranchLink*) link;

			if(bl->flipped()) {
				glColor3d(0,0,1);
			} else {
				glColor3d(0,1,0);
			}

			const BAtom *aTail = bl->tailBRef().bAtom();
			const BAtom *aHead = bl->headBRef().bAtom();

			glDrawArrow(aTail->x,aHead->x,0.02*dw.winSize);
		}
	}


	if(brTail.node) {
		// Color of the arrow changes if we pass over a link
		if(brHead.node) {
			glColor3d(1.0,1.0,0.0);
		} else {
			glColor3d(1.0,0.0,0.0);
		}

		// Display the link we are currently drawing...
		glDrawArrow(brTail.node->x(brTail.bIndex),dragPoint,0.02*dw.winSize);
	}
}

void BranchToolHandler::close() {
	dw.redraw();
}

void BranchToolHandler::init() {
	eDisplayMode = DISPLAY_MREP;
	dw.redraw();
}


void BendToolHandler::init()
{
	eDisplayMode = DISPLAY_MREP;
	dw.redraw();
}

void ScaleRotateToolHandler::init()
{
	eDisplayMode = DISPLAY_MREP;
	dw.redraw();
}

void SplineToolHandler::init()
{
	eDisplayMode = DISPLAY_SPLINE;
	dw.redraw();
}

void StretchToolHandler::init()
{
	eDisplayMode = DISPLAY_MREP;
	dw.redraw();

}

void TranslateToolHandler::init()
{
	eDisplayMode = DISPLAY_MREP;
	dw.redraw();
}

void SelectionToolHandler::init()
{
	eDisplayMode = DISPLAY_MREP;
	dw.redraw();
}
