#include "graph2d.h"
#include "g2DUI.h"

void Graph2D::createState(Vector xMin,Vector xMax,int res,Graph2DState &s) {
	// Set variables
	s.xMin = xMin;
	s.xMax = xMax;
	s.res = res;

	// Compute the graph
	ui->sldProgress->value(0);

	// Position the computing window
	ui->winComputing->position(ui->winMain->x()+ui->winGraph->x()+(ui->winGraph->w() - ui->winComputing->w())/2,
		ui->winMain->y()+ui->winGraph->y()+(ui->winGraph->h() - ui->winComputing->h())/2);
	ui->winComputing->show(); 



	s.F.setSize(res,res);

	// Position to evaluate
	Vector x = xMin;

	// Step sizes
	Vector step = (xMax-xMin) / (res-1);

	// Value xMax and xMin (always keep zero in the range)
	s.vMax = 0;s.vMin = 0;

	// Compute the function
	for(int i=0;i<res;i++) {
		for(int j=0;j<res;j++) {
			double fx = f.evaluate(x);

			// Update xMax and xMin
			if(fx > s.vMax)
				s.vMax = fx;
			else if(fx < s.vMin)
				s.vMin = fx;

			// Store in matrix
			s.F(i,j) = fx;

			// Update x
			x(1) += step(1);
		}

		// Update x
		x(1) = xMin(1);
		x(0) += step(0);

		// Update progress bar
		ui->sldProgress->value(i / (res-1.0));
		Fl::check();
	}

	// Hide progress bar
	ui->winComputing->hide();
}

void Graph2D::paintPixels() {
	// Current state reference
	Graph2DState &s = states[state];

	// Deallocate old pixel array
	if(pixels != NULL)
		delete pixels;

	// Allocate new pixel array
	pixels = new rgb[s.res*s.res];

	// Intensity multiplier
	double absMax = s.vMax > -s.vMin ? s.vMax : -s.vMin;   

	// Set pixel values
	int offset = 0;
	for(int j=0;j<s.res;j++) {
		for(int i=0;i<s.res;i++) {
			double v = s.F(i,j);
			v = 255*ui->winPalette->getIntensity(v/absMax);
			if(v < 0) {
				pixels[offset].b = (unsigned char) 255;
				pixels[offset].r = (unsigned char) (255 + v);
				pixels[offset].g = (unsigned char) (255 + v);
			}
			else {
				pixels[offset].r = (unsigned char) 255;
				pixels[offset].b = (unsigned char) (255 - v);
				pixels[offset].g = (unsigned char) (255 - v);
			}
			offset++;
		}
	}

	ui->winGraph->redraw();
}

void Graph2D::redraw() {
	ui->winGraph->redraw();
}

void Graph2D::setActiveState(int st) {
	// No setting invalid states
	assert(st < states.size() && st >= 0);

	// Set the state
	state = st;

	// Current state reference
	Graph2DState &s = states[state];

	// Set range controls
	ui->inpXMin->value(s.xMin(0));
	ui->inpYMin->value(s.xMin(1));
	ui->inpXMax->value(s.xMax(0));
	ui->inpYMax->value(s.xMax(1));

	// Set resolution control
	int h = (int) (::log(s.res) / log(2.0));
	ui->inpRes->value(h);
	ui->outRes->value(1 << h);

	// Pixel multiplier
	double absMax = (s.vMax > -s.vMin) ? s.vMax : -s.vMin;

	// Enable the back and forth controls
	if(state > 0)
		ui->bBack->activate();
	else
		ui->bBack->deactivate();

	if(state < states.size()-1)
		ui->bNext->activate();
	else
		ui->bNext->deactivate();

	// Special mark is no longer visible
	ui->winGraph->markVisible = false;
	mark(false);

	// Draw the pixels
	paintPixels();
}

Graph2D::Graph2D(Function &function,Vector xMin,Vector xMax,int res) :
f(function)
{
	// Make an assertion that the vectors are 2D
	assert(xMin.size() == 2 && xMax.size() == 2);

	// Initial state is 0
	state = 0;

	// Add a state to the state array
	Graph2DState blankState;
	states.push_back(blankState);

	// No pixels yet
	pixels = NULL;

	// No callback yet
	markCallback = false;

	// Fet flags
	axesVisible = true;
	gridVisible = gridLabeled = false;
	dragMode = true;

	// Create a user interface
	ui = new Graph2DUI();
	ui->graph = this;
	ui->makeWindow();
	ui->winGraph->ui = ui;
	ui->winPalette->ui = ui;
	ui->winMain->show();
	ui->winGraph->show();
	ui->winPalette->show();

	// Create a state
	createState(xMin,xMax,res,states[state]);

	// Set current state as active
	setActiveState(state);
}

// Destructor
Graph2D::~Graph2D() {
	if(pixels)
		delete pixels;
	ui->winMain->hide();
	delete ui;
}

// Compute the graph based on given settings
void Graph2D::compute() {
	// Read xMin and xMax
	Vector xMin(2),xMax(2);
	xMin(0) = ui->inpXMin->value();
	xMin(1) = ui->inpYMin->value();
	xMax(0) = ui->inpXMax->value();
	xMax(1) = ui->inpYMax->value();

	// Read resolution
	int res = (int) ui->outRes->value();

	// Remove all states after the present state
	states.erase(states.begin() + state + 1,states.end());
	states.push_back(Graph2DState());
	createState(xMin,xMax,res,states[state+1]);
	setActiveState(state+1);
}

void Graph2D::zoom(double x0,double y0,double x1,double y1) {
	// Read xMin and xMax
	Vector xMin(2),xMax(2);
	xMin(0) = (x0 < x1) ? x0 : x1;
	xMin(1) = (y0 < y1) ? y0 : y1;
	xMax(0) = (x0 < x1) ? x1 : x0;
	xMax(1) = (y0 < y1) ? y1 : y0;

	// Read resolution
	int res = (int) ui->outRes->value();

	// Remove all states after the present state
	states.erase(states.begin() + state + 1,states.end());
	states.push_back(Graph2DState());
	createState(xMin,xMax,res,states[state+1]);
	setActiveState(state+1);
}

void Graph2D::zoomOut() {
	Vector xMin = 0.5*(3.0*states[state].xMin-states[state].xMax);
	Vector xMax = 0.5*(3.0*states[state].xMax-states[state].xMin);
	zoom(xMin(0),xMin(1),xMax(0),xMax(1));
}

void Graph2D::zoomIn() {
	Vector xMin = 0.25*(3.0*states[state].xMin+states[state].xMax);
	Vector xMax = 0.25*(3.0*states[state].xMax+states[state].xMin);
	zoom(xMin(0),xMin(1),xMax(0),xMax(1));
}


void Graph2D::showAxis(bool mode) {
	axesVisible = mode;
	ui->winGraph->redraw();
}

// Show grid
void Graph2D::showGrid(bool mode) {
	gridVisible = mode;
	ui->winGraph->redraw();
}

// Show grid labels
void Graph2D::showLabels(bool mode) {
	gridLabeled = mode;
	ui->winGraph->redraw();
}

inline void vpToF(double vx,double vy,double &x,double &y,Graph2DState &s) {
	x = s.xMin(0) + (s.xMax(0)-s.xMin(0))*(1.0+vx) / 2.0;
	y = s.xMin(1) + (s.xMax(1)-s.xMin(1))*(1.0-vy) / 2.0;
}

inline void fToVP(double x,double y,double &vx,double &vy,Graph2DState &s) {
	vx = 2.0 * (x - s.xMin(0)) / (s.xMax(0)-s.xMin(0)) - 1.0;
	vy = 2.0 * (y - s.xMin(1)) / (s.xMax(1)-s.xMin(1)) - 1.0;
}

void Graph2D::mark(bool visible,double x,double y)
{
	if(visible) {
		ui->outX->value(x);
		ui->outY->value(y);
		ui->outFxy->value(f.evaluate(Vector(2,x,y)));
	}
	else {
		ui->outX->value(0);
		ui->outY->value(0);
		ui->outFxy->value(0);
	}
	if(markCallback)
		markCallback(visible,x,y);
}
/***********************************************************
Window paint routine
***********************************************************/
void Graph2DWindow::draw() {
	if (!valid()) {      
		glViewport(0,0,w(),h());
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glPushMatrix();
		glMatrixMode(GL_MODELVIEW);
		glFlush();
	}    

	// Fill the window
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glFlush();

	// Draw the pixels
	if(ui && ui->graph && ui->graph->states.size()) {
		// Get current state
		Graph2DState &s = ui->graph->states[ui->graph->state];

		// Paint pixels if any
		if(ui->graph->pixels) {
			// Figure out the zoom factor
			glPixelZoom(1.0*w()/s.res,1.0*h()/s.res);
			glDrawPixels(s.res,s.res,GL_RGB,GL_UNSIGNED_BYTE,ui->graph->pixels);
		}

		// Figure out the amount of grid shown
		double xGridStep = pow(10,floor(log10(s.xMax(0)-s.xMin(0))));
		double yGridStep = pow(10,floor(log10(s.xMax(1)-s.xMin(1))));
		double xGridStart = s.xMin(0) - fmod(s.xMin(0),xGridStep);
		double yGridStart = s.xMin(1) - fmod(s.xMin(1),yGridStep);

		// Convert to viewport coordinates
		fToVP(xGridStart,yGridStart,xGridStart,yGridStart,s);
		xGridStep *= 2.0 / (s.xMax(0)-s.xMin(0));
		yGridStep *= 2.0 / (s.xMax(1)-s.xMin(1));


		// Paint grid if enabled
		if(ui->graph->gridVisible) {
			glColor3d(0.0,0.75,0.0);
			glBegin(GL_LINES);

			for(double gx=xGridStart;gx<=1.0;gx+=xGridStep) {
				glVertex2d(gx,-1);glVertex2d(gx,1);
			}
			for(double gy=yGridStart;gy<=1.0;gy+=yGridStep) {
				glVertex2d(-1,gy);glVertex2d(1,gy);
			}
			glEnd();
		}

		// Paint axes if enabled
		if(ui->graph->axesVisible) {
			double axesX,axesY;
			fToVP(0,0,axesX,axesY,s);

			glColor3d(0.5,1.0,0.5);
			glBegin(GL_LINES);

			// Paint the axis themselves
			glVertex2d(axesX,-1);glVertex2d(axesX,1);
			glVertex2d(-1,axesY);glVertex2d(1,axesY);

			// Paint tickmarks
			for(double gx=xGridStart;gx<=1.0;gx+=xGridStep) {
				glVertex2d(gx,axesY-0.01);glVertex2d(gx,axesY+0.01);
			}
			for(double gy=yGridStart;gy<=1.0;gy+=yGridStep) {
				glVertex2d(axesX-0.01,gy);glVertex2d(axesX+0.01,gy);
			}

			glEnd();
		}
		if(ui->graph->gridLabeled) {

		}        

		// If dragging mode, draw a box for that
		if(ui->graph->dragMode && dragging) {
			glColor3d(1,1,0);
			glBegin(GL_LINE_LOOP);
			glVertex2d(dragX0,dragY0);
			glVertex2d(dragX0,dragY1);
			glVertex2d(dragX1,dragY1);
			glVertex2d(dragX1,dragY0);
			glEnd();
		}

		// If a special value is clicked 
		if(markVisible) {
			glColor3d(0,0,0);
			glBegin(GL_LINES);
			// drawCircle(markX,markY,0.02*sqrt(2.0));
			glVertex2d(markX-0.02,markY-0.02);
			glVertex2d(markX+0.02,markY+0.02);
			glVertex2d(markX-0.02,markY+0.02);
			glVertex2d(markX+0.02,markY-0.02);
			glEnd();
		}

		// Call the client draw method
		glPushMatrix();
		//glTranslated(s.xMin(0),s.xMin(1),0);
		//glScaled(1.0 / (s.xMax(0)-s.xMin(0)),1.0 / (s.xMax(1)-s.xMin(1)),1);
		glTranslated(-1-2*s.xMin(0),-1-2*s.xMin(1),0);
		glScaled(2 / (s.xMax(0)-s.xMin(0)),2 / (s.xMax(0)-s.xMin(0)),1);
		ui->graph->drawChild();
		glPopMatrix();

	}


	// Flush this crap
	glFlush();
}

/***********************************************************
Window handle routine
***********************************************************/
int Graph2DWindow::handle(int event) {    
	// If not initialized, get out of here
	if(!ui || !ui->graph || ui->graph->states.size() == 0)
		return 0;

	// Event coordinates in viewport coords
	double ex = 2.0 * Fl::event_x() / w() - 1.0;
	//double ey = 2.0 * Fl::event_y() / h() - 1.0;
	double ey = 1.0 - 2.0 * Fl::event_y() / h();

	// Quick state access
	Graph2DState &s = ui->graph->states[ui->graph->state];

	switch(event) {    
   case FL_PUSH:
	   if(ui->graph->dragMode) {
		   dragging = true;
		   dragX0 = ex;dragY0 = ey;
		   dragX1 = ex;dragY1 = ey;
	   }
	   return 1;

   case FL_DRAG:      
	   if(ui->graph->dragMode) {
		   dragX1 = ex;dragY1 = ey;
		   redraw();
	   }
	   return 1;

   case FL_RELEASE:          
	   if(ui->graph->dragMode && fabs(ex-dragX0) > 0.03 && fabs(ey-dragY0) > 0.03) {
		   dragging = false;
		   dragX1 = ex;dragY1 = ey;

		   // Compute actual zoom coordinates
		   double x0,y0,x1,y1;
		   vpToF(dragX0,-dragY0,x0,y0,s);
		   vpToF(dragX1,-dragY1,x1,y1,s);

		   // Zoom in!
		   ui->graph->zoom(x0,y0,x1,y1);
	   }
	   else if(Fl::event_key() == FL_Button+1) {
		   markVisible = true;
		   markX = ex;markY = ey;
		   double x,y;
		   vpToF(markX,-markY,x,y,s);

		   ui->graph->mark(true,x,y);
		   redraw();
	   }
	   else {
		   markVisible = false;
		   ui->graph->mark(false);
		   redraw();
	   }
	   // ... mouse up event ...      
	   return 1;
   case FL_KEYBOARD:
	   // ... keypress, key is in Fl::event_key(), ascii in Fl::event_text()
	   return 1;    
   default:
	   // tell fltk that I don't understand other events      
	   return 0;    
	}  
}


/***********************************************************
Palette Window paint routine
***********************************************************/
void Graph2DPalette::draw() {
	if (!valid()) {      
		glViewport(0,0,w(),h());
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glPushMatrix();
		glMatrixMode(GL_MODELVIEW);
		glFlush();
	}    

	// Fill the window
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glFlush();

	glDisable(GL_CULL_FACE);
	glShadeModel(GL_SMOOTH);

	// Draw the pixel ramp

	glBegin(GL_TRIANGLES);
	glNormal3d(0,0,1);

	glColor3d(0,0,1);
	glVertex2d(-1,-1);
	glColor3d(1+M(0,0),1+M(0,0),1.0);
	glVertex2d(M(0,0),-1);
	glVertex2d(M(0,0),M(0,1)*2.0 - 1.0);

	glVertex2d(M(0,0),M(0,1)*2.0 - 1.0);
	glColor3d(1,1,1);
	glVertex2d(0,1);
	glVertex2d(0,M(0,1)*2.0 - 1.0);
	glEnd();

	glBegin(GL_QUADS);
	glColor3d(1+M(0,0),1+M(0,0),1.0);
	glVertex2d(M(0,0),-1);
	glVertex2d(M(0,0),M(0,1)*2.0 - 1.0);
	glColor3d(1,1,1);
	glVertex2d(0,M(0,1)*2.0 - 1.0);
	glVertex2d(0,-1.0);
	glEnd();

	// Now the right side
	glBegin(GL_TRIANGLES);
	glColor3d(1,0,0);
	glVertex2d(1,-1);
	glColor3d(1.0,1-M(1,0),1-M(1,0));
	glVertex2d(M(1,0),-1);
	glVertex2d(M(1,0),M(1,1)*2.0 - 1.0);

	glVertex2d(M(1,0),M(1,1)*2.0 - 1.0);
	glColor3d(1,1,1);
	glVertex2d(0,1);
	glVertex2d(0,M(1,1)*2.0 - 1.0);
	glEnd();

	glBegin(GL_QUADS);
	glColor3d(1.0,1-M(1,0),1-M(1,0));
	glVertex2d(M(1,0),-1);
	glVertex2d(M(1,0),M(1,1)*2.0 - 1.0);
	glColor3d(1,1,1);
	glVertex2d(0,M(1,1)*2.0 - 1.0);
	glVertex2d(0,-1.0);
	glEnd();

	// Flush this crap
	glFlush();
}

/***********************************************************
Palette Window handle routine
***********************************************************/
int Graph2DPalette::handle(int event) {
	// If not initialized, get out of here
	if(!ui || !ui->graph)
		return 0;

	// Event coordinates in viewport coords
	double ex = 2.0 * Fl::event_x() / w() - 1.0;
	double ey = 1.0 - 2.0 * Fl::event_y() / h();

	switch(event) {    
   case FL_PUSH:
	   return 1;

   case FL_DRAG:      
	   return 1;

   case FL_RELEASE:          
	   if(ex < 0) {
		   M(0,0) = ex;
		   M(0,1) = (ey+1) / 2;
	   }
	   else {
		   M(1,0) = ex;
		   M(1,1) = (ey+1) / 2;
	   }
	   redraw();
	   ui->graph->repaint();

	   return 1;
   case FL_KEYBOARD:
	   // ... keypress, key is in Fl::event_key(), ascii in Fl::event_text()
	   return 1;    
   default:
	   // tell fltk that I don't understand other events      
	   return 0;    
	}  
}

double Graph2DPalette::getIntensity(double x) {
	double z = 0;
	if(x < 0) {
		if(x < M(0,0)) {
			z  = (x+1)*M(0,1)/(1+M(0,0));
		}
		else {
			z = M(0,1) + (x-M(0,0))*(1-M(0,1))/(-M(0,0));
		}
		return z-1;
	}
	else {
		if(x > M(1,0)) {
			z = M(1,1) + (x - M(1,0))*(-M(1,1))/(1-M(1,0));
		}
		else {
			z = 1 + x*(M(1,1)-1)/M(1,0);
		}
		return 1-z;
	}
}

void Graph2DPalette::reset() {
	M(0,0) = -1;
	M(0,1) = 0;
	M(1,0) = 1;
	M(1,1) = 0;
	ui->graph->repaint();
	redraw();
}
