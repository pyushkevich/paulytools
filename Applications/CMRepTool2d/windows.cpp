#include "dslwin.h"
#include "tools.h"
#include "dsltool.h"
#include "ispace.h"
#include "misc.h"
#include "Fl/gl.h"
#include "Fl/Fl.H"
#include "GL/glu.h"
#include "ui.h"

#define MAX_DOUBLE 1e100

#include <minmax.h>

#include "SplineMRep.h"

extern RegularSplineSample *splineSample;

extern void initGL();

/**********************************************************************
* Implementation for DisplayWindow
**********************************************************************/
DisplayWindow::DisplayWindow(int X, int Y, int W, int H, const char *L) 
: Fl_Gl_Window(X, Y, W, H, L) 
{
	// Set zoom/pan
	winPosition;
	winSize = 1.0;

	// We are not dragging
	dragging = false;

	// Set tool to default (do nothing) tool
	tool = new SelectionToolHandler(*this);

	// Deal with selection
	haveSelection = false;
	selTopLeft;
	selBotRight;

	glTextureName = 0; // dlModel = dlSelection = -1;
	glDLName = -1;

	needTextureUpdate = true;
}

// There is a core drawing procedure that works only when a core exists
extern void drawCore();

void DisplayWindow::draw() {
	if (!valid()) {
		initGL();

		// Set up the viewport
		glViewport(0,0,w(),h());
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		// gluPerspective(60.0, (GLdouble)w/h, 0.01, 100);
		// glPushMatrix();
		
		glMatrixMode(GL_MODELVIEW);
		// glLoadIdentity();
		// glFlush();

	}

	if(needTextureUpdate) {

		createTexture();
		needTextureUpdate = false;
		glDLName = -1;

	}

	// Fill the window
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glFlush();

	// Set up matrices for 2D display  
	glPushMatrix();

	// Multiply by a transform matrix that maps cartesian square (-1,-1):(1,1) to a image square (0,0):(1,1);
	glTranslated(-1,1,0);
	glScaled(2,-2,1);	

	// Draw the center of the viewport.  
	/* 
	glColor3d(0.4,0.4,0.4);
	glBegin(GL_LINES);
	glVertex2d(0.45,0.5);
	glVertex2d(0.55,0.5);
	glVertex2d(0.5,0.45);
	glVertex2d(0.5,0.55);
	glEnd();
	*/


	glPushMatrix();
	glScaled(1.0/winSize,1.0/winSize,1.0);
	glTranslated(-winPosition.x,-winPosition.y,0.0);

	glPushMatrix();
	glTranslated(0,0,1.0);

	// Display
	drawGrid();

	// Display the image polygon
	drawImage();

	glPopMatrix();

	// Display the model
	switch(eDisplayMode) {
	case DISPLAY_MREP : 
		::draw(&model,false,false);
		break;

	case DISPLAY_PEM : 
		::draw(&model,false,true);
		break;

	case DISPLAY_BEND :	
		break;

	case DISPLAY_BRANCH : 
		::drawBranchMode(&model);
		break;

	case DISPLAY_SPLINE : 
		drawSpline(splineSample);
		break;
	}

	// If selection is enabled, display the selection box
	if(haveSelection && eDisplayMode == DISPLAY_MREP) {
		double border = 0.01*winSize;

		glEnable(GL_LINE_STIPPLE);
		glLineStipple(3,0x5555);
		glColor3d(0.0,1.0,1.0);
		glBegin(GL_LINE_LOOP);
		glVertex2d(selTopLeft.x-border,selTopLeft.y-border);
		glVertex2d(selTopLeft.x-border,selBotRight.y+border);
		glVertex2d(selBotRight.x+border,selBotRight.y+border);
		glVertex2d(selBotRight.x+border,selTopLeft.y-border);
		glEnd();
		glDisable(GL_LINE_STIPPLE);
	}

	// Draw the core if any
	drawCore();

	// Call the tool's display code
	tool->display();

	// Pop the display matrix
	glPopMatrix();
	glPopMatrix();

	glFlush();
}

inline Vector2D DisplayWindow::screenToSpace(double x,double y) {
	return winPosition + Vector2D(x/w(),y/h()) * winSize;
}

Vector2D DisplayWindow::spaceToScreen(double x,double y) {
	Vector2D v1 = (Vector2D(x,y) - winPosition) / winSize;
	v1.x *= w();
	v1.y *= h();
	return v1;
}

int DisplayWindow::handle(int eid) {
	static MouseEvent *event = NULL;
	Vector2D xScreen((double)Fl::event_x()/w(),(double)Fl::event_y()/h());
	Vector2D xSpace = screenToSpace(Fl::event_x(),Fl::event_y());

	switch(eid) {
	case FL_PUSH:
		// A new event has occured.  Create a new event object
		if(event) 
			delete event;
		event = new MouseEvent();

		// Set the start settings for the event
		event->xRawScreen = event->sRawScreen = Vector2D(Fl::event_x(),Fl::event_y());
		event->xScreen = event->sScreen = xScreen;
		event->xSpace = event->sSpace = xSpace;
		event->alt = (Fl::event_state() & FL_ALT) != 0;
		event->ctrl = (Fl::event_state() & FL_CTRL) != 0;
		event->shift = (Fl::event_state() & FL_SHIFT) != 0;
		event->button = Fl::event_button();

		// Call the tool's response method
		tool->press(*event);

		return 1;

	case FL_DRAG:
		// Add to the event's drag length
		event->dragLength += xScreen.distanceTo(event->xScreen);

		// Update the current position
		event->xRawScreen = Vector2D(Fl::event_x(),Fl::event_y());
		event->xScreen = xScreen;
		event->xSpace = xSpace;

		// We are dragging
		dragging = true;	

		// Call the tool's response method
		tool->drag(*event);

		return 1;

	case FL_RELEASE:

		// Update the current position
		event->xRawScreen = Vector2D(Fl::event_x(),Fl::event_y());
		event->xScreen = xScreen;
		event->xSpace = xSpace;

		// We are not dragging
		dragging = false;

		// Determine if click or drag
		if(event->dragLength < 0.01 && event->xScreen.distanceTo(event->sScreen) < 0.01)
			event->click = true;

		// Call the tool's response method
		tool->release(*event);

		return 1;

	case FL_SHORTCUT:
	case FL_KEYBOARD:
		{
			char c = *Fl::event_text();
			if(c >= '1' && c <= '9')
				uiRunTool((ToolEnum)(c-'0'));
			else if(c==' ') {
				eToolId = (ToolEnum)((eToolId+1) % 11);
				uiRunTool(eToolId);
			}		
		}
	}

	return 0;
}

void DisplayWindow::setTool(ToolHandler *tool) {
	this->tool->close();

	// delete this->tool;

	this->tool = tool;
	this->tool->init();

	redraw();
}


void DisplayWindow::createTexture() {
	// Allocate a texture name
	if(glTextureName <= 0) {
		glGenTextures(1,&glTextureName);
	}

	// Get an image
	if(imageSpace==NULL)
		return;

	// Get a byte array from the image
	Array2D<unsigned char> bytes;
	imageSpace->getImageBytes(bytes);

	// Load the texture into graphics memory
	if(bytes.width() > 0 && bytes.height() > 0) {

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D,glTextureName);

		int interp = regOptions.cmpStringValue("display.gl.interpolation","linear") ? GL_LINEAR : GL_NEAREST;

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, interp);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, interp);

		//glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
		//glTexImage2D(GL_TEXTURE_2D,0,1,cimage_xdim(baseImage),cimage_ydim(baseImage),0,
		//	GL_LUMINANCE,GL_UNSIGNED_BYTE,baseBytes.getData());
		gluBuild2DMipmaps(GL_TEXTURE_2D,1,bytes.width(),bytes.height(),
			GL_LUMINANCE,GL_UNSIGNED_BYTE,bytes.getData());

		glDisable(GL_TEXTURE_2D);
	}

}

double myfmod(double d,double z) {
	return floor(d/z) * z;
}

void DisplayWindow::drawGrid() {
	if(imageSpace == NULL ||  imageSpace->width() <= 0) {
		double g = pow(10,ceil(log10(winSize)));
		double gLast = 0;
		double fact = 0.5;
		double clr = 1.0;
		int nLines = 0;

		vector<double> gArray;
		while(winSize / g < 40) {
			gArray.push_back(g);
			gLast = g;
			g = fact * g;			
			fact = 0.1 / fact;
			clr  *= 0.7;
		}

		glBegin(GL_LINES);

		for(int i=(int)gArray.size()-1;i >= 0;i--) {
			glColor3d(clr,clr,clr); 
			g = gArray[i];

			for(double t = myfmod(winPosition.x,g);t < winPosition.x + winSize;t += g) {
				// if(gLast == 0 || (t - myfmod(t,gLast)) > 0.1*g) {
				glVertex2d(t,winPosition.y);
				glVertex2d(t,winPosition.y+winSize);
				//}
			}

			for(t = myfmod(winPosition.y,g);t < winPosition.y + winSize;t += g) {
				//if(gLast == 0 || (t - myfmod(t,gLast)) > 0.1*g) {
				glVertex2d(winPosition.x,t);
				glVertex2d(winPosition.x+winSize,t);
				//}
			}

			clr /= 0.7;
		}

		glEnd();
	}
}

void DisplayWindow::drawImage() {
	if(imageSpace != NULL && imageSpace->width() > 0) {

		if(glDLName < 0) {
			glDLName = glGenLists(1);

			glNewList(glDLName,GL_COMPILE);

			glPushAttrib(GL_ALL_ATTRIB_BITS);
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D,glTextureName);

			glColor3d(1,1,1);
			glNormal3d(0,0,1);
			glBegin(GL_QUADS);
			glTexCoord2f(0.0,0.0);glVertex2d(0,0);
			glTexCoord2f(1.0,0.0);glVertex2d(1,0);
			glTexCoord2f(1.0,1.0);glVertex2d(1,1);
			glTexCoord2f(0.0,1.0);glVertex2d(0,1);
			glEnd();
			glPopAttrib();		

			glEndList();
		}

		


		/*
		glPushAttrib(GL_TEXTURE_BIT);
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D,glTextureName);

		glColor3f(1,1,1);
		glNormal3d(0,0,1);
		glBegin(GL_QUADS);
		glTexCoord2f(0.0,0.0);glVertex2d(0,0);
		glTexCoord2f(1.0,0.0);glVertex2d(1,0);
		glTexCoord2f(1.0,1.0);glVertex2d(1,1);
		glTexCoord2f(0.0,1.0);glVertex2d(0,1);
		glEnd();

		glPopAttrib();		*/
		glCallList(glDLName);
	}
}

/*************************************************************
Scale Space Window Controls
***********************************************************/
void ScaleSpaceWindow::draw() {
	if (!valid()) {
		// Set up the viewport
		glViewport(0,0,w(),h());
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		// Set up the image
		iw = w() - w()%4;
		ih = h() - h()%4;

		image.resize(iw,ih);
		setScale(scale);

		// gluPerspective(60.0, (GLdouble)w/h, 0.01, 100);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}

	if(imageSpace==NULL)
		return;

	// Fill the window
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glFlush();

	glDrawPixels(iw,ih,GL_LUMINANCE,GL_UNSIGNED_BYTE,image.getData());

	glFlush();
}

int ScaleSpaceWindow::handle(int) {
	return 0;
}

void ScaleSpaceWindow::setScale(double newScale) {
	scale = newScale;

	if(imageSpace==NULL)
		return;

	double imax = imageSpace->getImageMax();
	double imin = imageSpace->getImageMin();
	double iRange = imax-imin;

	Vector2D p;

	// Create a new image that could fit the blurred image
	Array2D<double> tmpImage((int)ceil(1.0*iw/zoom),(int)ceil(1.0*ih/zoom));

	// Compute the steps in which the position in the image will be changing
	double yStep = - ((double)imageSpace->height()) / tmpImage.height();
	double xStep = ((double)imageSpace->width()) / tmpImage.width();
	p.y = ((double)imageSpace->height()) * (tmpImage.height() - 1) / tmpImage.height();

	// Compute bounds on the image that we get back
	double fMax = -MAX_DOUBLE;
	double fMin = MAX_DOUBLE;

	for(int j=0;j<tmpImage.height();j++) {
		p.x = 0.0;
		for(int i=0;i<tmpImage.width();i++) {
			double value = 0;

			switch(inSSDisplay->value()) {
		 case 0 : 
			 // Image itself
			 tmpImage(i,j) = imageSpace->getPixel(p.x,p.y);
			 break;
		 case 1 : 
			 // Blur value
			 tmpImage(i,j) = imageSpace->getValue(p,scale);
			 break;
		 case 2 : 
			 // Get gradient magnitude
			 tmpImage(i,j) = imageSpace->getGradient(p,scale).twoNorm();
			 break;
		 case 3 : 
			 // Laplacean filter
			 tmpImage(i,j) = imageSpace->applyLaplaceanKernel(p,scale);
			 break;
		 case 4 : 
			 // LPP
			 tmpImage(i,j) = imageSpace->applyLPPKernel(p,scale);
			 break;
			}

			fMin = (fMin < tmpImage(i,j)) ? fMin : tmpImage(i,j);
			fMax = (fMax > tmpImage(i,j)) ? fMax : tmpImage(i,j);

			p.x += xStep;
		}

		p.y += yStep;
	}

	// Calculate the image intensity multiplier
	cout << "fMax = " << fMax << endl;
	// cout << "fMin = " << fMin << endl;
	double iMult = 255.0 / (fMax - fMin);
	iMult = 1;

	// Now compute the byte image, stretching f values to intensity range.  
	if(zoom > 1) {
		int jTmp=0;
		for(int j=0;j<ih;j+=zoom) { 
			// Compute the row with image data
			int iTmp = 0;

			for(int i=0;i<iw;i+=zoom) {
				// Compute pixel value
				image(i,j) = (unsigned char) ((tmpImage(iTmp,jTmp)-fMin)*iMult);

				// Copy to the next few pixels
				for(int i1=i+1;(i1<i+zoom) && (i1 < iw);i1++)
					image(i1,j) = image(i,j);

				iTmp++;
			}

			for(int j1=j+1;(j1<j+zoom) && (j1<ih);j1++) {
				for(int i=0;i<iw;i++) 
					image(i,j1) = image(i,j);
			}

			jTmp++;
		}
	}

	else {
		for(int j=0;j<ih;j++) 
			for(int i=0;i<iw;i++)
				image(i,j) = (unsigned char) ((tmpImage(i,j)-fMin)*iMult);
	}
}

void ScaleSpaceWindow::setZoom(int newZoom) {
	zoom = newZoom;
	setScale(scale);
}

/*************************************************************
Profile Window Methods
***********************************************************/
void IProfileWindow::draw() { 
	int i;

	if (!valid()) {
		// Set up the viewport
		glViewport(0,0,w(),h());
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}

	// Fill the window
	glClearColor(0.75,0.75,0.75,0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glFlush();
	glClearColor(0,0,0,0);

	static int dl = -1;

	if(dl < 0) {
		dl = glGenLists(1);
		glNewList(dl,GL_COMPILE);

		// Vertical marks
		for(double t=-0.6666;t<1.00;t+=0.3334) {
			if(fabs(t) < 0.1)
				glColor3d(0.25,0.25,0.25);
			else
				glColor3d(0.65,0.65,0.65);

			glBegin(GL_LINES);
			glVertex2d(t,-1.1);
			glVertex2d(t,1.1);
			glEnd();
		}

		// Horizonal marks      
		glBegin(GL_LINES);
		for(t=1;t<10000;t*=10) {
			if(t==1)
				glColor3d(0.25,0.25,0.25);
			else
				glColor3d(0.65,0.65,0.65);

			double l = (log10(t)/3) * 2.0 - 1.0;
			glVertex2d(-1,l);
			glVertex2d(1,l);
		}
		glEnd();

		glEndList();
	}

	glPushMatrix();
	glScaled(0.95,0.95,1.0);   

	glCallList(dl);


	if(imageSpace==NULL || !display || regOptions.getBooleanValue("display.profile.show",1) == 0) {
		glPopMatrix();
		return;
	}

	samples = regOptions.getIntValue("display.profile.samples",40);

	// Steps to take and from where to take them
	Vector2D xStep = ba.n * (6 * ba.scale / (samples-1));
	BAtom baSample = ba;
	baSample.x -= xStep * (samples/2);

	// Place to store the medialness
	static vector<double> mData(samples);

	// Compute the medialness and intensity
	for(i=0;i<samples;i++) {
		// Now compute medialness at the point.
		mData[i] = -iMatchComputer.compute(baSample,*imageSpace) / 6;
		baSample.x += xStep;
	};

	// Now draw the medialness graph
	double uStep = 2.0 / (samples-1);
	double u=-1;

	glLineWidth(2);
	glBegin(GL_LINE_STRIP);
	for(i=0;i<samples;i++) {
		u+=uStep;

		double mval = mData[i]; //(mData[i]-mMin) / (mMax-mMin);
		glColor3d(1.0,mval,0.0);
		glVertex2d(u,2.0*mval-1);
	}
	glEnd();

	glLineWidth(1.0);
	glPopMatrix();

	glFlush();
}

/*************************************************************
Image Gradient Display methods
***********************************************************/
void IGradientWindow::draw() {
	int w = this->w();
	int h = this->h();

	if (!valid()) {
		// Set up the viewport
		glViewport(0,0,w,h);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}

	// Fill the window
	glClearColor(0.75,0.75,0.75,0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glFlush();
	glClearColor(0,0,0,0);

	static int dl = -1;
	if(dl < 0) {
		dl = glGenLists(1);
		glNewList(dl,GL_COMPILE);

		// outside box
		glColor3d(0,0,0);
		glBegin(GL_LINE_LOOP);
		glVertex2d(-0.99,-0.99);
		glVertex2d(-0.99,0.99);
		glVertex2d(0.99,0.99);
		glVertex2d(0.99,-0.99);
		glEnd();

		glColor3d(0.0,0.65,0.0);

		// Outside circle
		glBegin(GL_LINE_LOOP);
		for(double t=0;t<360;t+=10)
			glVertex2d(0.8*dcos(t),0.8*dsin(t));
		glEnd();

		// Inside circle
		glBegin(GL_LINE_LOOP);
		for(t=0;t<360;t+=10)
			glVertex2d(0.4*dcos(t),0.4*dsin(t));
		glEnd();

		// X,Y axis
		glColor3d(0.0,0.65,0.0);
		glBegin(GL_LINES);
		glVertex2d(-1,0);
		glVertex2d(1,0);
		glVertex2d(0,-1);
		glVertex2d(0,1);
		glEnd();

		glEndList();
	}

	// Scale space is required
	if(imageSpace==NULL || !display) {
		glCallList(dl);
		return;
	}

	// Now, compute the gradient vector
	Vector2D pImage = ba.x.scaledBy(imageSpace->width(),imageSpace->height());
	double sImage = imageSpace->width()*ba.scale;

	Vector2D g = imageSpace->getGradient(pImage,sImage);
	double gMag = g.normalize();

	// Get the type of image that's displayed
	string type = regOptions.getStringValue("display.magnifier.image","pixels");

	if(type != "nothing") {
		// Create an array of pixels to display on the screen
		float *pixels = new float[w*h];
		memset(pixels,0,sizeof(float));

		// Populate the array with blur values from the image
		if(type == "pixels" || type == "blurred" || type == "gradient") {

			Vector2D xStart(pImage.x - sImage*2.5,pImage.y + sImage*2.5);
			Vector2D xStep(sImage*5.0/w,-sImage*5.0/h);
			int offset = 0;
			Vector2D xIndex = xStart;

			for(int j=0;j<h;j++) {
				if(xIndex.y >= 0 && xIndex.y < imageSpace->height()) { 
					for(int i=0;i<w;i++) {
						if(xIndex.x >= 0 && xIndex.x < imageSpace->width()) {
							switch(type[0]) {
							case 'p':
								pixels[offset] = imageSpace->getPixel(xIndex.x,xIndex.y) / 256.0;
								break;
							case 'b':
								pixels[offset] = imageSpace->getValue(xIndex,sImage) / 256.0;
								break;
							case 'g':
								pixels[offset] = imageSpace->getGradient(xIndex,sImage).twoNorm() / 256.0;
								break;
							}
						}
						offset++;
						xIndex.x += xStep.x;
					}
				}
				else {
					offset += w;
				}
				xIndex.y += xStep.y;
				xIndex.x = xStart.x;
			}
		}

		else if (type == "likelihood") {

			BAtom baStart,baSample = ba;
			baStart.x.x -= ba.scale*2.5;
			baStart.x.y += ba.scale*2.5;
			baSample = baStart;

			Vector2D xStep(ba.scale*5.0/w,-ba.scale*5.0/h);
			int offset = 0;

			for(int j=0;j<h;j++) {
				if(baSample.x.y >= 0 && baSample.x.y < 1.0) { 
					for(int i=0;i<w;i++) {
						if(baSample.x.x >= 0 && baSample.x.x < 1.0) {
							pixels[offset] = iMatchComputer.compute(ba,*imageSpace);
						}
						offset++;
						baSample.x.x += xStep.x;
					}
				}
				else {
					offset += w;
				}
				baSample.x.y += xStep.y;
				baSample.x.x = baStart.x.x;
			}
		}

		// Draw pixels on the screen
		glDrawPixels(w,h,GL_LUMINANCE,GL_FLOAT,pixels);
		delete pixels;
	}

	// Draw grid on top of the image
	glCallList(dl);

	// Draw Lines to gradient and back
	//glEnable(GL_LINE_SMOOTH);

	glBegin(GL_LINES);
	glColor3d(1,1,0);
	glVertex2d(g.x,-g.y);
	glVertex2d(-g.x,g.y);
	glVertex2d(0,0);
	glColor3d(1,0,0);
	glVertex2d(ba.n.x,-ba.n.y);
	glVertex2d(0,0);
	glEnd();
	//glDisable(GL_LINE_SMOOTH);

	// Draw a yellow circle for the gradient
	drawOutlineCircle(Vector2D(0.8*g.x,-0.8*g.y),0.1,Gl_Color::rgb(1,1,0));
	drawOutlineCircle(Vector2D(0.8*-g.x,0.8*g.y),0.1,Gl_Color::rgb(1,1,0));

	// Draw a red circle for the primitive
	drawOutlineCircle(Vector2D(0.8*ba.n.x,0.8*-ba.n.y),0.1,Gl_Color::rgb(1,0,0));

	glFlush();
}



void DisplayWindow::drawSpline(RegularSplineSample *sample)
{

	double ct = 0.7;

	glPushAttrib(GL_LINE_BIT);
	glLineWidth(1.0);
	glEnable(GL_LINE_SMOOTH);

	list<SplineSamplePoint>::iterator it,it1,itFirst,itLast,itEnd;
		
	// Frequency of sail vectors
	int sailFreq = regOptions.getIntValue("display.spline.sailVectorFrequency",10);

	for(unsigned int iCurve=0;iCurve<sample->samples.size();iCurve++) {
		for(unsigned int iSeg=0;iSeg<sample->samples[iCurve].size();iSeg++) {

			unsigned int j;
			
			// A spline curve
			list<SplineSamplePoint> &atoms = sample->samples[iCurve][iSeg].points;
			itFirst = atoms.begin();
			itEnd = itLast = atoms.end();
			itLast--;

			// Fat line for the medial axis
			glLineWidth(5.0);		
			glColor3d(ct,1.0,ct);
			
			//glColor3d(0,0.5 + 0.5*(i%2),0);
			//glBegin(GL_POINTS);


			glBegin(GL_LINE_STRIP);			
			for(it = itFirst;it!=itEnd;it++) {
				glVertex2d(it->atom.x());
			}
			glEnd();
			
			// Thin line for the boundary
			glLineWidth(5.0);
			glColor3d(1,0.3,0);

			// Compute the boundary 
			glBegin(GL_LINES);		
			for(it = itFirst;it!=itLast;) {
				if(it->cv[0] < 0)
					glColor3d(ct,ct,1.0);
				else 
					glColor3d(1,ct,ct);
				
				glVertex2d(it->atom.x(MAtom::LEFT));
				glVertex2d((++it)->atom.x(MAtom::LEFT));
			}
			
			// glEnd();
			// glBegin(GL_LINES);		

			for(it = itFirst;it!=itLast;) {
				if(it->cv[1] < 0) 
					glColor3d(ct,ct,1.0);
				else 
					glColor3d(1,ct,ct);				
				glVertex2d(it->atom.x(MAtom::RIGHT));
				glVertex2d((++it)->atom.x(MAtom::RIGHT));
			}
			glEnd();

			// Thin lines
			glLineWidth(1.0);

			// Lines to the boundary
			glBegin(GL_LINES);	
			for(it = itFirst,j=0;it!=itEnd;it++,j++) {
				if(j % sailFreq == 0 || it==itLast) {
					
					if(it->cv[0] < 0 && it != itLast)
						glColor3d(0,ct,1.0);
					else
						glColor3d(1,ct,0);
					
					glVertex2d(it->atom.x(MAtom::LEFT));			
					glVertex2d(it->atom.x());
					/*
					if(it->cv[1] < 0)
						glColor3d(0,ct,1.0);
					else
						glColor3d(1,ct,0);
					*/
					glVertex2d(it->atom.x());
					glVertex2d(it->atom.x(MAtom::RIGHT));			
				}
			}
			glEnd();
		}
		
		// Draw the control points
		for(int i=0;i<sample->spline->curves[iCurve]->size();i++) {
			glLineWidth(4.0);
			glColor3d(ct,1,ct);
			double r = 0.008*winGL->winSize;
			// glColor3d(1,1,0);
			MySMLVec3f C = sample->spline->curves[iCurve]->getControlVec(i);
			//drawOutlineCircle(Vector2D(C.x,C.y),0.008*winGL->winSize,Gl_Color::green,Gl_Color::yellow);
			glBegin(GL_LINES);
			glVertex2d(C.x - r,C.y - r);
			glVertex2d(C.x + r,C.y + r);
			glVertex2d(C.x + r,C.y - r);
			glVertex2d(C.x - r,C.y + r);
			glEnd();

		}
	}

	glPopAttrib();
}
