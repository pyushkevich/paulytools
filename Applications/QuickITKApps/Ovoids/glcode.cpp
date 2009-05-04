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

#include <FL/Fl_Window.H>
#include <FL/glut.H>
#include <GL/glu.h>
#include "fracviewer.h"
#include "blobmodl.h"
#include "implicit.h"
#include "plib.h"

Fl_Window *win3d;
PList *iSurfList = new PList();

double xCursor=1,yCursor=1,zCursor=1;

int dlAxis;

extern Model *model;

// Current point coordinates and data
PointData *currentPoint = NULL;

extern void fillCurrentPointInfo();

bool drawSpheres = false;

/**
 * Draw the current point (point user selected by hitting space bar).
 */
void drawCurrentPoint() {
	if(!currentPoint)
		return;

	glPushMatrix();
	glTranslated(currentPoint->x[0],currentPoint->x[1],currentPoint->x[2]);

	GLUquadricObj *sphere = gluNewQuadric();
	glColor3d(0,0,1);
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	gluSphere(sphere,0.02,10,10);
	gluDeleteQuadric(sphere);

	if(drawSpheres) {
		GLUquadricObj *disk;
		glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

		glPushMatrix();
		disk = gluNewQuadric();
		glTranslated(-currentPoint->f[2][0]/currentPoint->kappa1,
						 -currentPoint->f[2][1]/currentPoint->kappa1,
					    -currentPoint->f[2][2]/currentPoint->kappa1);
		
		gluSphere(disk,fabs(1/currentPoint->kappa1),160,160);
		gluDeleteQuadric(disk);
		glPopMatrix();

		glPushMatrix();
		disk = gluNewQuadric();
		glTranslated(-currentPoint->f[2][0]/currentPoint->kappa2,
			          -currentPoint->f[2][1]/currentPoint->kappa2,
						 -currentPoint->f[2][2]/currentPoint->kappa2);

		gluSphere(disk,fabs(1/currentPoint->kappa2),160,160);
		gluDeleteQuadric(disk);
		glPopMatrix();
	}

	glBegin(GL_LINES);

	if(!currentPoint->umbillic) {
		glVertex3d(0,0,0);
		glVertex3d(currentPoint->f[0][0]/10,currentPoint->f[0][1]/10,currentPoint->f[0][2]/10);

		glVertex3d(0,0,0);
		glVertex3d(currentPoint->f[1][0]/10,currentPoint->f[1][1]/10,currentPoint->f[1][2]/10);
	}

	glVertex3d(0,0,0);
	glVertex3d(currentPoint->f[2][0]/10,currentPoint->f[2][1]/10,currentPoint->f[2][2]/10);

	glEnd();
	glPopMatrix();
}

extern GLfloat EyeDist, EyeEl, EyeAz;

void myDisplay() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glFlush();
	
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glPushMatrix();  
	agvViewTransform();
	
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();  

	
	
	// Place the light source at my location
	glPushMatrix();

	glRotatef(-EyeAz, 0, 1, 0);
	glRotatef(-EyeEl, 1, 0, 0);
	glTranslatef(0,0,EyeDist);

	//GLfloat light_position[] = { 10.0f, 10.0f, 10.0f, 1.0f };
	GLfloat light_position[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	glPopMatrix();

	// Cetnter on current point
	glPushMatrix();

	if(currentPoint) {
		glTranslated(-currentPoint->x[0],-currentPoint->x[1],-currentPoint->x[2]);
	}
	
	glColor3d(1,1,1);
	glCallList(dlAxis);

	for(int i=0;i<	iSurfList->getSize();i++) {
		((ISurface *)iSurfList->get(i))->display();
	}

	// Display the current point
	drawCurrentPoint();

	glPopMatrix();

	//glClearColor(1,1,1,1);
	glPopMatrix();
	
	//glutSwapBuffers();
	glFlush();
}

void myReshape(int w,int h) {
	glViewport(0,0,w,h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, (GLdouble)w/h, 0.01, 100);
	glPushMatrix();
	glMatrixMode(GL_MODELVIEW);
	glFlush();
}

void myGLInit() {
	// Lights, camera, action
	GLfloat light_ambient[] = { 0.4f, 0.4f, 0.4f, 1.0f };
	GLfloat light_diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat light_position[] = { 50.0f, -100.0f, 50.0f, 0.0f };

	GLfloat lmodel_ambient[] = { 0.1f, 0.1f, 0.1f, 1.0f };
	
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);

	glEnable(GL_NORMALIZE);
	glShadeModel(GL_FLAT);

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE);

	// Specular stuff
	//GLfloat specref[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	//glMaterialfv(GL_FRONT,GL_SPECULAR,specref);
	//glMaterialf(GL_FRONT,GL_SHININESS,128);
	

	glFlush();
}

void selectPoint(int x,int y) {
	double x1,y1,z1;
	double x2,y2,z2;

	GLdouble modelmx[16],projmx[16];
	GLint viewmx[4];

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();  
	
	if(currentPoint) {
		glTranslated(-currentPoint->x[0],-currentPoint->x[1],-currentPoint->x[2]);
	}

	glGetDoublev(GL_MODELVIEW_MATRIX,modelmx);
	glGetDoublev(GL_PROJECTION_MATRIX,projmx);
	glGetIntegerv(GL_VIEWPORT,viewmx);
	
	gluUnProject(x,viewmx[3]-y,0,modelmx,projmx,viewmx,&x1,&y1,&z1);
	gluUnProject(x,viewmx[3]-y,1,modelmx,projmx,viewmx,&x2,&y2,&z2);

	glPopMatrix();

	// Now we have a ray from eyepoint to somewhere.  Take small steps along that ray
	// until we hit the object
	if(currentPoint)
		delete currentPoint;
	currentPoint = NULL;
	
	double length = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
	double step = 0.05;

	double xStep = step*(x2-x1)/length;
	double yStep = step*(y2-y1)/length;
	double zStep = step*(z2-z1)/length;

	double lastValue = model->getFunction(x1,y1,z1);
	double stepValue = lastValue;
	
	for(int i=0;i<1000;i++) {
		x1+=xStep;y1+=yStep;z1+=zStep;

		stepValue = model->getFunction(x1,y1,z1);

		if(stepValue <= 0 && lastValue > 0 || stepValue > 0 && lastValue <= 0) {
			currentPoint = new PointData;

			// Root trap here.
			model->rootTrap(x1-xStep,y1-yStep,z1-zStep,
				             x1,y1,z1,
								 currentPoint->x[0],currentPoint->x[1],currentPoint->x[2],0.00000001);
			
			// Compute the point data
			currentPoint->compute(model);

			// Display was updated
			glutPostRedisplay();

			// Call UI method to fill out the current point information
			fillCurrentPointInfo();

			return;
		}
		
		lastValue = stepValue;
	}
}

void myKeyFunc(uchar c,int x,int y) {
	switch(c) {
		// Select the point under the cursor
		case ' ' : selectPoint(x,y);break;
		default:
			agvHandleKeys(c,x,y);
	}
}


/**
 * Start up glut window
 */
void initGLDisplay() {
	win3d = new Fl_Window(0,0,512,512,"3D Display");

	// glut will die unless parent window visible
	// this will cause Glut window to be a child
	win3d->show();				
	win3d->begin();			

	glutInitWindowSize(512, 512);

	// place it inside parent window
	glutInitWindowPosition(0,0); 

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE);

	glutCreateWindow("Title");
	win3d->end();

	win3d->resizable(glut_window);

	// We don't currently have our own idle function
	agvInit(1); 

	glutReshapeFunc(myReshape);
	glutDisplayFunc(myDisplay);

	glutKeyboardFunc(myKeyFunc);

	dlAxis = glGenLists(1);
	agvMakeAxesList(dlAxis);

	myGLInit();
}

