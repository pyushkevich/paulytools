#include "sketch.h"
#include "SketchWindow.h"
#include "sketchui.h"

#ifndef M_PI
#define M_PI 3.141592
#endif

void SketchWindow::draw() {
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

	// Draw the outline
	if(sketch && sketch->outline.spline) {
		// Draw the points and the normal vectors
		glColor3d(1,0.7,0.7);
		glBegin(GL_LINE_LOOP);
		for(double t=0.0;t<1.0;t+=0.01) {
			Vector2D x;
			sketch->outline.getInterpolation(t,x);
			glVertex2d(x.x,x.y);
		}
		glEnd();
/*
		glColor3d(0.6,0.6,0.6);
		glBegin(GL_LINES);
		for(i=0;i<sketch->outline.size();i++) {
			glVertex2d(sketch->outline[i](0),sketch->outline[i](1));
			glVertex2d(sketch->outline[i](0) + 0.03 * sketch->outline[i](2),
				sketch->outline[i](1) + 0.03 * sketch->outline[i](3));
		}
		glEnd(); */
	}

	// Draw the outline
	else if(sketch && sketch->sketch.vertices.size()) {
		// Draw the points and the normal vectors
		glColor3d(1,1,1);
		glBegin(GL_LINE_STRIP);
		for(int i=0;i<sketch->sketch.vertices.size();i++) {
			glVertex2d(sketch->sketch.vertices[i].x.x,sketch->sketch.vertices[i].x.y);
		}
		glEnd();
	}

	// Flush this crap
	glFlush();
}


int SketchWindow::handle(int event) {
	// If not initialized, get out of here
	if(!sketch)
		return 0;

	// Event coordinates in viewport coords
	double ex = 2.0 * Fl::event_x() / w() - 1.0;
	double ey = 1.0 - 2.0 * Fl::event_y() / h();

	switch(event) {    
   case FL_PUSH:
	   sketch->start(ex,ey);
	   return 1;

   case FL_DRAG:      
	   sketch->plot(ex,ey);
	   sketch->ui->winSketch->redraw();
	   return 1;

   case FL_RELEASE:          
	   sketch->end(ex,ey);
	   return 1;

   case FL_KEYBOARD:
	   // ... keypress, key is in Fl::event_key(), ascii in Fl::event_text()
	   return 1;    
   default:
	   // tell fltk that I don't understand other events      
	   return 0;    
	}  

}


void SketchTool::start(double x,double y) {
	sketch.vertices.clear();
	outline.reset();

	sketch.addVertex(Vector2D(x,y),Vector2D(0,0));
}

Vector estimateNormal(Vector x,Vector y) {
	Vector t = y - x;
	return Vector(2,-t(1),t(0));
}

void SketchTool::plot(double x,double y) {
	Vector2D X(x,y);
	Vector2D dx = sketch.vertices.back().x - X;

	if(dx.twoNorm() > 0.002)
		sketch.addVertex(X,Vector2D(0,0));
}

double weight(double z) {
	return exp(-z*z)/sqrt(2*M_PI);
}

void SketchTool::end(double x,double y) {
	plot(x,y);
/*
	// Number of points
	int n = ui->inpSegments->value();

	// Distance standard deviation
	double sdd = 0.1;

	// Compute the outline from the sketch
	for(int i=0;i<n;i++) {
		int k = 1.0 * i * sketch.size() / n;

		double wt = weight(0);
		double wsum = wt;

		Vector normal(2,0.0,0.0);
		Vector point = wt*sketch[k];

		for(int w=1;w<=3;w++) {
			Vector left = sketch[(k-w) % sketch.size()];
			Vector right = sketch[(k+w) % sketch.size()];
			double dl = (left-sketch[k]).twoNorm();
			double dr = (right-sketch[k]).twoNorm();

			wt = weight(dl/sdd);
			point += wt*left;
			normal += wt*estimateNormal(left,sketch[k]);
			wsum += wt;

			wt = weight(dr/sdd);
			point += wt*right;
			normal += wt*estimateNormal(sketch[k],right);
			wsum += wt;
		}

		normal.normalize();
		point /= wsum;

		outline.push_back(Vector(4,point(0),point(1),normal(0),normal(1)));
	}
*/
	// Compute the sketch
	outline.buildFromPB(sketch,(int) ui->inpSegments->value());
	sketch.vertices.clear();

	// Refresh
	ui->winSketch->redraw();
}

void getUserSketch(ContinuousBoundary &cb) {
	SketchTool tool;
	SketchUI ui;

	tool.ui = &ui;

	// Create a user interface
	ui.makeWindow();
	ui.winSketch->sketch = &tool;

	ui.winMain->set_modal();
	ui.winMain->show();
	ui.winSketch->show();

	// Run until closed
	while(ui.winMain->shown())
		Fl::check();

	/*
	PolygonBoundary pb;
	for(int i=0;i<tool.outline.size();i++) {
		pb.addVertex(Vector2D(tool.outline[i](0),tool.outline[i](1)),
			Vector2D(tool.outline[i](2),tool.outline[i](3)));
	}

	cb.buildFromPB(pb,tool.outline.size()); */

	// spline.buildFromPB(pb,tool.outline.size());
	PolygonBoundary pb;
	tool.outline.constructPolygonBoundary(pb,(int) (ui.inpSegments->value()*4));
	cb.buildFromPB(pb,(int)ui.inpSegments->value());
}

