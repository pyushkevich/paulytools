#include "misc.h"
#include "tools.h"
#include "dslwin.h"
#include "dsltool.h"
#include "ui.h"
#include "ispace.h"

#include "Fl/gl.h"
#include "GL/glu.h"

/*
class MRepColorScheme {
public:
	// Colors
	Gl_Color mAtom,mLink,bLink,bCurve;
	Gl_Color bAtomLeft,bAtomRight,bAtomEnd;

	MRepColorScheme() {
		mAtom = Gl_Color::hsv(hue,1,0,0);
		bAtom = Gl_Color::hsv(hue,0.6,1,1);
		mLink = Gl_Color::hsv(hue,1,1,0.6);
		bLink = Gl_Color::hsv(hue,1,1,0.6);
		bCurve = Gl_Color::hsv(hue,0.6,1,1);
	}

	void batom(Gl_Color c) {
		bAtomLeft = bAtomRight = bAtomEnd = c;
	}

	static MRepColorScheme monochome(int hue) {
		MRepColorScheme sch;
		sch.batom(Gl_Color::hsv(hue,0.6,1,1));

		sch.mAtom = Gl_Color::hsv(hue,1,1,1);
		sch.mLink = Gl_Color::hsv(hue,1,1,0.6);
		sch.bLink = Gl_Color::hsv(hue,1,1,0.6);
		sch.bCurve = Gl_Color::hsv(hue,0.6,1,1);
		return sch;
	}

	static MRepColorScheme standard() {
		return sch;
	}
};
*/

// Computes the radius used to draw the primitive
double computeRadius(const MAtom &atom) {
   return min(0.012*winGL->winSize,atom.r()/3.0);
}


void computeSelectionBox() {
	MNodeGroup *ng = model.selection();
	
	if(ng->size() > 0) {
      ng->getExtents(winGL->selTopLeft,winGL->selBotRight);
      winGL->haveSelection = true;
   }
   else {
      winGL->haveSelection = false;
	}
}

bool locateBAtom(MouseEvent &e,BRef &outBRef) {
   for(MNode *node=model.first();node!=NULL;node=model.next()) {
		// The radius of the center circle
      double rad = computeRadius(*node);
         
      Vector2D pos[4];
      pos[MAtom::MIDDLE] = node->x();      
      pos[MAtom::LEFT] = node->x(MAtom::LEFT);
      pos[MAtom::RIGHT] = node->x(MAtom::RIGHT);
		pos[MAtom::HEAD] = node->x() + node->n() * 2 * rad;
         
		// See if we are inside the circle
		for(int target=MAtom::MIDDLE;target<=MAtom::HEAD;target++) {
			double r = (target==MAtom::MIDDLE) ? rad : rad*0.6;
			if((e.xSpace-pos[target]).twoNorm() < r) {
				outBRef.node = node;
				outBRef.bIndex = target;
				return true;
			}
		}         
   }
   
   outBRef.node = NULL;
   return false;
}



void drawSingleAtom(const MAtom &atom) {
	static Gl_Color clrLeft  = Gl_Color::rgb(0.0,0.4,0.8);
	static Gl_Color clrRight = Gl_Color::rgb(0.0,0.8,0.4);
	static Gl_Color clrMiddle = Gl_Color::rgb(1.0,0.0,0.0);

	double r  = computeRadius(atom);

	glColor3d(.2,.2,1);
	glBegin(GL_LINE_STRIP);
	glVertex2d(atom.x(MAtom::LEFT));
	glVertex2d(atom.x(MAtom::MIDDLE));
	glVertex2d(atom.x(MAtom::RIGHT));
	glEnd();

	glColor3d(1,0,0);
	drawOutlineCircle(atom.x(MAtom::LEFT),r*0.6,clrLeft);
	drawOutlineCircle(atom.x(MAtom::RIGHT),r*0.6,clrRight);
	drawOutlineCircle(atom.x(MAtom::MIDDLE),r*0.8,clrMiddle);
}

void drawPrimitiveEditMode(MNode *node) {
   double r  = computeRadius(*node);
	
	static Gl_Color cLeft = Gl_Color::rgb(0.6,0.4,0.8,0.4);
	static Gl_Color cRight = Gl_Color::rgb(0.4,0.6,0.8,0.4);

	double maxAlpha = 0;
	for(int i=MAtom::LEFT;i<=MAtom::RIGHT;i++) {
		Gl_Color clr = (i==MAtom::LEFT) ? cLeft:cRight;

		// What's the alpha
		double alpha = (node->bAtom(i).bnd) ? 1 : 0.3;

		// Compute the max alpha for displaying the middle
		maxAlpha = alpha > maxAlpha ? alpha : maxAlpha;

		// Draw a connection to the center
		glBegin(GL_LINES);
		
		clr.apply();
		glVertex2d(node->x(MAtom::MIDDLE));
		glVertex2d(node->x(i));
		
		glEnd();

		// Draw the cirle
		drawOutlineCircle(node->x(i),r*0.6,clr.alpha(alpha));
	}

	glColor4d(0.2,0.8,0.2,0.4);
	glBegin(GL_LINES);
	glVertex2d(node->x() + node->n()*2*r);
	glVertex2d(node->x());
	glEnd();

	drawOutlineCircle(node->x() + node->n()*2*r,r*0.6,Gl_Color::rgb(0.0,1.0,0.0,maxAlpha));
   
	if(node->selected()) 
      drawOutlineCircle(node->x(),r*0.8,Gl_Color::rgb(1.0,1.0,0.0,maxAlpha));
   else
      drawOutlineCircle(node->x(),r*0.8,Gl_Color::rgb(1.0,0.0,0.0,maxAlpha));   
}

void draw(MNode *node) {
   double r  = computeRadius(*node);
   
   // Draw the lines 
	double maxAlpha = 0;
	for(int i=MAtom::LEFT;i<=MAtom::TAIL;i++) {
		// Heads are drawn only if we are at the end of the thingy
		if(i == MAtom::HEAD && node->hasSuccessor())
			continue;

		// Tails are drawn only if we are at the end of the thingy
		if(i == MAtom::TAIL && node->hasPredecessor())
			continue;

		// What's the alpha
		double alpha = (node->bAtom(i).bnd) ? 1 : 0.4;

		// Compute the max alpha for displaying the middle
		maxAlpha = alpha > maxAlpha ? alpha : maxAlpha;

		// Draw a connection to the center
		glBegin(GL_LINES);
	   glColor4d(0.2,0.2,1.0,0.4);
		glVertex2d(node->x(MAtom::MIDDLE));
		glVertex2d(node->x(i));
		glEnd();

		// Draw the cirle
		drawOutlineCircle(node->x(i),r*0.6,Gl_Color::rgb(0.4,0.8,0.0,alpha));
	}

   // Draw middle
	if(node->selected()) 
      drawOutlineCircle(node->x(),r*0.8,Gl_Color::rgb(1.0,1.0,0.0,maxAlpha));
   else
      drawOutlineCircle(node->x(),r*0.8,Gl_Color::rgb(1.0,0.0,0.0,maxAlpha));
}
void drawBranch(MNode *node,int bIndex) {
   double r = computeRadius(*node);

   double alpha = (node->bAtom(bIndex).tShape == BAtom::NOT_IN_SHAPE) ? 0.4 : 1.0;
   drawOutlineCircle(node->x(bIndex),r*0.8,Gl_Color::rgb(0.8,0.8,0.8,alpha));
}

void drawBranch(MNode *node) {
   double r  = computeRadius(*node);
   
   // Draw the lines 
   glBegin(GL_LINES);
   
	glColor4d(0.5,0.5,0.5,0.3);

   glVertex2d(node->x(MAtom::LEFT));
   glVertex2d(node->x());
   glVertex2d(node->x(MAtom::RIGHT));
   glVertex2d(node->x());
   
   if(!node->hasPredecessor()) {
      glVertex2d(node->x(MAtom::TAIL));
      glVertex2d(node->x());
   }
   if(!node->hasSuccessor()) {
      glVertex2d(node->x(MAtom::HEAD));
      glVertex2d(node->x());
   }

   glEnd();
   
   // Draw circles
   drawBranch(node,MAtom::MIDDLE);
   drawBranch(node,MAtom::LEFT);
   drawBranch(node,MAtom::RIGHT);

   if(!node->hasPredecessor()) {
      drawBranch(node,MAtom::TAIL);
   }
   if(!node->hasSuccessor()) {
      drawBranch(node,MAtom::HEAD);
   }
}



inline void setBndColor(bool fast,double polarity) {
	if(fast) {
		glColor4d(1.0,1.0,0.0,fabs(polarity));
	}
	else if(polarity < 0) {
		glColor4d(0.0,1.0,1.0,-polarity);
	}
	else {
		glColor4d(0.5,0.0,1.0,polarity);
	}
}


void drawFigureOval(MFigure *figure) {
	Vector2D c,mj,mn;
	figure->fitEllipse(c,mj,mn);
	
	glBegin(GL_LINES);
	glColor3d(0.6,1,0.2);
	glVertex2d(c);
	glVertex2d(c+mj);
	glColor3d(1,0.6,0.2);
	glVertex2d(c);
	glVertex2d(c+mn);
	glEnd();
}

// Draw a model
void draw(MGraph *graph,bool fast,bool pem) {
	MFigureItr f;
	MShapeItr s;
	BoundaryItr b;
	
	// Draw the boundaries
	for(s=graph->shapes().begin();s!=graph->shapes().end();s++) {
		MShape *shape = *s;
		for(b=shape->bounds().begin();b!=shape->bounds().end();b++) {
			Boundary *bound = *b;

			glBegin(GL_LINE_STRIP);
			for(int i=0;i<bound->size();i++) {
				int j = (i+1) % bound->size();

				// Compute the length of the boundary piece
				double len = bound->bRef(i).bAtom()->x.distanceTo(bound->bRef(j).bAtom()->x);

				// Compute the step size
				double nSteps = 200 * len / winGL->winSize;
				nSteps = nSteps > 256 ? 256 : (nSteps < 2 ? 2 : nSteps);
				double step = (fast) ? 4.0 / nSteps : 1.0 / nSteps;

				for(double t=i;t < i+1;t+=step) {
					BAtom ba = bound->estimate(t);
					setBndColor(fast,ba.polarity);
					glVertex2d(ba.x);
					// glVertex2d(ba.x+ba.n * ba.scale);
					//glVertex2d(ba.x);
				}
			}

			// Draw last segment - because my laptop doesn't do line-loop
			BAtom ba = bound->estimate(0);
			setBndColor(fast,ba.polarity);
			glVertex2d(ba.x);

			glEnd();
		}
	}

	// Draw the figures
	for(f=graph->figures().begin();f!=graph->figures().end();f++) {
		MFigure *figure = *f;

		// Compute the interpolation of the figure		
		for(int i=0;i<figure->size()-1;i++) {
			MAxelet &axis = figure->getMedialAxis(i);
	
			// Check if traceable
			//if(axis.isTraceable())
				glColor4d(0.8,0.1,0.1,0.8);
			//else
			//	glColor4d(1.0,1.0,1.0,1.0);


			double step = fast ? 0.1 : 0.01;
			glBegin(GL_LINE_STRIP);
			glVertex2d(figure->node(i+1)->x());
			for(double t=step;t < 1.0;t+=step) {
				// Get an estimation of the atom
				const MAtom &a = axis.getEstimation(t);
				if(a.r() > 0)
					glVertex2d(a.x());
			}
			glVertex2d(figure->node(i)->x());
			glEnd();
		}
		
/*
		// Connect the atoms
		glColor4d(0.8,0.1,0.1,0.4);
		glBegin(GL_LINE_STRIP);
		for(int i=0;i<figure->size();i++) {
			glVertex2d(figure->node(i)->x());
		}
		glEnd();
*/
		// Draw each atom
		for(i=0;i<figure->size();i++) {
			if(pem)
				drawPrimitiveEditMode(figure->node(i));
			else 
				draw(figure->node(i));
		}

		// drawFigureOval(figure);
	}
}

void drawBranchMode(MGraph *graph) {
	MFigureItr f;
	MShapeItr s;
	BoundaryItr b;
	
	// Draw the boundaries
   glColor3d(0.8,0.8,0.8);
	for(s=graph->shapes().begin();s!=graph->shapes().end();s++) {
		MShape *shape = *s;
		for(b=shape->bounds().begin();b!=shape->bounds().end();b++) {
			Boundary *bound = *b;			

			glBegin(GL_LINE_STRIP);
			for(int i=0;i<bound->size();i++) {
				int j = (i+1) % bound->size();

				// Compute the length of the boundary piece
				double len = bound->bRef(i).bAtom()->x.distanceTo(bound->bRef(j).bAtom()->x);

				// Compute the step size
				double nSteps = 200 * len / winGL->winSize;
				nSteps = nSteps > 256 ? 256 : (nSteps < 2 ? 2 : nSteps);
				double step = 1.0 / nSteps;

				for(double t=i;t < i+1;t+=step) {
					BAtom ba = bound->estimate(t);
					glVertex2d(ba.x);
					//glVertex2d(ba.x+ba.n * ba.scale);
					//glVertex2d(ba.x);
				}
			}

			// Draw last segment - because my laptop doesn't do line-loop
			BAtom ba = bound->estimate(0);
			glVertex2d(ba.x);

			glEnd();
		}
	}

	// Draw the figures
	for(f=graph->figures().begin();f!=graph->figures().end();f++) {
		MFigure *figure = *f;

		// Connect the atoms
		glColor4d(0.8,0.8,0.8,0.3);
		glBegin(GL_LINE_STRIP);
		for(int i=0;i<figure->size();i++) {
			glVertex2d(figure->node(i)->x());
		}
		glEnd();

		// Draw each atom
		for(i=0;i<figure->size();i++) {
			drawBranch(figure->node(i));
		}
	}

}

// Returns true if just one figure selected
MFigure *getSingleSelectedFigure() {
   if(model.figures().size() == 1) {
		model.select();
		return model.figures().front();
	}

	MNodeGroup *ng = model.selection();
	if(ng->size() == 0)
		return NULL;

	MFigure *f = ng->node(0)->figure();
	model.deselect();
	f->select();
   return f;
}


