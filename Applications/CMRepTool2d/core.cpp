#include "ui.h"
#include "likehood.h"
#include "reg.h"
#include "undo.h"
#include <SimplexMethod.h>
#include <EvolutionaryStrategy.h>
#include <ConjugateGradientMethod.h>

#ifdef WIN32
#include <limits>
#else 
#include <climits>
#include <limits.h>
#endif

#include <algorithm>

#undef max
#undef min

/****************************************************************
 Core Tracker 
 ****************************************************************/
list <MAtom> core;

// Left and right cores...
// DSLFigure2D coreLeft,coreRight;

// Direction of current tracking effort.  1 = forward,
// -1 = backward, 0 stopped
int coreDirection;

// Flag that is set top true if user interrupts core creation
bool coreInterrupt = false;

// We save a copy of the regular undo buffer
UndoBuffer<MGraph> *savedUndoBuffer = NULL;
MGraph savedModel;

void toggleWidget(Fl_Widget *w,bool toggle) {
   if(toggle)
      w->activate();
   else
      w->deactivate();
}

void toggleNonCoreWidgets(bool toggle) {
   // Disable undo/redo (otherwise user can screw things up)
   toggleWidget(bUndo,toggle);
   toggleWidget(bRedo,toggle);

   // Disable all menus 
   toggleWidget(mbMenuBar,toggle);

   // Disable all tools that are useless for one primitive
   toggleWidget(bSelect,toggle);
   toggleWidget(bBend,toggle);
   toggleWidget(bStretch,toggle);
   // toggleWidget(bTrackCore,toggle);
   bTrackCore->value((toggle) ? 0 : 1);

   // Clear the core
   core.clear();

   // Redraw
   winMain->redraw();
}

void uiTrackCore(Fl_Button *b, void*) {
   // Display core wizard window
   uiScript.uiCommand("show winCoreWizard");
   
   // Disallow double entry in the method
   if(savedUndoBuffer)
      return;

   // Store the undo stack in a backup thingy
   savedUndoBuffer = undoBuffer;
   undoBuffer = new UndoBuffer<MGraph>(50);
   savedModel = model;
   updateUndoMenus();
   
   // Save the current model and replace it with a single figure / single primitive model
   model.beginTopologyEdit();
	model.removeAll();

	// Create a seed atom in the middle
	MAtom seed(winGL->winPosition.x + winGL->winSize/2,
				  winGL->winPosition.y + winGL->winSize/2,
				  winGL->winSize/10,90,90);
   seed.polarity(0,MAtom::HEAD);
	seed.polarity(0,MAtom::TAIL);
	vector<MAtom> avec(1,seed);
	
	// Add the seed atom to the model
	model.addFigure(avec);

	// Rebuild the model
	model.endTopologyEdit();

   // Switch to the PEM tool.
   uiRunTool(TOOL_PEM);   

   // Disable all controls that interfere with core tracking
   toggleNonCoreWidgets(false);
}

class CoreTrackProblem : public Function {
public:
   // Source primitive
   const MAtom &src;
   MAtom pred;
   
   // Distance between the primitives
   double d;

   double evaluate(const Vector &v) {
      // Apply vector to get a primitive
      MAtom guess;

      // Generate the guess
      applyVector(guess,v);

      // Compute a penalty term
      double pp = computePPrior(guess);      

      // Iteration spent
      evaluationCost++;
      
      // If the prior term is too small, go on
      if(pp == 0) 
         return 1e100;
      
      // Compute medialness at the primitive
      double medness = iMatchComputer.compute(guess.bAtom(MAtom::LEFT),*imageSpace) + 
							  iMatchComputer.compute(guess.bAtom(MAtom::RIGHT),*imageSpace);
         
      return medness - log(pp);
   }

   void applyVector(MAtom &trg,const Vector &v) {
      // Direction to next primitive
		Vector2D dx(d*dcos(v(0)),-d*dsin(v(0)));
      trg.x(src.x() + dx);

      // Axial angle
      trg.faDeg(v(1));

      // Radius of next primitive
      trg.r(fabs(v(2)));      

      // Object angle of the next primitive
      trg.oaDeg(v(3));
   }

   void computePrediction() {
      pred = src;
		
		// Position goes along the normal
      pred.translateBy(src.n()*d);

      // Radius changes by 
      double r = src.r() - d*cos(src.oa());
      pred.r(r < 0 ? 0 : r);      
   }

   Vector getStartPosition() {
      Vector v(4);
      v(0) = src.faDeg();
      v(1) = pred.faDeg();
      v(2) = pred.r();
      v(3) = src.oaDeg();
      return v;
   }

   double computePPrior(MAtom &a) {
      // Some vectors for us to work with
      Vector2D d0 = a.x() - src.x();
      Vector2D d1 = a.x(MAtom::LEFT) - src.x(MAtom::LEFT);
      Vector2D d2 = a.x(MAtom::RIGHT) - src.x(MAtom::RIGHT);
      Vector2D d1_ = d1 / d1.twoNorm();
      Vector2D d2_ = d2 / d2.twoNorm();

      Vector2D dp0 = pred.x() - src.x();
      Vector2D dp1 = pred.x(MAtom::LEFT) - src.x(MAtom::LEFT);
      Vector2D dp2 = pred.x(MAtom::RIGHT) - src.x(MAtom::RIGHT);
      Vector2D dp1_ = dp1 / dp1.twoNorm();
      Vector2D dp2_ = dp2 / dp2.twoNorm();

      // Cosine of angle between n0 and d0.  We want it to be one, or very close to one.
      double cosAlpha = d0.dotProduct(dp0) / (d*d);
      double alpha = acos(cosAlpha)*180/M_PI;

      // Dot product between tangent angles at boundaries and vectors to boundary prims
      double bpCosAlpha1 = d1_.dotProduct(dp1_);
      double bpCosAlpha2 = d2_.dotProduct(dp2_);

      double bpAlpha1 = acos(bpCosAlpha1)*180/M_PI;
      double bpAlpha2 = acos(bpCosAlpha2)*180/M_PI;

      double penalty = 1.0;

      // Penalties for angles            
      /*
      penalty *= penaltyFunction(alpha,0,5,10);
      penalty *= penaltyFunction(p.getTO(),pred.getTO(),5,15);
      penalty *= penaltyFunction(p.getTA(),pred.getTA(),5,10);
      penalty *= penaltyFunction(log(p.getR() / pred.getR()),0,log(1.2));
      penalty *= penaltyFunction(bpAlpha1,0,5,10);
      penalty *= penaltyFunction(bpAlpha2,0,5,10);
      */
      
      penalty *= penaltyFunction(alpha,0,10,20);
      penalty *= penaltyFunction(a.oaDeg(),pred.oaDeg(),0,5);
      penalty *= penaltyFunction(a.faDeg(),pred.faDeg(),10,20);
      penalty *= penaltyFunction(log(a.r() / pred.r()),0,log(1.1),log(1.2));
      penalty *= penaltyFunction(bpAlpha1,0,15,30);
      penalty *= penaltyFunction(bpAlpha2,0,15,30);
      penalty *= penaltyFunction(a.oaDeg(),90,10,30);

/*      
      penalty *= penaltyFunction(alpha,0,40,45);
      penalty *= penaltyFunction(p.getTO(),pred.getTO(),40,45);
      penalty *= penaltyFunction(p.getTA(),pred.getTA(),40,45);
      penalty *= penaltyFunction(log(p.getR() / pred.getR()),log(1.2),log(1.5));
      penalty *= penaltyFunction(bpAlpha1,0,40,45);
      penalty *= penaltyFunction(bpAlpha2,0,40,45);
*/      

      return penalty;
   }

   // dToNext - distance to next primitive, may be negative depending on tracking direction
   CoreTrackProblem(const MAtom &source,double dToNext) : src(source) {
      d = dToNext;
      computePrediction();
   }
};

/*
class CoreStartProblem : public Function {
public:
   // Constructor
   CoreStartProblem(DSLPrimitive2D &prim,ImageSpace &inSpace);
   
   // Evaluate function
   double evaluate(const Vector &x);
   
   // Apply/unapply vector to a primitive
   DSLPrimitive2D &applyVector(const Vector &x);
   
private:
   // Primitive that we will be working with
   DSLPrimitive2D prim0,primT;

   // Radius - unit for all computations
   double r;

   // Tangent vectors at boundary primitives
   Vector2D t0[2];

   // Image space that we will use
   ImageSpace &space;

   // Compute the penalty term
   double computePenalty();
};

CoreStartProblem::CoreStartProblem(DSLPrimitive2D &prim,ImageSpace &inSpace) :
prim0(prim),space(inSpace)
{
   r = prim0.getR();
   t0[0] = prim0.getNormal(1).getNormal();
   t0[1] = prim0.getNormal(2).getNormal();
}

DSLPrimitive2D &CoreStartProblem::applyVector(const Vector &x) {
   primT = prim0;
   
   primT.translateBy(x(0)*r,x(1)*r);
   primT.scaleBy(pow(2,x(2)));
   primT.rotateBy(x(3));
   primT.setTO(primT.getTO()+x(4));
   return primT;
}

double CoreStartProblem::evaluate(const Vector &x) {
   applyVector(x);
   evaluationCost++;
   
   double penalty = computePenalty();

   if(penalty < 0.001) {
      return penalty;
   }
   else {
      double likelihood = iMatchComputer.cmpDiscMedAt(primT,space,false);
      return likelihood * penalty;
   }
}

double CoreStartProblem::computePenalty() {
   double penalty = 1.0;

   // The penalty is based on displacements of boundary primitives.
   // Compute the normals and stuff
   
   for(int i=1;i<=2;i++) {
      // Displacement
      Vector2D dx = primT.getPosition(i) - prim0.getPosition(i);

      // Movement along normal
      double mn = fabs(dx.dotProduct(prim0.getNormal(i)));

      // Movement along tangent
      double mt = fabs(dx.dotProduct(t0[i-1]));
      
      // Cosine of angle between normals
      double dn = 1.0 - prim0.getNormal(i).dotProduct(primT.getNormal(i));

      // Compute penalty
      penalty *= penaltyFunction(mn,0,0.7*r);
      penalty *= penaltyFunction(mt,0,0.2*r);
      penalty *= penaltyFunction(dn,0,1.0 - dcos(20.0));
   }

   return penalty;
}

void corePlaceSeed() {
   // Core seed
   DSLPrimitive2D seed = model.getFigure(0).getPrimitive(0);

   // A common solution space
   CoreStartProblem problem(seed,*imageSpace);
   GaussianSS ss(Vector(5),Vector(5,0.25,0.25,0.5,30.0,30.0));
   EvolutionaryStrategy cgm(problem,ss,2,4,SELECTION_MuPlusLambda);
   
   // NumericalFunction nProblem(problem,0.0001);
   // ConjugateGradientMethod cgm(nProblem,Vector(5));

   while(!cgm.isFinished() && problem.getEvaluationCost() < 1000)  
      cgm.performIteration();

   model.getFigure(0).update(0,problem.applyVector(cgm.getBestEverX()));
}
*/

void toggleExtractionControls(bool toggle) {
   toggleWidget(bCoreStart,toggle);
   toggleWidget(bCoreStop,!toggle);
   toggleWidget(bCoreClear,toggle);
   toggleWidget(bCoreFinish,toggle);
   toggleWidget(bCoreCancel,toggle);
   toggleWidget(inCoreStart,toggle);
   toggleWidget(inCoreEnd,toggle);
   toggleWidget(inCoreSample,toggle);
   toggleWidget(inCoreInterval,toggle);
}

void uiCoreStart(Fl_Button*, void*) {
   // Enable stop button
   toggleExtractionControls(false);

   // Optimize the primitive first.
	runNodeOptimization();

   // Add primitive to the core.  The model better have an atom...   
   if(model.nodes().size() != 1)
		throw "The model must have one atom for core tracking!";
	core.clear();
	core.push_back(*model.nodes().front());

   // Start in forward direction
	coreDirection = 1;

	// Get some registry values first
	double coreInterval = regOptions.getDoubleValue("coreTracker.interval",1.0) / imageSpace->width();

   while(coreDirection) {
      // Distance to next value (sign given by direction, magnitude user spcified)
      double dNext = coreDirection * coreInterval;

      // Create a core tracking problem
      CoreTrackProblem ctp(core.back(),dNext);

      // Setup the solution space
      GaussianSS ss(ctp.getStartPosition(),Vector(4,2.0,2.0,fabs(dNext/10),2.0));
      NumericalFunction nProblem(ctp,0.0001);

		// Create an optimizer
		NumericalMethod *method = createOptAlg(nProblem,ss);

      // Run for a while
		while(ctp.getEvaluationCost() < 5000 && !method->isFinished()) {
         method->performIteration();
      }

      // Run the event loop
      if (Fl::ready())
         Fl::check();               

      // If an interruption happened, reverse directions
      if(coreInterrupt) {
         coreInterrupt = false;
         if(coreDirection == 1) {
            // Switch direction
            coreDirection = -1;

            // Swap the order of primitives in the core
				reverse(core.begin(),core.end());
         }
         else {
            coreDirection = 0;
         }
      }
      else {
         // Add solution to the end of core
         MAtom soln;
         ctp.applyVector(soln,method->getBestEverX());
         core.push_back(soln);

         // Cores have changed, so we need to redisplay
         winGL->redraw();
      }

		delete method;
   }

   // OK, so we are done extracting now.  
   // Enable all buttons
   toggleExtractionControls(true);

   // Update the core start and core end and sampling control
   inCoreStart->maximum(core.size()-1);
   inCoreEnd->maximum(core.size()-1);
   inCoreSample->maximum(core.size());

   inCoreStart->value(0);
   inCoreEnd->value(core.size()-1);
   if(inCoreSample->value() > core.size())
      inCoreSample->value(core.size());
}

void uiCoreStop(Fl_Button*, void*) {
   // Stop tracking in given direction
   coreInterrupt = true;
}

void uiCoreClear(Fl_Button*, void*) {
   // Clear the currently extracted core
   core.clear();

   // Update the coreStart,coreEnd controls
   toggleWidget(inCoreStart,false);
   toggleWidget(inCoreEnd,false);
   toggleWidget(inCoreSample,false);
   toggleWidget(bCoreClear,false);

   // Redisplay
   winGL->redraw();
}

void uiCoreEndChange(Fl_Value_Slider*, void*) {
   // Set the max for start to be the end
   inCoreStart->maximum(inCoreEnd->value());
   inCoreEnd->minimum(inCoreStart->value());
   if(inCoreStart->value() > inCoreStart->maximum())
      inCoreStart->value(inCoreStart->maximum());
   if(inCoreEnd->value() < inCoreEnd->minimum())
      inCoreEnd->value(inCoreEnd->minimum());

   // Compute the number of primitives in core
   int nCorePrims = inCoreEnd->value()-inCoreStart->value()+1;

   // Set up the samples slider
   if(inCoreSample->value() == inCoreSample->maximum() || inCoreSample->value() > nCorePrims) {
      inCoreSample->value(nCorePrims);
   }
   inCoreSample->maximum(nCorePrims);

   // Redisplay
   winGL->redraw();
   winCoreWizard->redraw();
}

// This method draws the core using open GL
void drawCore() {   
   glPointSize(3.0);
   
   glBegin(GL_POINTS);
	
	// This guy iterates through the atoms
	list<MAtom>::const_iterator coreItr = core.begin(),start,stop;

   for(unsigned int i=0;i<core.size();i++) {
      if(inCoreStart->active() && inCoreStart->value() <= i &&  inCoreEnd->value() >= i)
         glColor3d(1,1,0);
      else
         glColor3d(1,0,0);
      
      glVertex2d(coreItr->x());
      
      glColor3d(0.6,0.4,0.8);
      glVertex2d(coreItr->x(MAtom::LEFT));

      glColor3d(0.4,0.6,0.8);
      glVertex2d(coreItr->x(MAtom::RIGHT));

		// Remember the start and stop atoms
		if(i == inCoreStart->value())
			start = coreItr;
		if(i == inCoreEnd->value())
			stop = coreItr;

		coreItr++;
   }

   glEnd();

   // Draw primitives at end points
   if(inCoreStart->active()) {
      drawSingleAtom(*start);
      drawSingleAtom(*stop);
   }

   glPointSize(1.0);   
}

void uiCoreCancel(Fl_Button*, void*) {
   // Clear the core
   uiCoreClear(NULL,NULL);
   
   // Restore the old undo buffer
   delete undoBuffer;
   undoBuffer = savedUndoBuffer;
   savedUndoBuffer = NULL;
   updateUndoMenus();

   // Get the model back too
   model = savedModel;

   // Hide the wizard
   winCoreWizard->hide();

   // Redraw
   onModelOrSelectionChange();
	winGL->redraw();
   
   // Enable all features disabled by the core tracker
   toggleNonCoreWidgets(true);
}

void uiCoreFinish(Fl_Button*, void*) {
   // Sample the core to form a figure.  I use nearest neighbour sampling, but I could be using
   // interpolation.  The latter would be too slow though.   
   int pStart = inCoreStart->value();
   int pEnd = inCoreEnd->value();
   int k = inCoreSample->value();
	int n = pEnd - pStart + 1;

	// Going from n atoms to k atoms
	vector <MAtom> figure;
	figure.reserve(k);
 
   if((int) core.size() < pEnd) {
      uiCoreCancel(NULL,NULL);
      return;
   }

   double sStep = 1.0 * n / k;
	double nextAdd = pStart;
	double p = 0;
	for(list<MAtom>::const_iterator i = core.begin();i!=core.end();i++) {
		if(p >= nextAdd) {
			figure.push_back(*i);
			nextAdd += sStep;
			if(figure.size() == k)
				break;
		}
		++p;
	}

   // Clear the core
   uiCoreClear(NULL,NULL);
   
   // Restore the old undo buffer
   delete undoBuffer;
   undoBuffer = savedUndoBuffer;
   savedUndoBuffer = NULL;

   // Get the model back too
   model = savedModel;
   pushModel();

	// Add the figure to the model
	model.beginTopologyEdit();
   model.addFigure(figure);
	model.endTopologyEdit();

   updateUndoMenus();

   // Enable all features disabled by the core tracker
   toggleNonCoreWidgets(true);

   // Hide the wizard
   winCoreWizard->hide();

   // Redraw
   onModelOrSelectionChange();
   winGL->redraw();
}
