#include <SimplexMethod.h>
#include <EvolutionaryStrategy.h>
#include <ConjugateGradientMethod.h>
#include <mreps2D.h>
#include <interp2D.h>
#include <limits>
#include "likehood.h"
#include "dsltool.h"
#include "undo.h"
#include "ui.h"
#include "reg.h"
#include "time.h"
#include "SplineMrep.h"

#define MAX_DOUBLE 1e100

#include <minmax.h>
#include <algorithm>

inline double log2(double x) {
   static double m = 1.0 / log(2.0);
   return log(x) * m;
}

Timer tPenalty,tEval,tTotal;
extern Timer tBoundlet;

class RegProblem : public Function {
public:
	RegProblem(MNodeGroupIF &group,ImageSpace &inSpace,Vector penalty);
	
	// Evaluate function
	double evaluate(const Vector &v);
	
	// Apply/unapply vector to a model
	void applyVector(const Vector &x);
	
private:
   MNodeGroupIF &ng;
   ImageSpace &space;
   
   Vector oldX;
   Vector penaltySigma;
};


RegProblem::RegProblem(MNodeGroupIF &group,ImageSpace &inSpace,Vector penalty) : 
ng(group),space(inSpace) 
{
	penaltySigma = penalty;
   oldX = Vector(4);
}

void RegProblem::applyVector(const Vector &x) {
   // Remember the current transform
   Vector dx = x - oldX;
	
   // Scale
   ng.scaleBy(pow(2.0,dx(2)));
	
   // Rotate 
   ng.rotateBy(dx(3));
	
   // Translate
   ng.translateBy(Vector2D(dx(0),dx(1)));
	
   // Remember the current displacement
   oldX = x;
}

double RegProblem::evaluate(const Vector &v) {
	// Apply the vector to the model
	applyVector(v);
	
   // Update the eval cost
   evaluationCost++;
	
   // Compute penalty
   double penalty = 1.0;
	
   for(int i=0;i<v.size();i++) {
      penalty *= penaltyFunction(v(i),0,penaltySigma(i));
   }
	
   if(penalty < 0.001)
      return penalty;
	
	// Compute the boundariness along the model
	return penalty * iMatchComputer.compute(ng,space);
   // return -fabs(computeGradientMatch(model,space));
}

/************************************************************************
Node-by-Node optimization
************************************************************************/
class NodePrior {
public:
	virtual ~NodePrior() {
	};
	  
	// Computes prior probability of two nodes
	virtual double priorProbability(MNode *ref,MNode *trg) = 0;
	  
	// Returns the smallest value of prior probability for which the optimizer should
	// still bother to compute the likelihood term
	virtual double minPP() {
		return 0;
	}
	  
	// Returns the solution space in which to run this optimization. 
	virtual SolutionSpace *createSS(MNode *ref);
};

// Utility structure for neighbour priors 
struct NeighborPriorTerms {
   // Some vectors
   Vector2D v12l,v12r;	

   // Some distances
	double d12l,d12r,d12;

   // Some constraints
	double c[7];

   // Compute the prior terms
   void compute(MNode *n1,MNode *n2);
};


// This prior computes penalties based on the relationships between neighboring atoms
// in reference and target models
class NeighborPrior : public NodePrior {
	NeighborPriorTerms npt[2];

   double computePair(NeighborPriorTerms &ref,NeighborPriorTerms &trg);
public:
	bool penalizeSliding;

   NeighborPrior(MNode *reference);

	virtual double priorProbability(MNode *ref,MNode *trg);
	virtual SolutionSpace *createSS(MNode *ref);
};


// This prior computes penalties based on the relationships between neighboring atoms
// in reference and target models
class BoundaryPrior : public NodePrior {
	double compute(const BAtom &p,const BAtom &x,const BAtom &s);
public:
	double priorProbability(MNode *ref,MNode *trg);
	SolutionSpace *createSS(MNode *ref);
};


// This prior computes penalties based on how an atom moves relative to it's previous
// position
class AutoPrior : public NodePrior {
	// These values control the prior.  They are user-definable
	double nrmShift,tngShift,nrmRotn;
	  
public:
	virtual double priorProbability(MNode *ref,MNode *trg);
	  
	// Set the user-defined constraints
	void setMovementLimits(double normalShift,double tangentShift,double normalRotation);
};

class NodeProblem : public Function {
private:
	ImageSpace &m_image;
	NodePrior &m_prior;
	MNode *m_r,*m_t;
	  
	bool moveAlong;

public:
	// Constructor
	NodeProblem(MNode *ref,MNode *trg,ImageSpace &image,NodePrior &prior);		
	double evaluate(const Vector &x);

	void applyVector(const Vector &v) {
		if(moveAlong || m_r->isEnd()) {
			m_t->setData(v,m_t->isEnd());
		}
		else {
			Vector w = v;
			Vector2D xt(w(0),w(1));
			xt = m_r->x() - m_r->n() * ((xt - m_r->x()).dotProduct(m_r->n()));
			w(0) = xt.x;
			w(1) = xt.y;
			m_t->setData(w,m_t->isEnd());
		}
	}

	void setMoveAlong(bool b) {
		moveAlong = b;
	}
};

class MultiNodeProblem : public Function {
private:
	ImageSpace &image;
	vector<NodePrior *> prior;

	// Reference and target node groups
	MNodeGroupIF *ref,*tar;

public:
	// Constructor
	MultiNodeProblem(MNodeGroupIF *ref,MNodeGroupIF *tar,ImageSpace &img,vector<NodePrior *> pri) : image(img),prior(pri)
	{
		this->ref = ref;
		this->tar = tar;
	}

	// Apply a Euclidean vector to a model
	void applyVector(const Vector &v) {
		int k=0;
		for(MNode *n=tar->first();n != NULL;n = tar->next()) {
			k = n->setData(v,n->isEnd(),k);
		}
	}

	// Evaluation fuction 
	double evaluate(const Vector &v) {
		// Apply the vector to the group
		applyVector(v);
		
		// Keeps track of prior probability
		double logPP = 0;

		// Apply prior to each pair of nodes
		MNode *nt = tar->first(), *nr = ref->first();
		int i = 0;
		while(nt && nr) {
			double pp = prior[i++]->priorProbability(nr,nt);
			if(pp <= 0.0) 
				return MAX_DOUBLE;
			logPP += log(pp);

			nt = tar->next();
			nr = ref->next();
		}
				
		// Update counter
		evaluationCost++;
		
		// Compute the medialness   
		double lHood = iMatchComputer.compute(*tar,image);
		
		// Return bayesian match term
		return lHood - logPP;
	}

};

SolutionSpace *NodePrior::createSS(MNode *ref) {
	Vector mean,sd;
	
	// Mean comes from reference node
	mean = ref->getData(ref->isEnd());
	
	// Do we care about ends?
	if(ref->isEnd())
		sd = Vector(6,0.01,0.01,0.005,M_PI/6,M_PI/6,0.1);      
	else
		sd = Vector(5,0.01,0.01,0.005,M_PI/6,M_PI/6);
	
	return new GaussianSS(mean,sd);
}

double AutoPrior::priorProbability(MNode *ref,MNode *trg) {
	double pp = 1.0;
	
	// The penalty is based on displacements of boundary primitives.
	// Compute the normals and stuff
	for(int i=MAtom::MIDDLE;i<=MAtom::RIGHT;i++) {
		// Displacement
		Vector2D dx = trg->x(i) - ref->x(i);
		pp *= penaltyFunction(dx.twoNorm(),0,trg->r()/5,trg->r()/2);
/*
		// Movement along normal
		double mn = fabs(dx.dotProduct(ref->n(i)));
		
		// Movement along tangent
		double mt = fabs(dx.dotProduct(ref->n(i).getNormal()));
		
		// Cosine of angle between normals
		double dn = 1.0 - ref->n(i).dotProduct(trg->n(i));
		
		// Compute penalty
		pp *= penaltyFunction(mn,0,nrmShift*ref->r());
		pp *= penaltyFunction(mt,0,tngShift*ref->r());
		pp *= penaltyFunction(dn,0,nrmRotn);
	*/
	}
	
	return pp;
}

void AutoPrior::setMovementLimits(double normalShift,double tangentShift,double normalRotation) {
	nrmShift = normalShift;
	tngShift = tangentShift;
	nrmRotn = 1.0 - dcos(normalRotation);
}

void NeighborPriorTerms::compute(MNode *n1,MNode *n2) {
   // Vectors
   v12l = n2->x(MAtom::LEFT) - n1->x(MAtom::LEFT);
   v12r = n2->x(MAtom::RIGHT) - n1->x(MAtom::RIGHT);
   
   // Distances
   d12l = v12l.twoNorm();
   d12r = v12r.twoNorm();
	d12 = n1->x().distanceTo(n2->x());

   // Constraint 1: d12/r2
   c[0] = log2(d12l / n2->r());
	c[1] = log2(d12r / n2->r());
   
   // Constraint 2: r1/r2
   c[2] = log2(n1->r() / n2->r());

   // Constraint 3 nrm(v12) dot n1
   //c[3] = acos(v12l.dotProduct(n1->n()) / d12l);
   //c[4] = acos(v12r.dotProduct(n1->n()) / d12r);
	c[3] = n2->slope(MAtom::LEFT) - n1->slope(MAtom::LEFT);
	c[4] = n2->slope(MAtom::RIGHT) - n1->slope(MAtom::RIGHT);

   // Constraint 4 narrowing rate
   c[5] = acos((n2->r() - n1->r()) / d12 - cos(n1->oa()));

   // Change in object angle
   c[6] = n2->oaDeg() - n1->oaDeg();
}

double NeighborPrior::computePair(NeighborPriorTerms &ref,NeighborPriorTerms &trg) {
   //static Vector sigma1(7,log2(2.0),log2(2.0),log2(2.0),M_PI/12,M_PI/12,M_PI/12,5.0);
   //static Vector sigma2(7,log2(4.0),log2(4.0),log2(4.0),M_PI/6,M_PI/6,M_PI/6,10.0);
	//static Vector sigma1(7,log2(1.1),log2(1.1),log2(1.0),M_PI/12,M_PI/12,M_PI/12,5.0);
   //static Vector sigma2(7,log2(2.0),log2(2.0),log2(2.0),M_PI/6,M_PI/6,M_PI/6,10.0);
	static Vector sigma1(7,log2(1.1),log2(1.1),log2(1.0),M_PI/12,M_PI/12,M_PI/12,5.0);
   static Vector sigma2(7,log2(2.0),log2(2.0),log2(2.0),M_PI/6,M_PI/6,M_PI/6,10.0);

   // Now, get penalties
	double pp = 1;
	
   for(int i=0;i<7 && pp > 0.0;i++) {
      pp *= penaltyFunction(trg.c[i],ref.c[i],sigma1(i),sigma2(i));
   }
	
	return pp;
}

SolutionSpace *NeighborPrior::createSS(MNode *ref) {
   Vector mean,sd;
	
	// Mean comes from reference node
	mean = ref->getData(ref->isEnd());
	
	// Do we care about ends?
   if(ref->isEnd())
      sd = Vector(6,0.01,0.01,0.005,40.0,20.0,0.1);      
   else
		sd = Vector(5,0.01,0.01,0.005,40.0,20.0);
	
	return new GaussianSS(mean,sd);
}

double BoundaryPrior::compute(const BAtom &p,const BAtom &x,const BAtom &s) {
	Vector2D t = -x.tangent();
	
	Vector2D px = x.x-p.x;
	Vector2D sx = s.x-x.x;

	px.normalize();
	sx.normalize();

	Vector2D T = (sx+px)/2;
	
	return t.dotProduct(T);
	/*
	double c1 = -t.dotProduct(s.x - x.x) / s.x.distanceTo(x.x);
	double c2 = -t.dotProduct(x.x - p.x) / p.x.distanceTo(x.x);
	double c3 = -s.tangent().dotProduct(s.x - x.x) / s.x.distanceTo(x.x);
	double c4 = -p.tangent().dotProduct(x.x - p.x) / p.x.distanceTo(x.x);

	c1 = c1 > 0.25 ? 1 : 0;
	c2 = c2 > 0.25 ? 1 : 0;
	c3 = c3 > 0.25 ? 1 : 0;
	c4 = c4 > 0.25 ? 1 : 0;

	return c1*c2*c3*c4;
	*/
}

double BoundaryPrior::priorProbability(MNode *ref,MNode *trg) {
	double pp = 1.0;
	int power = 0;

	for(int i = MAtom::LEFT;i<=MAtom::TAIL;i++) {
		const BAtom &ba = trg->bAtom(i);
		if(ba.bnd) {
			int idx = (int)ba.tShape;

			const BAtom &p = ba.bnd->estimate(ba.bnd->size() + idx - 1);
			const BAtom &n = ba.bnd->estimate(idx + 1);

			pp *= compute(p,ba,n);
			
			power++;
		}
	}

	return pow(pp,1.0/power);
}

SolutionSpace *BoundaryPrior::createSS(MNode *ref) {
   Vector mean,sd;
	
	// Mean comes from reference node
	mean = ref->getData(ref->isEnd());
	
	// Do we care about ends?
   if(ref->isEnd())
      sd = Vector(6,0.01,0.01,0.005,40.0,20.0,0.1);      
   else
		sd = Vector(5,0.01,0.01,0.005,40.0,20.0);
	
	return new GaussianSS(mean,sd);	
}


double NeighborPrior::priorProbability(MNode *ref,MNode *trg) {
   NeighborPriorTerms nptTrg[2];

   double pp = 1;
   if(trg->hasPredecessor()) {
      nptTrg[0].compute(trg->getPredecessor(),trg);
      pp *= computePair(npt[0],nptTrg[0]);
   }
   if(pp > 0.0 && trg->hasSuccessor()) {
      nptTrg[1].compute(trg,trg->getSuccessor());
      pp *= computePair(npt[1],nptTrg[1]);
   }
	
	//
	if(pp > 0.0 && trg->isEnd()) {
		if(trg->endRC() < 1.0) 
			pp *= penaltyFunction(trg->endRC(),1,0.01);
		else
			pp *= penaltyFunction(trg->endRC(),1.2,0.5);
	}

	// Special case for singleton nodes
	if(!trg->hasPredecessor() && !trg->hasSuccessor()) {
		// pp * = computeSingleton(ref,trg);
	}
	
   // Extra terms applied to the penalty to keep interatom spacing similar...
	if(pp > 0.0 && trg->hasPredecessor() && trg->hasSuccessor()) {
		double sigma11,sigma12,sigma21,sigma22;
		if(penalizeSliding) {
			sigma11 = log(1.1);
			sigma12 = log(1.2);
			sigma21 = log(1.1);
			sigma22 = log(1.2);
		} else {
			sigma11 = log(1.6);
			sigma12 = log(3.2);
			sigma21 = log(1.6);
			sigma22 = log(3.2);
		}

      double d12t = trg->getSuccessor()->x().distanceTo(trg->getPredecessor()->x());
		double d12r = ref->getSuccessor()->x().distanceTo(ref->getPredecessor()->x());
	
		pp *= penaltyFunction(log2(nptTrg[0].d12/nptTrg[1].d12),log2(npt[0].d12/npt[1].d12),sigma11,sigma12);
		pp *= penaltyFunction(log2((nptTrg[0].d12 + nptTrg[1].d12)/d12t),
                            log2((npt[0].d12 + npt[1].d12)/d12r),sigma21,sigma22);
	}
	
   return pp;
}

NeighborPrior::NeighborPrior(MNode *reference) {
	penalizeSliding = true;
   if(reference->hasPredecessor()) {
      npt[0].compute(reference->getPredecessor(),reference);
   }
   if(reference->hasSuccessor()) {
      npt[1].compute(reference,reference->getSuccessor());
   }
	
}

NodeProblem::NodeProblem(MNode *ref,MNode *trg,ImageSpace &image,NodePrior &prior) :
m_r(ref),m_t(trg),m_image(image),m_prior(prior) {
	moveAlong = true;
};

double NodeProblem::evaluate(const Vector &v) {
	// Apply vector to the atom
	applyVector(v);
	  	
   // tPenalty.start();
   
   // Compute penalty
	double pp = m_prior.priorProbability(m_r,m_t);
	  
	// tPenalty.stop();

   // If prior probability is very small, return just its log and dont compute the 
	// likelihood
	if(pp <= 0.0)
		return MAX_DOUBLE;
	else if(pp < m_prior.minPP())
		return -log(pp);
	
	// Update counter
	evaluationCost++;
   
   // Compute the medialness   
	// tEval.start();
   double lHood = iMatchComputer.compute(m_t,m_image);
   // tEval.stop();

   return lHood - log(pp);
}

inline int msElapsed(clock_t tStart,clock_t tEnd) {
	return (tEnd-tStart) * 1000 / CLOCKS_PER_SEC;
}

// This will read the settings and determine what node prior to use
NodePrior *createNodePrior(MNode *reference) {
	// Read the name of the prior from the registry
	const char *priorName = regOptions.getStringValue("optimizer.nodes.prior","neighbor");
	
	if(0==strcmp(priorName,"auto")) {
		// Configure the auto prior
		AutoPrior *ap = new AutoPrior;
		ap->setMovementLimits(
			regOptions.getDoubleValue("optimizer.nodes.autoPrior.normalShift",0.5),
			regOptions.getDoubleValue("optimizer.nodes.autoPrior.tangentShift",0.5),
			regOptions.getDoubleValue("optimizer.nodes.autoPrior.normalRotation",0.5));
		return ap;
	}
	else if(0==strcmp(priorName,"boundary")) {
		// Configure the auto prior
		return new BoundaryPrior;
	}
	else {
		NeighborPrior *np = new NeighborPrior(reference);
		np->penalizeSliding = regOptions.getBooleanValue("optimizer.nodes.neighborPrior.penalizeSliding",true);
		return np;
	}
}

// This method uses all the settings in the regOptions registry to create
// an optimization algorithm.  It returns a pointer to the algorithm that should be
// deleted later.  The parameter must be a NumericalFunction because of possibility of a
// ConjugateGradient method being used.  The second parameter is the solution space.
NumericalMethod *createOptAlg(NumericalFunction &problem,SolutionSpace &ss) {
	// Read the name of the prior from the registry
	const char *algName = regOptions.getStringValue("optimizer.algorithm","es");
	
	if(0 == strcmp(algName,"simplex")) {
		return new SimplexMethod(problem,ss);
	}
	else if(0 == strcmp(algName,"conjgrad")) {
		return new ConjugateGradientMethod(problem,ss.getMean());
	}
	else {
		double mu = regOptions.getIntValue("optimizer.es.mu",2);
		double lambda = regOptions.getIntValue("optimizer.es.lambda",4);
		bool plus = regOptions.getBooleanValue("optimizer.es.plus",true);
		int sel = plus ? SELECTION_MuPlusLambda : SELECTION_MuCommaLambda;
		
		EvolutionaryStrategy *es = new EvolutionaryStrategy(problem,ss,mu,lambda,sel);
		es->setSigmaFactor(0.2);
		es->setNth(0,ss.getMean());

		Vector ds(ss.getMean().rows());
		ds.setAll(0.6);
		es->setDeltaSigma(ds);

		return es;
	} 
}

bool optRunning = false;
bool optMustQuit = false;

void optimizeNode(MNode *ref,MNode *trg) {
   // Reset the timers
   tTotal.reset();
   tPenalty.reset();
   tEval.reset();
   tBoundlet.reset();
   
   // Create a prior
	NodePrior *prior = createNodePrior(ref);
	
	// Create a problem
	NodeProblem np(ref,trg,*imageSpace,*prior);

	// Fix along the axis movement if needed
	np.setMoveAlong(regOptions.getBooleanValue("optimizer.nodes.moveAlong",true));

	// Create a solution space for this problem
	SolutionSpace *ss = prior->createSS(ref);
	
	// Create a numerical problem
	NumericalFunction nProblem(np,0.0001);

	// Create an optimization method
	NumericalMethod *method = createOptAlg(nProblem,*ss);
	
	// We will update the display every so often based on elapsed time
	// The optimization will run for a given amount of time as well
	int msUpd = regOptions.getIntValue("optimizer.nodes.icm.msPerUpdate",250);
	int msRun = regOptions.getIntValue("optimizer.nodes.icm.msPerNode",2000);
	
	// Some timers
	clock_t tStart,tUpd,tNow;
	tStart = tUpd = tNow = clock();
	
	// Loop until the time is up or algorithm is finished
	while(!method->isFinished() && msElapsed(tStart,tNow) < msRun && (!optMustQuit)) { 
      // Take a timer lap
      tTotal.nextLap();

		// Run an iteration of the method
		method->performIteration();         
		
		// See if enough time passed for a screen update
		if(msElapsed(tUpd,tNow) > msUpd) {
			// Apply the best ever vector to the node
			np.applyVector(method->getBestEverX());
			
			// Redraw stuff
			onModelOrSelectionChange();
			winGL->redraw();
			
			// Run the event loop
			if (Fl::ready()) {
				Fl::check();
			}
			
			tUpd = tNow;
		}
		
		// Check the time again
		tNow = clock();
	}

   tTotal.stop();

	// Update the model and display
	np.applyVector(method->getBestEverX());
	onModelOrSelectionChange();
	winGL->redraw();

   // Print timing results
   printf("Node timings:\n----------------\n");
   

   printf("Iterations            %d\n",tPenalty.nLaps());
   printf("Calls to eval         %d\n",tEval.nLaps());
   
   printf("Total time:           %lg\n",tTotal.ms());
   printf("Time per iteration    %lg\n",tTotal.ms()/tPenalty.nLaps());
   
   printf("Total penalty time    %lg\n",tPenalty.ms());
   printf("Time in eval          %lg\n",tEval.ms());
   printf("Other time            %lg\n",tTotal.ms() - (tEval.ms() + tPenalty.ms()));

   printf("One evaluation        %lg\n",tEval.ms()/tEval.nLaps());
   
   printf("Boundlet Count        %d\n",tBoundlet.nLaps());
   printf("Total bnd time        %lg\n",tBoundlet.ms());
   printf("Avg bnd time          %lg\n",tBoundlet.ms()/tBoundlet.nLaps());

	delete method;
	delete ss;
	delete prior;
}


// This method will optimize all atoms in a selected set together
void runMultiNodeOpt() {
	// First thing we do is check if we are already optimizing, in which case the call to
	// this function is supposed to stop the optimization
	if(optRunning) {
		optMustQuit = true;
		return;
	}
	else {
		optRunning = true;
		optMustQuit = false;
	}

	// Save the model for undo
	pushModel();

	// First, make a reference model
	MGraph ref = model;

	// Get selections from both models
	MNodeGroup *ngRef = ref.selection();
	MNodeGroup *ngTar = model.selection();

	// List of priors and solution spaces
	vector <NodePrior *>priors;
	vector <SolutionSpace *>spaces;
	int ssSize = 0;

	// Compute the spaces and priors
	for(MNode *n=ngRef->first();n!=NULL;n = ngRef->next()) {
		// Create a prior
		priors.push_back(createNodePrior(n));

		// Get a solution space chunk
		spaces.push_back(priors.back()->createSS(n));
		ssSize += spaces.back()->dimensions();
	}

	// Compute a combined solution space
	Vector ssMean(ssSize),ssSigma(ssSize);
	unsigned int i,j=0;
	for(i=0;i<spaces.size();i++) {
		ssMean.insertMatrix(j,0,spaces[i]->getMean());
		ssSigma.insertMatrix(j,0,spaces[i]->getStandardDeviations());
		j += spaces[i]->getMean().rows();
		delete spaces[i];
	}
	GaussianSS ss(ssMean,ssSigma / spaces.size());
	
	// Create a problem
	MultiNodeProblem np(ngRef,ngTar,*imageSpace,priors);

	// Create a numerical problem
	NumericalFunction nProblem(np,0.0001);

	// Create an optimization method
	NumericalMethod *method = createOptAlg(nProblem,ss);
	
	// We will update the display every so often based on elapsed time
	// The optimization will run for a given amount of time as well
	int msUpd = regOptions.getIntValue("optimizer.nodes.icm.msPerUpdate",250);
	int msRun = regOptions.getIntValue("optimizer.nodes.icm.msPerNode",2000) * 100 * ngRef->size();
	
	// Some timers
	clock_t tStart,tUpd,tNow;
	tStart = tUpd = tNow = clock();
	
	// Loop until the time is up or algorithm is finished
	while(!method->isFinished() && msElapsed(tStart,tNow) < msRun && (!optMustQuit)) { 
		// Run an iteration of the method
		method->performIteration();         
		
		// See if enough time passed for a screen update
		if(msElapsed(tUpd,tNow) > msUpd) {
			// Apply the best ever vector to the node
			np.applyVector(method->getBestEverX());
			
			// Redraw stuff
			onModelOrSelectionChange();
			winGL->redraw();
			
			// Run the event loop
			if (Fl::ready()) {
				Fl::check();
			}

			// Print the best ever solution
			cout << "Best: " << method->getBestEverValue() << "\n";
			cout << "Dist: " << ss.getFeasibility(method->getBestEverX()) << "\n\n";
			
			tUpd = tNow;
		}
		
		// Check the time again
		tNow = clock();
	}

	// Update the model and display
	np.applyVector(method->getBestEverX());
	onModelOrSelectionChange();
	winGL->redraw();

	// Delete the prior
	for(i=0;i<priors.size();i++)
		delete priors[i];
	
	// Delete the method
	delete method;

	optRunning = false;
	
}

void runNodeOptimization() {	
	// First thing we do is check if we are already optimizing, in which case the call to
	// this function is supposed to stop the optimization
	if(optRunning) {
		optMustQuit = true;
		return;
	}
	else {
		optRunning = true;
		optMustQuit = false;
	}

	// Save the model for undo
	pushModel();

	// Now lets get a node group from the model
	MNodeGroup *ng = model.selection();

	// See how many atoms there are here
	if(ng->size() == 1) {
		// Only one node in the group.  This means we don't need to do any permutations.  Just
		// optimize the one node.  What is done here is a bad idea : I make a copy of just this one
		// node.  Copying nodes should probably be disallowed...
		MNode ref = *(ng->node(0));

		// Now, optimize the node and be done
		optimizeNode(&ref,ng->node(0));
	}
	else if(ng->size() > 1) {
		// First, make a reference model
		MGraph ref = model;
		MNodeGroup *ngRef = ref.selection();

		// See how many times we want to do this
		int icmIter = regOptions.getIntValue("optimizer.nodes.icm.iterations",5);

		for(int i=0;i < icmIter;i++) {
         // First, make a reference model
		   // MGraph ref = model;
		   MNodeGroup *ngRef = ref.selection();

			// Create a permutation vector
			vector<int> perm(ng->size());
			for(int i=0;i<ng->size();i++)
				perm[i] = i;
			random_shuffle(perm.begin(),perm.end());

			// Run through this permutation
			for(int j=0;j<ng->size();j++) {
				// Grab a pair of nodes to optimize
				optimizeNode(ngRef->node(perm[j]),ng->node(perm[j]));

				// If we are told to stop, do so
				if(optMustQuit)
					break;
			}

			// If we are told to stop, do so
			if(optMustQuit)
				break;
		}
	}

	optRunning = false;
}


void runRegistration() {
	// First thing we do is check if we are already optimizing, in which case the call to
	// this function is supposed to stop the optimization
	if(optRunning) {
		optMustQuit = true;
		return;
	}
	else {
		optRunning = true;
		optMustQuit = false;
	}

	// Either work with the model or with the selection
	MNodeGroupIF *ng = (model.selection()->size()) ? (MNodeGroupIF*) model.selection() : (MNodeGroupIF*) &model;
	
   // Define the solution space
	double sSpace = regOptions.getDoubleValue("optimizer.rigid.sigmaSpace",0.05);
	double sAngle = regOptions.getDoubleValue("optimizer.rigid.sigmaAngle",45);
	double sScale = regOptions.getDoubleValue("optimizer.rigid.sigmaScale",0.5);
   GaussianSS ss(Vector(4),Vector(4,sSpace,sSpace,sScale,sAngle));
	
	// Create a problem
	RegProblem p(*ng,*imageSpace,ss.getStandardDeviations());
   NumericalFunction np(p,0.00001);

	// Get an optimization algorithm
	NumericalMethod *mth = createOptAlg(np,ss);
   
   // Save the model for undo purposes
   pushModel();
	
	// See how long to run for
	int msTotal = regOptions.getIntValue("optimizer.rigid.msTotal",5000);
   int msUpdate = regOptions.getIntValue("optimizer.rigid.msUpdate",250);
   Timer tUpdate,tTotal;

   // Go register
	int iter = 0;
   double lastBestEver = 0;
   while(!mth->isFinished() && tTotal.ms() < msTotal)  {
      tUpdate.nextLap();
      tTotal.nextLap();

      mth->performIteration();
      if(mth->getBestEverValue() < lastBestEver) {
         printf("Value[%d] = %lg\tCost = %lg\n",iter,
				mth->getBestEverValue(),mth->getFunction().getEvaluationCost());
			
			Vector x = mth->getBestEverX();
         p.applyVector(x);
         
         // Apply model
         onModelOrSelectionChange();
			
         // Redraw model (respond to events?)
         winGL->redraw();         
			
         // printf("I = %d\tF = %lg\n",iter,mth->getBestEverValue());
			
         lastBestEver = mth->getBestEverValue();
      }
		
      if(tUpdate.ms() < msUpdate) {
         tUpdate.reset();

			// Run the event loop
			if (Fl::ready()) {
				Fl::check();
				
				// If the procedure has just been entered, the allow entry flag would be set to true.
				// That means user wants to cancel
				if (optMustQuit) 
					break;
			}
		}

      iter++;
   }
	
  	p.applyVector(mth->getBestEverX());	
   onModelOrSelectionChange();
   winGL->redraw();
   
	delete mth;
	optRunning = false;	
}


/*
void computeEndBlumEstimate(DSLPrimitive2D &p,const DSLPrimitive2D &src,bool fwd) {
   DSLPrimitive2D q;
   double wp = 0.5;
   double wq = 0.5;

   double mult = fwd ? -1 : 1;
   double d = (p.getPosition() - src.getPosition()).twoNorm();
   
   // Position goes along the normal
   q.setPosition(src.getPosition() + src.getNormal(0)*mult*d);

   // TA stays same
   q.setTA(src.getTA());

   // Radius changes by 
   double r = src.getR()-mult*d*dcos(src.getTO());
   q.setR(r < 0 ? 0 : r);      

   // TO stays same
   q.setTO(src.getTO());

   p.setPosition(p.getPosition()*wp + q.getPosition()*wq);
   p.setR(p.getR()*wp + q.getR()*wq);
   p.setTA(p.getTA()*wp + q.getTA()*wq);
   p.setTO(p.getTO()*wp + q.getTO()*wq);
}

void computeMiddleBlumEstimate(DSLPrimitive2D &p,const DSLPrimitive2D &left,const DSLPrimitive2D &right) {
   double wL = 1.0 / 3;
   double wp = 1.0 / 3;
   double wR = 1.0 / 3;

   DSLPrimitive2D pL,pR;

   Vector2D dLR = right.getPosition() - left.getPosition();
   Vector2D d1 = p.getPosition() - left.getPosition();
   Vector2D d2 = right.getPosition() - p.getPosition();

   double DLR = dLR.normalize();
   double D1 = d1.normalize();
   double D2 = d2.normalize();

   double k1 = (D1/(D1+D2)) * DLR;
   double k2 = (D2/(D1+D2)) * DLR;

   // Find a position between the crossings on the normal vectors and line of positions where
   // distance to each primitive is proportional to distance to each primitive from p
   Vector2D p1 = left.getPosition() + left.getNormal(0) * (k1/(left.getNormal(0).dotProduct(dLR)));
   Vector2D p2 = right.getPosition() - right.getNormal(0) * (k2/(right.getNormal(0).dotProduct(dLR)));
   
   pL.setPosition(p1);
   pR.setPosition(p2);

   double z1 = (p.getPosition()-left.getPosition()).twoNorm();
   double z2 = (p.getPosition()-right.getPosition()).twoNorm();

   pL.setTA(left.getTA());
   pR.setTA(right.getTA());
   pL.setTO(left.getTO());
   pR.setTO(right.getTO());

   // Compute a radius estimate for both guys.
   double r1 = left.getR()-z1*dcos(left.getTO());
   double r2 = right.getR()+z2*dcos(right.getTO());
   pL.setR(r1 < 0 ? 0 : r1);
   pR.setR(r2 < 0 ? 0 : r2);

   // New primitive is the weigted sum of the three
   p.setPosition(pL.getPosition()*wL + p.getPosition()*wp + pR.getPosition()*wR);
   p.setR(pR.getR()*wL + p.getR()*wp + pR.getR()*wR);
   p.setTA(pR.getTA()*wL + p.getTA()*wp + pR.getTA()*wR);
   p.setTO(pR.getTO()*wL + p.getTO()*wp + pR.getTO()*wR);
}


// This method replaces every primitive in the model with its 'Blum estimate'
// This will probably be really cheesy
void computeBlumEstimate(DSLObject2D &model) {
   DSLObject2D src = model;

   for(int i=0;i<src.nfigures();i++) {
      DSLFigure2D &f = src.getFigure(i);
      if(f.nprimitives() < 2)
         continue;

      for(int j=0;j<f.nprimitives();j++) {
         if(f.getPrimitive(j).selected()) {
            // Here's the primitive
            DSLPrimitive2D p = f.getPrimitive(j);

            // See if it's the end primitive
            if(j==0) {
               computeEndBlumEstimate(p,f.getPrimitive(j+1),true);
            }
            else if(j==f.nprimitives()-1) {
               computeEndBlumEstimate(p,f.getPrimitive(j-1),false);
            }
            else {
               computeMiddleBlumEstimate(p,f.getPrimitive(j-1),f.getPrimitive(j+1));
            }

            model.getFigure(i).update(j,p);
         }
      }
   }
}

*/






// This method performs medial interpolation on a pair of nodes
MAtom getMedialInterpolation(MNode *s,MNode *e,double t) {
	static double treshold = 0.000001;

	// All atoms whould have same boundary
	dassert(s->bAtom(MAtom::LEFT).bnd = e->bAtom(MAtom::LEFT).bnd);
	dassert(s->bAtom(MAtom::RIGHT).bnd = e->bAtom(MAtom::RIGHT).bnd);
	dassert(s->bAtom(MAtom::LEFT).bnd = s->bAtom(MAtom::RIGHT).bnd);
	dassert(s->bAtom(MAtom::LEFT).tShape < e->bAtom(MAtom::LEFT).tShape);
	dassert(e->bAtom(MAtom::RIGHT).tShape < s->bAtom(MAtom::RIGHT).tShape);
	dassert(t < 1 && t > 0);

	// We will update the display every so often based on elapsed time
	// The optimization will run for a given amount of time as well
	int msRun = 3000,msRestart = 250;
	
	// Create a problem
	MIProblem problem(s,e,t,0);

	// Create a numerical problem
	NumericalFunction nProblem(problem,0.0001);

	// Mehtod used
	NumericalMethod *method = NULL;
	
	// Run the interpolator
	Vector ssMean(3,t,1-t);
	Vector ssSigma(3,0.1,0.1);
	GaussianSS ss(ssMean,ssSigma);
	
	// Some timers
	clock_t tStart,tNow,tRestart;
	tStart = tNow = clock();

	// Start with no help
	double help = 0.0;

	while(msElapsed(tStart,tNow) < msRun) {
		// Restart timer
		tRestart = clock();

		// Create an optimization method
		if(method != NULL) delete method;
		method = createOptAlg(nProblem,ss);
	
		// Loop until the time is up or algorithm is finished
		while(method->getBestEverValue() > treshold && msElapsed(tRestart,tNow) < msRestart) { 
			// Run an iteration of the method
			method->performIteration();         
				
			// Check the time again
			tNow = clock();
		}

		// Break out if everything finished
		if(method->getBestEverValue() < treshold)
			break;

		// Update help
		help += 0.1;
		printf("help = %lg\n",help);
		problem.setHelp(rand(help));
	}

	printf("x = %lg,%lg\tp = %lg\n",method->getBestEverX()(0),method->getBestEverX()(1),method->getBestEverValue());

	// Update the model and display
	MAtom a = problem.getAtom(method->getBestEverX());
	delete method;

	return a;
}


/***********************************************************
 Node groups for m-rep based optimization
 ***********************************************************/


SplineSampleImageMatch::SplineSampleImageMatch(RegularSplineSample *sample,ImageSpace *space) 
: match(sample->samples.size()) 
{
	this->sample = sample;
	this->space = space;

	// This is an addition for disabling r-proportional sampling.  This is a hack.
	scale = space->width() / regOptions.getDoubleValue("optimizer.spline.matchAperturePer256",256);

	for(unsigned int iCurve=0;iCurve<sample->samples.size();iCurve++) {
		match[iCurve].resize(sample->samples[iCurve].size());
		for(unsigned int iSeg=0;iSeg<sample->samples[iCurve].size();iSeg++) {
			integrateSegment(iCurve,iSeg);
		}
	}
}

double SplineSampleImageMatch::integrateImageMatch() {
	double totalLength = 0;
	double totalMatch = 0;

	for(unsigned int iCurve=0;iCurve<sample->samples.size();iCurve++) {
		for(unsigned int iSeg=0;iSeg<sample->samples[iCurve].size();iSeg++) {
			if(match[iCurve][iSeg].ts < sample->samples[iCurve][iSeg].ts) {
				integrateSegment(iCurve,iSeg);
			}
			totalLength += match[iCurve][iSeg].length;
			totalMatch += match[iCurve][iSeg].ival;
		}
	}

	return totalMatch / totalLength;
}

inline double SplineSampleImageMatch::computeMatch(MAtom &a) {
	
	// HACK!
	const BAtom &b1 = a.bAtom(MAtom::LEFT);
	const BAtom &b2 = a.bAtom(MAtom::RIGHT);
	
	return 
		space->applyDOGKernel(b1.x.scaledBy(space->width(),space->height()),b1.n,scale) + 
		space->applyDOGKernel(b2.x.scaledBy(space->width(),space->height()),b2.n,scale);
	
	/*
	return 
		space->applyDOGKernel(b1.x.scaledBy(space->width(),space->height()),b1.n,b1.scale*space->width()) + 
		space->applyDOGKernel(b2.x.scaledBy(space->width(),space->height()),b2.n,b2.scale*space->width());
		*/
}

void SplineSampleImageMatch::integrateSegment(int iCurve, int iSeg) {
	SplineSegment &segment = sample->samples[iCurve][iSeg];	
	list<SplineSamplePoint>::iterator b = segment.points.begin(),a=b++;

	// cout << iSeg << endl;

	// Initialize the integrals
	match[iCurve][iSeg].ival = 0;
	match[iCurve][iSeg].length = 0;


	// Compute the match at first point
	double aMatch = computeMatch(a->atom);

	// Compute along the whole curve
	while(b != segment.points.end()) {
		double bMatch = computeMatch(b->atom);
		double len = a->dn[1] + a->dn[2];

		// Update the integrals
		match[iCurve][iSeg].ival += 0.5 * len * (aMatch + bMatch);
		match[iCurve][iSeg].length += len;

		// Update the counters
		aMatch = bMatch;
		a++;
		b++;
	}
	
	// Update the timestamp
	match[iCurve][iSeg].ts = segment.ts;
}



SNGSegment::SNGSegment(SplineObject *spline,int iCurve,int iStart,int iEnd) 
: SplineNodeGroup(spline) 
{
	this->iCurve = iCurve;
	this->iStart = iStart;
	this->iEnd = iEnd;
	curve = spline->curves[iCurve];

	// Create a vector map
	int iVector = 0;
	for(int i=iStart;i<=iEnd;i++) {
		vmap.push_back(iVector++);
		vmap.push_back(iVector++);
		if(i!=1 && i!=curve->size()-2) {
			// R can't be modified at these points
			vmap.push_back(iVector++);	
		}
		else {
			vmap.push_back(-1);
		}
	}

	// This is how many elements we need
	vSize = iVector;
}

Vector SNGSegment::getVector() {
	Vector v(vSize);

	// Create a vector
	for(unsigned int i=0;i<vmap.size();i+=3) {
		MySMLVec3f c = curve->getControlVec(iStart+i/3);
		for(unsigned int d=0;d<3;d++) {
			if(vmap[i+d] >= 0)
				v(vmap[i+d]) = c.data()[d];
		}
	}

	return v;
}

Vector SNGSegment::getStepSize() {	
	// For step size calculations
	double sSpace = regOptions.getDoubleValue("optimizer.spline.sigmaSpace",0.02);
	double sRadius = regOptions.getDoubleValue("optimizer.spline.sigmaRadius",0.004);

	Vector step(vSize);
	
	for(unsigned int i=0;i<vmap.size();i+=3) {
		step(vmap[i]) = sSpace;
		step(vmap[i+1]) = sSpace;
		if(vmap[i+2] >= 0) {
			step(vmap[i+2]) = sRadius;
		}
	}
	return step;
}


void SNGSegment::applyVector(const Vector &v) {
	// Update using the vector
	for(unsigned int i=0;i<vmap.size();i+=3) {
		MySMLVec3f c = curve->getControlVec(iStart+i/3);
		for(int d=0;d<3;d++) {
			if(vmap[i+d] >= 0)
				c.data()[d] = v(vmap[i+d]);
		}
		curve->updateControl(iStart+i/3,c);
	}
}

void SNGSegment::updateSample(RegularSplineSample *sample) {
	sample->update();
}


SNGBranch::SNGBranch(SplineObject *spline,SplineBranch *branch) 
: SplineNodeGroup(spline) 
{
	this->branch = branch;
}

Vector SNGBranch::getVector() {
	Vector v(9);
	int vi = 0;

	MySMLVec3f C = branch->crv[0]->getControlVec(branch->ei[0]);		
	Vector2D CX(C.x,C.y);
		
	v(0) = C.x;
	v(1) = C.y;
	v(2) = C.z;
	for(int i=0;i<3;i++) {
		MySMLVec3f X = branch->crv[i]->getControlVec(branch->ri[i]);			
		Vector2D DX = Vector2D(X.x,X.y) - CX;
		double len = DX.twoNorm();
		if(len > 0) {
			v(3+i) = DX.getSlope();
			v(6+i) = DX.twoNorm();
		}
		else {
			v(3+i) = 0;
			v(6+i) = 0;
		}
	}
	return v;
}


Vector SNGBranch::getStepSize() {	
	// For step size calculations
	double sSpace = regOptions.getDoubleValue("optimizer.spline.sigmaSpace",0.02);
	double sRadius = regOptions.getDoubleValue("optimizer.spline.sigmaRadius",0.004);

	Vector v(9);	
	v(0) = v(1) = sSpace;
	v(2) = sRadius;
	v(3) = v(4) = v(5) = M_PI / 12;
	v(6) = v(7) = v(8) = sRadius;
	return v;
}

void SNGBranch::applyVector(const Vector &v) {
	// Center posision
	Vector2D CX(v(0),v(1));
	float CR = v(2);
	branch->crv[0]->updateControl(branch->ei[0],CX,CR);
	for(int i=0;i<3;i++) {
		float a = v(3+i);
		float d = v(6+i);
		Vector2D X(CX.x + d * cos(a),CX.y + d * sin(a));
		branch->crv[i]->updateControlX(branch->ri[i],X);
	}
}

void SNGBranch::updateSample(RegularSplineSample *sample) {
	sample->update();
}

/***********************************************************
 Spline Constraints
 ***********************************************************/
/*
class SplineConstraints {
public:
	SplineObject *spline;
	double computeConstraints()
};*/


double SplineMatchProblem::computeConstraints() {
	double cnsBoundary = 0.0;
	double cnsRadius = 0.0;
	double cnsBranch = 0.0;
	double len = 0.0;

	for(unsigned int iCurve=0;iCurve<sample->samples.size();iCurve++) {		
		for(unsigned int iSeg=0;iSeg<sample->samples[iCurve].size();iSeg++) {		

			list<SplineSamplePoint>::iterator it,it1,itFirst,itLast,itEnd;
			itFirst = sample->samples[iCurve][iSeg].points.begin();
			itEnd = itLast = sample->samples[iCurve][iSeg].points.end();
			itLast--;
			
			for(it = itFirst;it!=itLast;it++) {							
				if(it->cv[0] < 0) {
					cnsBoundary -= it->cv[0] / it->dn[0];
				}
				if(it->cv[1] < 0) {
					cnsBoundary -= it->cv[1] / it->dn[0];
				}

				len += it->dn[1];
				len += it->dn[2];
			}
			for(it=itFirst;it!=itEnd;it++) {
				if(it->atom.r() <= 0) {
					cnsRadius -= it->atom.r();
				}
			}
		}
	}

	for(unsigned int iBranch=0;iBranch<sample->spline->branches.size();iBranch++) {
		cnsBranch += max(0,sample->spline->branches[iBranch]->getMaxAngle());
	}

	return wgtPenaltyJacobian * cnsBoundary/len + 
		wgtPenaltyRadius * cnsRadius + 
		wgtPenaltyBranch * cnsBranch;
}

/**
 * For now, compute the match without any prior terms
 */
double SplineMatchProblem::evaluate(const Vector &x) {
	// Check the vector
	Vector v = x;
	for(int i=0;i<v.size();i++) {
		if(v(i)==std::numeric_limits<float>::quiet_NaN()) {
			cout << "inf" << endl;
			v(i) = 0;
		}
	}

	// Update the spline with the vector
	updateSpline(v);	

	// Create a sample
	// SplineSample *sample = new SplineSample(sng->spline,40);
	// RegularSplineSample *sample = new ArcSplineSample(sng->spline,40);

	// Compute penalty
	double penalty = computeConstraints();

	// Compute the prior
	double prior = computeRegularizationPrior();

	// Compute the image match
	double likehood = imatch->integrateImageMatch();

	// Print out the penalty, lhd
	// cout << "LHD: " << likehood << endl;
	// cout << "PEN: " << penalty << endl;

	// Return result
	return penalty + prior + likehood;
}
double SplineMatchProblem::computeRegularizationPrior() {
	// The prior
	double priorElastic = 0;
	double priorCurvature = 0;
	double meanAngle = 0;
	double length = 0;
	
	// Compute the markov regularization prior
	for(unsigned int iCurve=0;iCurve<sample->spline->curves.size();iCurve++) {	
	
		double p1 = 0,p2 = 0,arcLength = 0;
		int np = 0;

		SplineCurve *curve = sample->spline->curves[iCurve];

		// This is a new regularization prior that checks ds/dt at regular intervals 
		for(unsigned int iSeg=0;iSeg<sample->samples[iCurve].size();iSeg++) {
			
			list<SplineSamplePoint>::iterator it,itEnd;
			itEnd = sample->samples[iCurve][iSeg].points.end();
			it = sample->samples[iCurve][iSeg].points.begin();
			
			while(it!=itEnd) {
				double al = sqrt(it->M[1].x * it->M[1].x + it->M[1].y * it->M[1].y);
				double scale = 1.0/(al*al*al);
				double kappa = (it->M[1].x * it->M[2].y - it->M[1].y * it->M[2].y) * scale;
				double sigma = (it->M[1].x * it->M[2].x + it->M[1].y * it->M[2].y) * scale;

				p1 += sigma * sigma; // it->M[1].x * it->M[2].x + it->M[1].y * it->M[2].y;
				p2 += kappa * kappa;
				arcLength += al;
				np++;
				
				++it;
			}
		}

		// cout << "p1: " << p1 << endl;
		// cout << "p2: " << p2 << endl;
		// cout << "al: " << arcLength << endl;
		// priorElastic += p1 / (arcLength * arcLength);
		priorElastic += p1 * arcLength / (np * np);
		priorCurvature += p2 * arcLength / (np * np);
	 
		/*
		// The elastic prior isn't computed for next-to-terminal control points,
		// because of their special relationship with the last point
		for(iPoint=1;iPoint<curve->size()-2;iPoint++) {
			int i0 = (iPoint==1) ? 0 : iPoint;
			int i1 = (iPoint==curve->size()-3) ? curve->size()-1 : iPoint+1;

			MySMLVec3f X0 = curve->getControlVec(i0);
			MySMLVec3f X1 = curve->getControlVec(i1);
			MySMLVec3f Xi = X1-X0;
			double l2 = Xi.x * Xi.x + Xi.y * Xi.y;
			
			length += sqrt(l2);
			
			priorElastic += l2;
		}
		for(iPoint=1;iPoint<curve->size()-1;iPoint++) {
			MySMLVec3f X0 = curve->getControlVec(iPoint-1);
			MySMLVec3f X1 = curve->getControlVec(iPoint);
			MySMLVec3f X2 = curve->getControlVec(iPoint+1);
			MySMLVec3f Xii = (X2+X0) - (X1+X1);

			// Mean angle computation
			Vector2D dx1(X2.x-X1.x,X2.y-X1.y);
			Vector2D dx2(X0.x-X1.x,X0.y-X1.y);

			meanAngle += 1 + dx1.dotProduct(dx2) / (dx1.twoNorm() * dx2.twoNorm());

			priorCurvature += Xii.x * Xii.x + Xii.y * Xii.y;
			np++;
		}	
		*/
	}

	// meanAngle /= np;
	// We compute the curvature prior divided by the squared total length of the object.
	// This ensures scale invariance: double the object, get same prior.
	// cout << priorCurvature / (length * length) << endl;

	// cout << "Curvature Prior " << priorCurvature <<  endl;
	// cout << "Elastic Prior " << priorElastic <<  endl;
	
	// return (wgtPriorCurvature * priorCurvature + wgtPriorElastic * priorElastic) / (length*length);
	return wgtPriorCurvature * priorCurvature + wgtPriorElastic * priorElastic;
}

SplineMatchProblem::SplineMatchProblem(SplineNodeGroup *sng,ImageSpace *space)
{	
	this->sng = sng;
	this->space = space;
	
	sample = new ArcSplineSample(sng->spline,10);
	imatch = new SplineSampleImageMatch(sample,space);
	cout << "Made new match" << endl;

	// Get the weights from the registry
	wgtPenaltyJacobian	= regOptions.getDoubleValue("optimizer.spline.weight.jacobianPenalty",1.0e3);
	wgtPenaltyRadius	= regOptions.getDoubleValue("optimizer.spline.weight.radiusPenalty",1.0e9);
	wgtPenaltyBranch	= regOptions.getDoubleValue("optimizer.spline.weight.branchPenalty",1.0e2);
	wgtPriorCurvature	= regOptions.getDoubleValue("optimizer.spline.weight.regularizationPrior",1.0e2);
	wgtPriorElastic		= regOptions.getDoubleValue("optimizer.spline.weight.elasticPrior",1.0e2);
}

SplineMatchProblem::~SplineMatchProblem() {
	delete sample;
	delete imatch;
}

GaussianSS SplineMatchProblem::getSolutionSpace()
{
	Vector vMean = sng->getVector();
	Vector vSigma = sng->getStepSize();
	return GaussianSS(vMean,vSigma);
}


/**
 * This method will match a spline to an image, k control points at a time
 */
void runSplineMatch(SplineObject *mrep, ImageSpace *space) {
	unsigned int iGroup;

	// First thing we do is check if we are already optimizing, in which case the call to
	// this function is supposed to stop the optimization
	if(optRunning) {
		optMustQuit = true;
		return;
	}
	else {
		optRunning = true;
		optMustQuit = false;
	}

	// Get all the parameters from the registry
	int nHoodSize = regOptions.getIntValue("optimizer.spline.problemSize",2);
	int nIterations = regOptions.getIntValue("optimizer.spline.iterations",10);
	int msPerItn = regOptions.getIntValue("optimizer.spline.msPerIteration",500);
	int msUpdate = regOptions.getIntValue("optimizer.spline.msUpdate",250);
	int cgScaleFactor = regOptions.getIntValue("optimizer.conjgrad.gradMultiplier",0.001);

	for(int iIteration=0;iIteration<nIterations;iIteration++) {
		bool madeImprovement = false;

		// Create an array of all possible node groups
		vector<SplineNodeGroup*> groups;
		for(int iCurve=0;iCurve<(int)mrep->curves.size();iCurve++) {
			SplineCurve *curve = mrep->curves[iCurve];
			for(int iStart=0;iStart <= curve->size()-nHoodSize;iStart++) {
				groups.push_back(new SNGSegment(mrep,iCurve,iStart,iStart+nHoodSize-1));
			}
		}
		for(int iBranch=0;iBranch<(int)mrep->branches.size();iBranch++) {
			SNGBranch *br =  new SNGBranch(mrep,mrep->branches[iBranch]);
			Vector v = br->getVector();
			br->applyVector(v);
			Vector v1 = br->getVector();
			groups.push_back(br);
			// groups.push_back(new SNGBranch(mrep,mrep->branches[iBranch]));
		}

		// Create a permutation array for this iteration
		random_shuffle(groups.begin(),groups.end());

		// Optimize each neighborhood
		for(iGroup=0;iGroup < groups.size();iGroup++) {
			
			// Create a problem
			SplineMatchProblem p(groups[iGroup],space);

			// Create a gaussian solution space around the problem
			GaussianSS ss = p.getSolutionSpace();

			// Create an optimization method 
			NumericalFunction np(p,0.0001);
			np.setGradientScaleFactor(0.001);
			NumericalMethod *mth = createOptAlg(np,ss);
			
			cout << "Starting Best: " << mth->getBestEverValue() << endl;
   
			// Create update and iteration timers
			Timer tUpdate,tIteration;

			// Go register
			double fLastBestEver = mth->getBestEverValue();
			double fBefore = fLastBestEver;

			// cout << "START = " << fLastBestEver;
   
			while(!mth->isFinished() && tIteration.ms() < msPerItn)  {
				tUpdate.nextLap();
				tIteration.nextLap();

				// Perform an iteration
				mth->performIteration();				

				// Show best ever values
				if(mth->isFinished()) {
					cout << "*";
				}
				else if(mth->getBestEverValue() < fLastBestEver) {
					cout << "+"; // << /* mth->getBestEverValue()*/ << " ";
					fLastBestEver = mth->getBestEverValue();
				}
				else {
					cout << ".";
				}

				// Check if it is time to update the display
				if(tUpdate.ms() > msUpdate || mth->isFinished()) {
					// Display the best ever model
					p.updateSpline(mth->getBestEverX());

					cout << ":";

					// Update the screen display
					onModelOrSelectionChange();
					winGL->redraw();

					// Run the event loop
					if (Fl::ready()) {
						Fl::check();
					
						if (optMustQuit) break;
					}

					tUpdate.reset();
				}
			} // while method is still running

			// Apply the best result
			p.updateSpline(mth->getBestEverX());

			cout << endl << "Final Best: " << mth->getBestEverValue() << endl;

			if(fBefore > mth->getBestEverValue())
				madeImprovement = true;

			// Delete the method
			delete mth;

			if (optMustQuit) break;
		} // for each control point set

		// cout << endl;

		// Delete the groups
		for(iGroup=0;iGroup < groups.size();iGroup++) {
			delete groups[iGroup];
		}

		// No improvement in an iteration
		if(!madeImprovement)
			break;

		if (optMustQuit) break;
	} // for each iteration
	
	// Update the display
	onModelOrSelectionChange();
	winGL->redraw();

	// Optimization is no longer running
	optRunning = false;	
}

