#include "likehood.h"
#include "SplineMRep.h"

/*****************************************************************************
Compute medialness of a model
*****************************************************************************/
const int IMatchComputer::FILTER_POWER = 1;
const int IMatchComputer::FILTER_DOG = 2;
const int IMatchComputer::FILTER_LPP = 3;
const int IMatchComputer::FILTER_LAPLACEAN = 4;
const int IMatchComputer::FILTER_DOG_FIXED_SCALE=5;
const int IMatchComputer::FILTER_NONE = 0;
const int IMatchComputer::FILTER_DEFAULT = -1;

double IMatchComputer::compute(const BAtom &ba,ImageSpace &space,int filter) {
	// Which filter?
	if(filter == FILTER_DEFAULT)
		filter = dscrFilter;

	// Space size
	double w = space.width();
   double h = space.height();

	// Compute at this atom
	Vector2D x = ba.x.scaledBy(w,h);
	double s = ba.scale * w;

	double result = 0;

	switch(filter) {
	case FILTER_POWER : 
      result = -fabs(space.applyPowerKernel(x,ba.n,s,powerKernelParameter));
		break;

   case FILTER_DOG :
      // result = -fabs(space.applyDOGKernel(x,ba.n,s)); 
	   result = space.applyDOGKernel(x,ba.n,s); 
		break;

   case FILTER_LPP : 
      result = -fabs(space.applyLPPKernel(x,ba.n,s));
		break;

   case FILTER_LAPLACEAN :
      result = -fabs(space.applyLaplaceanKernel(x,s));
		break;
	}

	return ba.polarity * result;
}

int IMatchComputer::convertFilter(const char *filterName) {
	if(0 == strcmp(filterName,"none"))
		return FILTER_NONE;
	else if(0 == strcmp(filterName,"dog"))
		return FILTER_DOG;
	else if(0 == strcmp(filterName,"power"))
		return FILTER_POWER;
	else if(0 == strcmp(filterName,"lpp"))
		return FILTER_LPP;
	else if(0 == strcmp(filterName,"laplacean"))
		return FILTER_LAPLACEAN;
	return FILTER_NONE;
}

void IMatchComputer::readRegistry(Registry &folder) {
	contFilter = convertFilter(folder.getStringValue("continuous.kernel","power"));
	dscrFilter = convertFilter(folder.getStringValue("discrete.kernel","power"));
	contWeight = folder.getDoubleValue("continuous.kernelWeight",1.0);
	contSamples = folder.getIntValue("continuous.sampling",10);
	powerKernelParameter = folder.getDoubleValue("powerParameter",4.0);
}

double IMatchComputer::compute(MNode *node,ImageSpace &space) {
	double bdisc = 0;
	double bcont = 0;
	double tStep = 1.0 / (contSamples + 1.0);

	for(int i=MAtom::LEFT;i<=MAtom::TAIL;i++) {
		BAtom ba = node->bAtom(i);
		if(ba.bnd != NULL) {
			
			if(dscrFilter != FILTER_NONE) { 
				// Have a usable atom.  Get discrete boundariness
				bdisc += compute(ba,space);
			}
			
			if(contFilter != FILTER_NONE && contWeight > 0) {
				// Now get continuous bness				
				for(double t=tStep;t < 1.0;t+=tStep) {
					BAtom bPred = ba.bnd->estimate(ba.tShape - t);
					BAtom bSucc = ba.bnd->estimate(ba.tShape + t);
				   bcont += compute(bPred,space);
					bcont += compute(bSucc,space);
				}
			}
		}
	}

	return bdisc + contWeight*bcont;
}

double IMatchComputer::compute(MNodeGroupIF &ng,ImageSpace &space) {
	double b = 0;

	for(MNode *node=ng.first();node!=NULL;node=ng.next()) {
		b += compute(node,space);
	}

	return b;
}
/*
double IMatchComputer::compute(RegularSplineSample &sample,ImageSpace &space) {
	double integral = 0;
	double length = 0;

	// Repeat for each medial curve in the object
	for(unsigned int iCurve=0;iCurve<sample.samples.size();iCurve++) {
		
		MAtom *ma0 = &sample.samples[iCurve].front().atom;
		double f0l = compute(ma0->bAtom(MAtom::LEFT),space);
		double f0r = compute(ma0->bAtom(MAtom::RIGHT),space);
		
		list<SplineSamplePoint>::iterator it;
		for(it = sample.samples[iCurve].begin();it!=sample.samples[iCurve].end();it++) {
			
			MAtom *ma1 = &it->atom;
			double f1l = compute(ma1->bAtom(MAtom::LEFT),space);
			double f1r = compute(ma1->bAtom(MAtom::RIGHT),space);
		
			double ll = ma0->x(MAtom::LEFT).distanceTo(ma1->x(MAtom::LEFT));
			double lr = ma0->x(MAtom::RIGHT).distanceTo(ma1->x(MAtom::RIGHT));

			// Add the result for this segment
			integral += 0.5 * (f1l+f0l) * ll;
			integral += 0.5 * (f1r+f0r) * lr;
			length += lr + ll;
					
			ma0 = ma1;
			f0r = f1r;
			f0l = f1l;
		}
	}

	return integral / length;
}
*/