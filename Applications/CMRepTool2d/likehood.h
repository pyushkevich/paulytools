#ifndef _LIKELEHOOD_
#define _LIKELEHOOD_

/*****************************************************
This module contains code for computing likelihood
of a model, figure or primitive

We use DSL2D,ImageSpace and Crep2D here
****************************************************/
#include <mreps2D.h>
#include "ispace.h"
#include "SplineMrep.h"

class SplineObject;
class RegularSplineSample;

// This structure describes the medialness that is used
class IMatchComputer {
private:
	int convertFilter(const char *filterName);
	int fixedScale;

public:

	static const int FILTER_NONE;
	static const int FILTER_POWER;
	static const int FILTER_LPP;
	static const int FILTER_LAPLACEAN;
	static const int FILTER_DOG;
	static const int FILTER_DOG_FIXED_SCALE;
	static const int FILTER_DEFAULT;

	// Type of filter that is applied at each boundary primitive
	int dscrFilter;

	// Type of filter that is applied between boundary primtives
	int contFilter;

	// This parameter determines whether the kernel is continuous or discrete.
	// When set to 0, only discrete filters are used, if set to 1 only continuos
	// filters are used
	double contWeight;

	// Samples per boundlet
	double contSamples;

	// When power filter is used, this determines it's paramter
	double powerKernelParameter;

	IMatchComputer() {
		dscrFilter = FILTER_NONE;
		contFilter = FILTER_NONE;
		contWeight = 0.0;
		powerKernelParameter = 4.0;
		contSamples = 5.0;
	}

	// Read the values of all parameters from a registry folder
	void readRegistry(Registry &folder);

	double compute(const BAtom &ba,ImageSpace &space,int filter=FILTER_DEFAULT);
	double compute(MNode *node,ImageSpace &space);
	double compute(MNodeGroupIF &ng,ImageSpace &space);

	// double compute(SplineObject &mrep,ImageSpace &space);
	// double compute(RegularSplineSample &sample,ImageSpace &space);
};

class SplineSampleImageMatch {
protected:
	RegularSplineSample *sample;
	ImageSpace *space;

	struct Segment {
		double ival;
		double length;
		Timestamp ts;
	};

	vector< vector <Segment> > match;

	void integrateSegment(int iCurve,int iSeg);

	// This is a hack
	double scale;

	double computeMatch(MAtom &a);

public:
	SplineSampleImageMatch(RegularSplineSample *sample,ImageSpace *space);

	double integrateImageMatch();
};


#endif

