// SplineMrep.h: interface for the SplineMrep class.
//
//////////////////////////////////////////////////////////////////////
#ifndef SplineMRep_H
#define SplineMRep_H

#include <SMLVec4f.h>
#include <vector>
#include <vector2d.h>
#include <mreps2D.h>
#include <bspline.h>

using namespace std;

class ImageSpace;
class IMatchComputer;
class SplineBranch;
class SplineObject;

/*
class ControlRelationship {
public:
	void onChange(SplineMRep *mrep,int control)
private:

};
*/

class Timestamp {
private:
	long ts;
	static long t;
public:

	Timestamp() {
		update();
	}
/*
	Timestamp(long l) {
		ts = l;
	}
*/
	operator long &() {
		return ts;
	}

	void update() {
		ts = ++t;
	}
};



class SplineCurve : public BSpline1D {
protected:
	// Rho
	double rho;

	// Two slots that can be occupied by branches
	SplineBranch *br[2];

	// Private interpolation method
	void interpolateMedialAtom(const MySMLVec3f &M0,const MySMLVec3f &M1,MAtom &mAtom);

public:
	// Load from a folder
	SplineCurve(Registry *folder);
	virtual ~SplineCurve();

	// Save to a folder
	void save(Registry *folder);

	// An enum for the ends
	enum {TAIL = 0,HEAD = 1};

	// Get the spline position at a parameter value between 1 and P.size-2 and 
	void interpolateAxis(int seg,float t,MySMLVec3f &target) {		
		MySMLVec4f W[4];
		kv.basisJet(seg+3,t,W);
		interpolatePoint(seg,W,0,0,2,target.data());
	}

	// This method sets a control point without checking any constraints
	void setControlVec(int idx,const MySMLVec3f &x) {
		BSpline1D::setControl(idx,0,x.x);
		BSpline1D::setControl(idx,1,x.y);
		BSpline1D::setControl(idx,2,x.z);

		// Refresh timestamp
		tsPoint[idx].update();
	}

	// Set the position of a control point, and apply all branch/end constraints
	void updateControl(int idx, const MySMLVec3f &x);

	// Set the position of a control point, and apply all branch/end constraints
	void updateControl(int idx, const Vector2D &x, float r) {
		updateControl(idx,MySMLVec3f(x.x,x.y,r));
	}

	// Set the position of a control point
	void updateControlX(int idx, const Vector2D &x) {
		updateControl(idx,x,BSpline1D::getControl(idx,2));
	}

	// Set the r of a control point
	void updateControlR(int idx, float r) {
		updateControl(idx,Vector2D(BSpline1D::getControl(idx,0),BSpline1D::getControl(idx,1)),r);
	}

	// Enforce boundary or branch constraints
	void enforceConstraints(bool first=true,bool last=true);

	// Compute the medial atom
	void interpolateAtom(int seg,float t,MAtom &mAtom);

	// Compute the medial atom as well as derivatives of the spline curve
	void interpolateAtomJet(int seg,float t,MAtom &mAtom,MySMLVec3f &x,MySMLVec3f &xt,MySMLVec3f &xtt);

	// Interpolate the end cap boundary
	// void interpolateEndCap(int end,float t,BAtom &ba) const;

	MySMLVec3f getControlVec(int i) {
		return MySMLVec3f(BSpline1D::getControl(i,0),BSpline1D::getControl(i,1),BSpline1D::getControl(i,2));
	}

	// Size of what???
	int size() const {
		return m+1;
	}

	void setRho(double rho) {
		this->rho = rho;
		for(int i=0;i<=m;i++)
			tsPoint[i].update();
	}

	SplineCurve(int nPoints,double rho=0.25,double rDef=0.05);


	// Timestamp array for the control points
	Timestamp *tsPoint;

	friend class SplineObject;
	friend class SplineBranch;
	friend class RegularSplineSample;
};

class SplineBranch {
public:
	SplineCurve *crv[3];
	int ei[3];
	int ri[3];

	void enforceConstraints();

	double getMaxAngle();
};

class SplineObject {
public:
	// A collection of spline curves
	vector<SplineCurve*> curves;

	// A collection of branches
	vector<SplineBranch*> branches;

	// Is it possible to add a branch at a point?
	bool canAddBranch(SplineCurve *curve,int idx);
	
	// Add a branch to the collection (this adds a generic subfigure to an existing figure)
	SplineCurve *addBranch(SplineCurve *curve,int idx,int nPoints);

	// Add a branch to the collection (this adds a generic subfigure to an existing figure)
	SplineCurve *addBranch(SplineCurve *curve,int idx,SplineCurve *sub,int subPoint);

	// Add a non-connected curve
	SplineCurve *addCurve(int nPoints, float rho=0.25);
	SplineCurve *addCurve(SplineCurve *curve);

	// Remove a curve
	void removeCurve(int index);

	// Remove a branch
	void removeBranch(SplineBranch *branch);

	// Save to a registry
	void save(Registry &folder);

	// Load from registry
	void load(Registry &folder);

	// Overall timestamp
	Timestamp tsStructure;
	
private:

};

/*
 A sample point from a spline
 */
struct SplineSamplePoint {
	// The medial axis and its derivatives
	MySMLVec3f M[3];

	// The medial atom
	MAtom atom;

	// The segment and the t value
	double t;
	int seg;

	// Constraint violation?
	bool bad;

	// Squared distance to the next point along medial axis and boundary
	double dn[3];

	// Constraint violation terms for each interval
	double cv[2];

	// Image matches 
	double im[2];

	SplineSamplePoint(double t,int seg) {
		this->t = t;
		this->seg = seg;
		bad = false;
	}
};

class SplineSegment {
public:
	list<SplineSamplePoint> points;
	Timestamp ts;
};

class RegularSplineSample {
protected:
	// Initialize
	void init(int nPerSegment);

	// Update the samples associated with a stretch of the curve
	virtual void updateSegment(int iCurve,int iSeg);

	// Compute the timestamp affecting a segment
	long getTimestampForSegment(int iCurve,int iSegment);

	// Protected constructor
	RegularSplineSample(SplineObject *object);

	// Compute distances and related sample info
	void computeDistances(SplineSamplePoint &a,SplineSamplePoint &b);

public:
	SplineObject *spline;
	vector< vector< SplineSegment > > samples;

	RegularSplineSample(SplineObject *object,int nPerCurve);

	// Time stamp for the whole object
	Timestamp tsOverall;

	// Update the curves as necessary (assumes no structural changes)
	void update();
};

class ArcSplineSample : public RegularSplineSample {
protected:
	int nPerSegmentMin;
		
	void updateSegment(int iCurve,int iSeg);

public:
	ArcSplineSample(SplineObject *object,int nPerSegmentMin);		
};


#endif // SplineMRep_H
