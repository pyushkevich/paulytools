#include "mspline.h"


/*
class CompleteTrimCurve {
private:
	vector<float> u,v;
	
	// Patch sort function
	static int sortFunction(void *,void *);

public:
	void compute(SplineDataCache *sd);
};



void CompleteTrimCurve::compute(SplineDataCache *sd) {
	struct Link {
		list<MedialPoint *> mp;
	};

	// Construct the links
	list<Link> links;
	for(int i=0;i<sd->patch.size();i++) {
		PatchDataCache *pdc = sd->patch(i);
		int t0 = 0;
		for(int cv=0;cv<pdc->idxTrimCurvesIdx;cv++) {
			Link l;
			for(int t=t0;t<pdc->idxTrimCurves[cv];t++) {				
				l.mp.push_back(pdc->aux(t));
			}
			t0 = t;
			links.push_back(l);
		}
	}
	
	// Connect the links
	while(links.size() > 2) {
		Link l = links.pop_front();
		list<Links>::iterator it;
		for(it = links.begin();it!=links.end();it++) {
			Link *l = *it;
		}
	}
};

class CurveAxisProblem : public Function {
public:
	enum {m=10};

	CurveAxisProblem(MSpline *spline,SplineDataCache *splineData);
	~CurveAxisProblem();

	double evaluate(const Vector &v);
	void applyVector(const Vector &v);

	Vector getCurrentState();

private:
	MSpline *spline;
	SplineDataCache *splineData;

	DynamicBSpline2D *bsp;	

	// Window definition
	int w0,w1;
};

CurveAxisProblem::CurveAxisProblem(MSpline *spline,SplineDataCache *splineData) {
	this->spline = spline;
	this->splineData = splineData;	

	bsp = new DynamicBSpline2D(m,3,4,4);
	
	// Initialize the spline to a basic configuration
	for(int i=0;i<=m;i++) {
		bsp->setControl(i,0,0,
	}
	
}

*/