// SplineMrep.cpp: implementation of the SplineMrep class.
//
//////////////////////////////////////////////////////////////////////
#include "SplineMrep.h"
#include <algorithm>
#include <limits>

long Timestamp::t = 0;

void SplineCurve::interpolateMedialAtom(const MySMLVec3f &M0,const MySMLVec3f &M1,MAtom &mAtom) {
	// Position of the atom
	mAtom.x(Vector2D(M0.x,M0.y));
	mAtom.r(M0.z);
	mAtom.rho(rho);
	
	// The object angle and the axis angle
	double ds_dt = sqrt(M1.x*M1.x + M1.y*M1.y);
	
	// Because of numerical error, handle this carefully
	/*
	if(fabs(ds_dt-M1.z) < 0.00001) {
		mAtom.oa(acos(-1));
	}
	else if(fabs(ds_dt+M1.z) < 0.00001) {
		mAtom.oa(acos(1));
	}
	else {
		double dr_ds = M1.z / ds_dt;	
		mAtom.oa(acos(-dr_ds));
	}*/

	double dr_ds = M1.z / ds_dt;	
	if(dr_ds > 1)
		dr_ds = 1;
	if(dr_ds < -1)
		dr_ds = -1;
	mAtom.oa(acos(-dr_ds));
	mAtom.fa(atan2(-M1.y,M1.x));
}

// Compute the medial atom
void SplineCurve::interpolateAtom(int seg,float t,MAtom &mAtom) {

	// Get the basis
	MySMLVec4f W[4];
	kv.basisJet(seg+3,t,W);

	// Compute the zeroth and first derivatives
	MySMLVec3f M0,M1;
	interpolatePoint(seg,W,0,0,2,M0.data());
	interpolatePoint(seg,W,1,0,2,M1.data());
	
	interpolateMedialAtom(M0,M1,mAtom);
}

// Compute the medial atom
void SplineCurve::interpolateAtomJet(int seg,float t,MAtom &mAtom,MySMLVec3f &M0,MySMLVec3f &M1,MySMLVec3f &M2) {

	// Get the basis
	MySMLVec4f W[4];
	kv.basisJet(seg+3,t,W);

	// Compute the zeroth and first derivatives
	interpolatePoint(seg,W,0,0,2,M0.data());
	interpolatePoint(seg,W,1,0,2,M1.data());
	interpolatePoint(seg,W,2,0,2,M2.data());
	
	interpolateMedialAtom(M0,M1,mAtom);
}

/*
// Interpolate the end cap boundary
void SplineMRep::interpolateEndCap(int end,float t,BAtom &ba) const {
	// For now the end cap is just a circle
	MAtom ma;
	
	float a0,a1;
	
	// Get the medial atom at the first/last knot
	if(end == TAIL) {
		interpolateAtom(0,0,ma);
		a0 = ma.fa() + ma.oa();
		a1 = ma.fa() - ma.oa() + 2*M_PI;
	}
	else {
		interpolateAtom(n-4,1,ma);
		a0 = ma.fa() - ma.oa();
		a1 = ma.fa() + ma.oa();
	}
	
	// Interpolate the arc
	float a = a1*t + a0*(1-t);
	
	ba.n.x = cos(a);
	ba.n.y = -sin(a);
	ba.x = ma.x() + ba.n * ma.r();
	ba.scale = ma.r() * rho;
}
*/

SplineCurve::SplineCurve(int nPoints,double rho,double rDef)
: BSpline1D(nPoints-1,3,4)
{
	this->rho = rho;

	// Initialize timestamps
	tsPoint = new Timestamp[nPoints];	

	// Create a stock model
	for(int i=0;i<=m;i++) {
		setControlVec(i,MySMLVec3f(0.1+i*0.8/m,0.5,rDef));
	}

	// Set branches to nulls;
	br[0] = br[1] = NULL;

	// Enforce the end constraints
	enforceConstraints();
}


void SplineCurve::save(Registry *folder)
{
	folder->setIntValue("m",m);
	folder->setDoubleValue("rho",rho);
	for(int i=0;i<=m;i++) {
		MySMLVec3f P = getControlVec(i);
		folder->setDoubleValue("x[%d]",P.x,i);
		folder->setDoubleValue("y[%d]",P.y,i);
		folder->setDoubleValue("r[%d]",P.z,i);
	}
}

SplineCurve::SplineCurve(Registry *folder) :
BSpline1D(folder->getIntValue("m",0),3,4)
{
	rho = folder->getDoubleValue("rho",0.25);
	for(int i=0;i<=m;i++) {
		MySMLVec3f v;
		v.x = folder->getDoubleValue("x[%d]",0,i);
		v.y = folder->getDoubleValue("y[%d]",0,i);
		v.z = folder->getDoubleValue("r[%d]",0,i);
		setControlVec(i,v);
	}
	enforceConstraints();
}

SplineCurve::~SplineCurve() {
	delete tsPoint;
}

void SplineCurve::updateControl(int idx,const MySMLVec3f &x) {
	// Points next to ends never change their r value
	MySMLVec3f v = x;
	if(idx==1 || idx==m-1)
		v.z = BSpline1D::getControl(idx,2);

	// Check the branches
	if(idx==0 && br[0]) {
		for(int i=0;i<3;i++)
			br[0]->crv[i]->setControlVec(br[0]->ei[i],v);		
	}
	else if(idx==m && br[1]) {
		for(int i=0;i<3;i++)
			br[1]->crv[i]->setControlVec(br[1]->ei[i],v);		
	}
	else {
		setControlVec(idx,v);						
	}

	enforceConstraints(idx <= 1,idx >=m-1);

	if(br[0])
		br[0]->enforceConstraints();
	if(br[1])
		br[1]->enforceConstraints();
}

void SplineCurve::enforceConstraints(bool first,bool last) {	
	// Compute the dependencies
	MySMLVec4f W0[4],W1[4];
	MySMLVec3f M0,M1;
	
	

	// Compute the first end atom
	if(first) {
		kv.basisJet(3,0.0f,W0);
		interpolatePoint(0,W0,0,0,2,M0.data());
		interpolatePoint(0,W0,1,0,2,M1.data());
		
		float a = W0[1].x;
		float b = W0[1].y;
		float ds_dt = sqrt(M1.x*M1.x + M1.y*M1.y);
		float r1 = (ds_dt-a*M0.z) / b;

		BSpline1D::setControl(1,2,r1);
		tsPoint[1].update();
	}

	if(last) {
		kv.basisJet(m,1.0f,W1);	
		interpolatePoint(m-3,W1,0,0,2,M0.data());
		interpolatePoint(m-3,W1,1,0,2,M1.data());
		float a = W1[1].w;
		float b = W1[1].z;
		float ds_dt = sqrt(M1.x*M1.x + M1.y*M1.y);
		float r1 = -(ds_dt+a*M0.z) / b;

		BSpline1D::setControl(m-1,2,r1);
		tsPoint[m-1].update();
	}
}

void SplineBranch::enforceConstraints() {
	int i;

	// Interpolate the spline at the center and at three positions around it
	MySMLVec4f W[3][4];
	MySMLVec3f M0,M1[3];
	float ca[3],sa[3],dsdt[3];
	for(i=0;i<3;i++) {
		if(ei[i]==0) {
			crv[i]->kv.basisJet(3,0.0,W[i]);		
			crv[i]->interpolatePoint(0,W[i],0,0,2,M0.data());
			crv[i]->interpolatePoint(0,W[i],1,0,2,M1[i].data());
			dsdt[i] = sqrt(M1[i].x*M1[i].x+M1[i].y*M1[i].y);
			ca[i] = M1[i].x / dsdt[i];
			sa[i] = M1[i].y / dsdt[i];
		}
		else {
            crv[i]->kv.basisJet(crv[i]->m,1.0,W[i]);		
			crv[i]->interpolatePoint(crv[i]->m-3,W[i],0,0,2,M0.data());
			crv[i]->interpolatePoint(crv[i]->m-3,W[i],1,0,2,M1[i].data());
			dsdt[i] = sqrt(M1[i].x*M1[i].x+M1[i].y*M1[i].y);
			ca[i] = -M1[i].x / dsdt[i];
			sa[i] = -M1[i].y / dsdt[i];
		}		
	}

	// Now compute the radius at each adjacent point
	for(i=0;i<3;i++) {
		// Coefficients of r0 and r1 in cosTheta
		float a,b; 
		if(ei[i] == 0) {
			a = W[i][1].x;
			b = W[i][1].y;
		}
		else {
			a = -W[i][1].w;
			b = -W[i][1].z;
		}
		int i1 = (i+1) % 3, i2 = (i+2) % 3;
		float r1 = (dsdt[i]*(ca[i1]*ca[i2]+sa[i1]*sa[i2])-a*M0.z) / b;
		
		if(std::numeric_limits<float>::quiet_NaN() == r1) {
			cout << "NAN" << endl;
		}

		crv[i]->setControl(ri[i],2,r1);
			
		// Refresh timestamp
		crv[i]->tsPoint[ri[i]].update();
	}
}

double SplineBranch::getMaxAngle() {
	/*
	int i;

	float slope[3];
	Vector2D x[3];
	for(i=0;i<3;i++) {
		MySMLVec3f x3 = crv[i]->getControlVec(ri[i]);
		x[i].x = x3.x;
		x[i].y = x3.y;
	}

	for(i=0;i<3;i++) {
		float cs = (x[(i+1)%3] - x[i]).dotProduct(x[(i+2)%3]);
		if(cs < 0) {
			return cs;
		}
	}
	*/
	int i;

	float angle = 0;
	
	MySMLVec3f xc = crv[0]->getControlVec(ei[0]);
	Vector2D x[3];
	for(i=0;i<3;i++) {
		MySMLVec3f x3 = crv[i]->getControlVec(ri[i]) - xc;
		x[i].x = x3.x;
		x[i].y = x3.y;
		x[i].normalize();
	}

	for(i=0;i<3;i++) {
		angle += acos(x[i].dotProduct(x[(i+1)%3]));
	}

	return 2.0 * M_PI - angle;
}

SplineCurve *SplineObject::addCurve(int nPoints,float rho) {
	SplineCurve *nCurve = new SplineCurve(nPoints,rho);
	curves.push_back(nCurve);
	
	// Refresh timestamp
	tsStructure.update();

	return nCurve;
}

SplineCurve *SplineObject::addCurve(SplineCurve *curve) {
	curves.push_back(curve);
	
	// Refresh timestamp
	tsStructure.update();

	return curve;
}

// Is it possible to add a branch at a point?
bool SplineObject::canAddBranch(SplineCurve *curve,int idx) {
	// There must be at least three control points on each side
	return (idx >= 3 && idx <= curve->m - 3);
}
	
SplineCurve *SplineObject::addBranch(SplineCurve *curve,int idx,int nPoints) {
	SplineCurve *sub  = new SplineCurve(nPoints,curve->rho);

	// We use the medial atoms corresponding to the control point to 
	// design the third curve
	MAtom m;
	curve->interpolateAtom(idx,curve->kv.getParmAtKnot(idx+2),m);
	MySMLVec3f X(m.x().x,m.x().y,m.r());
	MySMLVec3f D(m.n().getNormal().x * m.r(),m.n().getNormal().y * m.r(),0);

	// Create the last figure
	for(int i=0;i<nPoints;i++) 
		sub->setControlVec(i,X+D*i);

	return addBranch(curve,idx,sub,0);
}

SplineCurve *SplineObject::addBranch(SplineCurve *curve,int idx,SplineCurve *sub,int subPoint) {
	int i,j;

	if(canAddBranch(curve,idx)) {
		// Create three new curves
		SplineCurve *cut1 = new SplineCurve(idx+1,curve->rho);
		SplineCurve *cut2 = new SplineCurve(curve->m+1-idx,curve->rho);

		// Set the control point values of the first two curves
		for(i=0;i<=idx;i++) 
			cut1->setControlVec(i,curve->getControlVec(i));
		for(i=idx;i<=curve->m;i++) 
			cut2->setControlVec(i-idx,curve->getControlVec(i));

		// Update the branches that already exist at the ends of the figure
		for(j=0;j<(int)branches.size();j++) {
			SplineBranch *br = branches[j];
			for(i=0;i<3;i++) {
				if(br->crv[i]==curve) {
					if(br->ei[i]==0) {
						br->crv[i] = cut1;
						cut1->br[0] = br;
					}
					else if(br->ei[i]==curve->m) {
						br->crv[i] = cut2;
						br->ei[i] = cut2->m;
						br->ri[i] = cut2->m-1;
						cut2->br[1] = br;
					}
				}
			}
		}
		
		// Create a new branch object describing all three things
		SplineBranch *branch = new SplineBranch();
		branch->crv[0] = cut1;
		branch->ei[0] = idx;
		branch->ri[0] = idx-1;
		branch->crv[1] = cut2;
		branch->ei[1] = 0;
		branch->ri[1] = 1;
		branch->crv[2] = sub;
		branch->ei[2] = subPoint;
		branch->ri[2] = (subPoint==0) ? 1 : subPoint-1;		
		branches.push_back(branch);		
		sub->enforceConstraints();
		branch->enforceConstraints();

		// Add branch to curve ends
		cut1->br[1] = branch;
		cut2->br[0] = branch;
		sub->br[subPoint ? 1 : 0] = branch;

		// Delete the old curve
		curves.erase(find(curves.begin(),curves.end(),curve));
		delete curve;

		curves.push_back(cut1);
		curves.push_back(cut2);

		// Add sub if it's not already there
		if(count(curves.begin(),curves.end(),sub) == 0)
			curves.push_back(sub);

		// Play with the branch to make sure that the radius is OK
		MySMLVec3f C = branch->crv[0]->getControlVec(branch->ei[0]);
		Vector2D Cx(C.x,C.y);
		for(int i=0;i<3;i++) {			
			MySMLVec3f F = branch->crv[i]->getControlVec(branch->ri[i]);
			Vector2D Fx(F.x,F.y);
			float tz = F.z + Fx.distanceTo(Cx);
			C.z = C.z > tz ? C.z : tz;
		}
		branch->crv[0]->updateControl(branch->ei[0],C);

		// Refresh timestamp
		tsStructure.update();

		// Return the branch
		return sub;
	}
	return NULL;
}

void SplineObject::removeCurve(int idx) {
	SplineCurve *curve = curves[idx];
	for(int b=0;b<=1;b++) {
		SplineBranch *br = curve->br[b];
		if(br) {
			removeBranch(br);
		}		
	}
	curves.erase(curves.begin() + idx);
	delete curve;

	// Refresh timestamp
	tsStructure.update();
}

void SplineObject::removeBranch(SplineBranch *branch) {
	for(int i=0;i<3;i++) {
		branch->crv[i]->br[branch->ei[i] ? 1 : 0] = NULL;
		branch->crv[i]->enforceConstraints();
	}
	branches.erase(find(branches.begin(),branches.end(),branch));
	delete branch;

	// Refresh timestamp
	tsStructure.update();
}

// Save to a registry
void SplineObject::save(Registry &folder) {
	unsigned int i;

	folder.setIntValue("curveCount",(int)curves.size());
	folder.setIntValue("branchCount",(int)branches.size());
	for(i=0;i<curves.size();i++) {
		SplineCurve *curve = curves[i];
		Registry &cr = folder.getSubFolder("curve[%d]",i);
		cr.setIntValue("controlCount",curve->size());
		cr.setDoubleValue("rho",curve->rho);
		for(int j=0;j<curve->size();j++) {
			MySMLVec3f cp = curve->getControlVec(j);
			cr.setDoubleValue("x[%d]",cp.x,j);
			cr.setDoubleValue("y[%d]",cp.y,j);
			cr.setDoubleValue("r[%d]",cp.z,j);
		}
	}
	for(i=0;i<branches.size();i++) {
		SplineBranch *br = branches[i];
		Registry &cr = folder.getSubFolder("branch[%d]",i);
		for(int j=0;j<3;j++) {
			int idx = find(curves.begin(),curves.end(),br->crv[j]) - curves.begin();
			cr.setIntValue("curve[%d]",idx,j);
			cr.setIntValue("ei[%d]",br->ei[j],j);
			cr.setIntValue("ri[%d]",br->ri[j],j);
		}
	}
}

// Load from registry
void SplineObject::load(Registry &folder) {
	int i;
		
	// Remove all curves and branches
	while(curves.size())
		removeCurve(0);
	branches.clear();
	
	// Load data
	int nc = folder.getIntValue("curveCount",0);
	int bc = folder.getIntValue("branchCount",0);

	for(i=0;i<nc;i++) {
		Registry &cr = folder.getSubFolder("curve[%d]",i);
		int cc = cr.getIntValue("controlCount",0);
		if(cc > 0) {
			double rho = cr.getDoubleValue("rho",0.25);
			SplineCurve *curve = new SplineCurve(cc,rho);
			for(int j=0;j<cc;j++) {
				MySMLVec3f cp;
				cp.x = cr.getDoubleValue("x[%d]",0,j);
				cp.y = cr.getDoubleValue("y[%d]",0,j);
				cp.z = cr.getDoubleValue("r[%d]",0,j);
				curve->setControlVec(j,cp);
			}
			curves.push_back(curve);
		}		
	}

	for(i=0;i<bc;i++) {
		Registry &cr = folder.getSubFolder("branch[%d]",i);
		SplineBranch *br = new SplineBranch();
		for(int j=0;j<3;j++) {
			br->crv[j] = curves[cr.getIntValue("curve[%d]",0,j)];
			br->ei[j] = cr.getIntValue("ei[%d]",0,j);
			br->ri[j] = cr.getIntValue("ri[%d]",0,j);
			br->crv[j]->br[br->ei[j] ? 1 : 0] = br;
		}
		branches.push_back(br);
	}

	// Refresh timestamp
	tsStructure.update();
}


RegularSplineSample::RegularSplineSample(SplineObject *spline) :
samples(spline->curves.size())
{
	this->spline = spline;
}


RegularSplineSample::RegularSplineSample(SplineObject *spline,int nPerSegment) :
samples(spline->curves.size())
{
	this->spline = spline;

	// Initialize the arrays
	init(nPerSegment);
}

long RegularSplineSample::getTimestampForSegment(int iCurve,int iSegment) {
	long tsMax = 0;
	for(int i=iSegment;i<iSegment+4;i++) {
		long ts = spline->curves[iCurve]->tsPoint[i];
		tsMax = ts > tsMax ? ts : tsMax;
	}
	return tsMax;
}

void RegularSplineSample::update() {
	// Repeat for each medial curve in the object
	for(unsigned int iCurve=0;iCurve<spline->curves.size();iCurve++) {

		// A shorthand
		SplineCurve *curve = spline->curves[iCurve];

		// Compute the atoms, never including the last one
		for(int iSeg = 0;iSeg < curve->size()-3;iSeg++) {
			
			// Get the current segment
			SplineSegment &segment = samples[iCurve][iSeg];

			// Get the new timestamp
			long tscp = getTimestampForSegment(iCurve,iSeg);
			
			// Update if needed 
			if(tscp > (long)segment.ts) {
				updateSegment(iCurve,iSeg);
				(long)segment.ts = tscp;
			}
		}
	}
}

void RegularSplineSample::init(int nPerSegment) {

	// Initialize the timestamp
	tsOverall = spline->tsStructure;

	// Repeat for each medial curve in the object
	for(unsigned int iCurve=0;iCurve<spline->curves.size();iCurve++) {

		// A shorthand
		SplineCurve *curve = spline->curves[iCurve];

		// Resize the lsit
		samples[iCurve].resize(curve->size()-3);
		
		// Compute the atoms, never including the last one
		for(int iSeg = 0;iSeg < curve->size()-3;iSeg++) {
			
			// Get the current segment
			SplineSegment &segment = samples[iCurve][iSeg];
		
			// Add the points to the current segment
			double t = curve->kv.getParmAtKnot(iSeg+3);
			double t1 = curve->kv.getParmAtKnot(iSeg+4);
			double tStep = (t1 - t) / nPerSegment;			
			for(int i=0;i<=nPerSegment;i++) {
				segment.points.push_back(SplineSamplePoint(t,iSeg));
				t+=tStep;
			}			

			// Clear the values at the last point (no need?)
			// segment.points.back().cv[0] = 0;
			// segment.points.back().cv[1] = 0;

			// Set the timestamp
			updateSegment(iCurve,iSeg);		
			segment.ts = getTimestampForSegment(iCurve,iSeg);
		}
	}
}
/*
void RegularSplineSample::updateSamples(int iCurve,int iStart,int nPoints) {
	SplineCurve *curve = spline->curves[iCurve];

	// Get the range of the affected points
	int segStart = max(0,iStart-3);
	int segEnd = min(curve->size()-3,iStart+nPoints);
	
	// Parse the list and update the affected points
	list<SplineSamplePoint> &lst = samples[iCurve];
	list<SplineSamplePoint>::iterator it;
	for(it=lst.begin();it!=lst.end();it++) {
		if(segStart <= it->seg && it->seg < segEnd) {
			curve->interpolateAtom(it->seg,it->t,it->atom);
		}
	}
}
*/
void RegularSplineSample::updateSegment(int iCurve,int iSeg) {
	SplineCurve *curve = spline->curves[iCurve];
	SplineSegment &segment = samples[iCurve][iSeg];

	// Parse the list and update the affected points
	list<SplineSamplePoint> &lst = segment.points;
	list<SplineSamplePoint>::iterator it;
	for(it=lst.begin();it!=lst.end();it++) {
		curve->interpolateAtomJet(iSeg,it->t,it->atom,it->M[0],it->M[1],it->M[2]);
	}
}


ArcSplineSample::ArcSplineSample(SplineObject *object,int nPerSegmentMin) 
: RegularSplineSample(object)
{
	this->nPerSegmentMin = nPerSegmentMin;
	init(nPerSegmentMin * 1.2);
}

void ArcSplineSample::updateSegment(int iCurve,int iSeg) {

	SplineCurve *curve = spline->curves[iCurve];
	SplineSegment &segment = samples[iCurve][iSeg];

	// Parse the list and update the affected points
	list<SplineSamplePoint> &lst = segment.points;
	list<SplineSamplePoint>::iterator a,b,c;

	// Find the starting and ending points for the operation
	for(a=lst.begin();a!=lst.end();a++) {
		curve->interpolateAtomJet(iSeg,a->t,a->atom,a->M[0],a->M[1],a->M[2]);
	}

	// Compute the total length and the minimum length
	double md[3] = {0.0,0.0,0.0};
	int nSamples = nPerSegmentMin;
	for(b=lst.begin(),a=b++;b!=lst.end();a++,b++) {
		// Compute the distances between atoms
		computeDistances(*a,*b);
		md[0] += a->dn[0]; md[1] += a->dn[1]; md[2] += a->dn[2];
	}
	md[0] /= nSamples; md[1] /= nSamples; md[2] /= nSamples;

	// The tolerance: the lowest t difference that we will allow during splitting
	const double tolerance = 1e-6;

	// Count the inserts and deletes
	int nInsert = 0,nDelete = 0;

	// This is a loop with two sequential pointers a, b
	for(b=lst.begin(),a=b++;b!=lst.end();) {
		double dt = fabs(b->t - a->t);

		// Check if the two samples are farther than apart than allowed but still within tolerance
		if((a->dn[0] > md[0] || a->dn[1] > md[1] || a->dn[2] > md[2]) && dt > tolerance) {
			// Perform the split
			c = lst.insert(b,SplineSamplePoint(0.5 * (a->t + b->t),a->seg));
			curve->interpolateAtomJet(c->seg,c->t,c->atom,c->M[0],c->M[1],c->M[2]);

			// Compute the new distances
			computeDistances(*a,*c);
			computeDistances(*c,*b);

			// Reassign b
			b=c;

			nInsert++;
		}
		else {
			a++;b++;
		}
	}

	md[0] /= 2;md[1] /= 2;md[2] /= 2;

	// Now trim the samples (remove the samples that are unnecessarily close to each other)
	for(c=lst.begin(),a=c++,b=c++;c!=lst.end();) {
		double d0 = a->dn[0] + b->dn[0];
		double d1 = a->dn[1] + b->dn[1];
		double d2 = a->dn[2] + b->dn[2];
		if(d0 < md[0] && d1 < md[1] && d2 < md[2]) {
			lst.erase(b);
			b = c;
			c++;
			computeDistances(*a,*b);

			nDelete++;
		}
		else {
			a++;b++;c++;
		}
	}

	// cout << "(" << nInsert << "/" << nDelete << "/" << lst.size() << ") ";
}

void RegularSplineSample::computeDistances(SplineSamplePoint &a,SplineSamplePoint &b) {
	Vector2D dm = b.atom.x(MAtom::MIDDLE) - a.atom.x(MAtom::MIDDLE);
	Vector2D dl = b.atom.x(MAtom::LEFT) - a.atom.x(MAtom::LEFT);
	Vector2D dr = b.atom.x(MAtom::RIGHT) - a.atom.x(MAtom::RIGHT);

	a.dn[0] = dm.twoNorm();
	a.dn[1] = dl.twoNorm();
	a.dn[2] = dr.twoNorm();
	a.cv[0] = dm.dotProduct(dl);
	a.cv[1] = dm.dotProduct(dr);
}

/*	
	// samples.reserve(spline->curves.size());
	for(unsigned int iCurve=0;iCurve<spline->curves.size();iCurve++) {
		// Create a new list
		
		
		// samples.push_back(vector<MAtom>());
		// vector<MAtom> &vec = samples.back();
		list<SplineSamplePoint> &lst = samples[iCurve];

		// Keep track of the total length of the spline
		double lTotal = 0;
		MAtom *aLast = NULL;

		// Create a curve for this spline
		for(int seg = 0;seg < curve->size()-3;seg++) {
			double t = curve->kv.getParmAtKnot(seg+3);
			lst.push_back(SplineSamplePoint(t,seg));
			curve->interpolateAtom(seg,t,lst.back().atom);

			if(aLast == NULL) {
				aLast = &lst.back().atom;
			}
			else {
				lTotal += aLast->x().distanceTo(lst.back().atom.x());
				aLast = &lst.back().atom;
			}
		}

		// Make sure we are using the right length
		double md = max(lTotal * maxDist,maxDist);

		// cout << endl << "lTotal = " << lTotal << endl;
		// cout << "md = " << lTotal << md;

		// Compute the last atom
		lst.push_back(SplineSamplePoint(1.0,curve->size()-4));
		curve->interpolateAtom(curve->size()-4,1.0,lst.back().atom);

		bool finished;
		do {
			finished = true;
			
			list<SplineSamplePoint>::iterator it0=lst.begin(), it1=++(lst.begin());
			while(it1 != lst.end()) {
				SplineSamplePoint &n0 = *it0;
				SplineSamplePoint &n1 = *it1;	
				MAtom &a0 = n0.atom;
				MAtom &a1 = n1.atom;
				if(a0.x().distanceTo(a1.x()) > md ||
					a0.x(MAtom::LEFT).distanceTo(a1.x(MAtom::LEFT)) > md ||
					a0.x(MAtom::RIGHT).distanceTo(a1.x(MAtom::RIGHT)) > md) 
				{
					it0 = lst.insert(it1,SplineSamplePoint(0.5 * (n0.t + n1.t),n0.seg));
					curve->interpolateAtom(it0->seg,it0->t,it0->atom);
					it1 = it0;
					it1++;
					finished = false;

					if(lst.size() > 1000) {
						cout << endl;
						cout << "OVERRUN " << endl << "md = " << lTotal << md << endl					;
						cout << "D0 " << a0.x().distanceTo(a1.x()) << endl;
						cout << "D1 " << a0.x(MAtom::LEFT).distanceTo(a1.x(MAtom::LEFT)) << endl;
						cout << "D2 " << a0.x(MAtom::RIGHT).distanceTo(a1.x(MAtom::RIGHT)) << endl;
						//finished = true;
						//break;
					}
				}
				else {
					it0++;
					it1++;
				}
			}
		} while(!finished);

		// cout << "arclst : " << lst.size() << endl;
	}
	*/
	/*
	// Compute arclength of each curve
	samples.reserve(spline->curves.size());
	for(unsigned int iCurve=0;iCurve<spline->curves.size();iCurve++) {

		// A shorthand
		SplineCurve *curve = spline->curves[iCurve];
		
		// Create a new list
		samples.push_back(vector<MAtom>());
		vector<MAtom> &vec = samples.back();
		// vec.reserve(curve->size()*nPerSegment+1);

		// Compute the arclength of the curve

		// Compute the atoms, never including the last one
		for(int seg = 0;seg < curve->size()-3;seg++) {
			double t = curve->kv.getParmAtKnot(seg+3);
			double tStep = (curve->kv.getParmAtKnot(seg+4)-t) / (nPerSegment);			
			for(int i=0;i<nPerSegment;i++) {
				vec.push_back(MAtom());
				curve->interpolateAtom(seg,t,vec.back());
				t+=tStep;
			}
		}

		// Compute the last atom
		vec.push_back(MAtom());
		curve->interpolateAtom(curve->size()-4,1.0,vec.back());
	}*/
// }

