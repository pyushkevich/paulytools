#include "mreps2D.h"
#include <support.h>
#include <BrentLinearMethod.h>


const int BAtom::NOT_IN_SHAPE = -1;

const int MAtom::MIDDLE = 0;
const int MAtom::LEFT   = 1;
const int MAtom::RIGHT  = 2;
const int MAtom::HEAD   = 3;
const int MAtom::TAIL   = 4;

const int MNode::SUCC      = 0;
const int MNode::PRED      = 1;
const int MNode::BRANCH   = 1;
const int MNode::CTF_START = 6;
// Careful - BRANCH+LEFT must make sense

// Some useful code for taking logs
const double M_LOG2 = log(2.0);
const double M_LOG2INV = 1.0/log(2.0);

inline double exp2(double x) {
    return exp(x*M_LOG2);
}

inline double log2(double x) {
    return log(x) * M_LOG2INV;
}

typedef list<BRef>::iterator BRefItr;

template <class T> void destroyList(list<T> &lst) {
    for(typename list<T>::iterator p=lst.begin();p!=lst.end();p++) {
	T t = *p;
	delete t;
    }
    lst.clear();
}

bool MAtom::inBoundingBox(const Vector2D &x1,const Vector2D &x2) const 
{
    return (x().x >= x1.x && x().x <= x2.x && x().y >= x1.y && x().y <= x2.y);
}

bool MAtom::inBoundingBox(double x0, double y0, double x1, double y1) const 
{
    return (x().x >= x0 && x().x <= x1 && x().y >= y0 && x().y <= y1);
}

/**
 * This method computes all the secondary attributes of the atom
 */
void MAtom::computeSecondary() {
    double fa = m_ba[MIDDLE].slope;
    m_ba[LEFT].slope = fa + m_oa;
    m_ba[RIGHT].slope = fa - m_oa;
    m_ba[HEAD].slope = fa;
    m_ba[TAIL].slope = M_PI + fa;


    double ca = cos(fa),sa = sin(fa),co = cos(m_oa),so = sin(m_oa);
    m_ba[MIDDLE].n.x = ca;
    m_ba[MIDDLE].n.y = -sa;
    m_ba[RIGHT].n.x = ca*co+sa*so;
    m_ba[RIGHT].n.y = -(sa*co-ca*so);
    m_ba[LEFT].n.x = ca*co-sa*so;
    m_ba[LEFT].n.y = -(sa*co+ca*so);
    m_ba[HEAD].n.x = ca;
    m_ba[HEAD].n.y = -sa;
    m_ba[TAIL].n.x = -ca;
    m_ba[TAIL].n.y = sa;

    double rEnd = m_r * m_endRC;
    m_ba[LEFT].x = m_ba[MIDDLE].x + m_ba[LEFT].n*m_r;
    m_ba[RIGHT].x = m_ba[MIDDLE].x + m_ba[RIGHT].n*m_r;
    m_ba[HEAD].x = m_ba[MIDDLE].x + m_ba[HEAD].n*rEnd;
    m_ba[TAIL].x = m_ba[MIDDLE].x + m_ba[TAIL].n*rEnd;

    m_ba[MIDDLE].scale = m_ba[LEFT].scale = m_ba[RIGHT].scale = m_r * m_rho;
    m_ba[HEAD].scale = m_ba[TAIL].scale = rEnd * m_rho;

    m_tsSecondary = m_tsPrimary;
}  

// Read atom from a folder
void MAtom::read(Registry &folder) {
    // Read in the standard values
    x(Vector2D(folder.getDoubleValue("x",0),folder.getDoubleValue("y",0)));
    r(folder.getDoubleValue("r", 0));
    m_ba[MIDDLE].slope = folder.getDoubleValue("ta", 0);
    m_oa = folder.getDoubleValue("to", M_PI/2);

    // Polarity is a little weird.  If there is a generic polarity flag, we read it.
    // and use it as a default when reading specific location polarities
    double p = folder.getDoubleValue("polarity",1.0);
    polarity(p,MIDDLE);
    polarity(folder.getDoubleValue("leftPolarity",p),LEFT);
    polarity(folder.getDoubleValue("rightPolarity",p),RIGHT);
    polarity(folder.getDoubleValue("endPolarity",p),HEAD);
    polarity(folder.getDoubleValue("endPolarity",p),TAIL);

    // Finally, read the end cap value
    endRC(folder.getDoubleValue("endRC", 1.0));
}

void MAtom::write(Registry &folder) {
    // Write the standard values
    folder.setDoubleValue("x",x().x);
    folder.setDoubleValue("y",x().y);
    folder.setDoubleValue("r",r());
    folder.setDoubleValue("ta",m_ba[MIDDLE].slope);
    folder.setDoubleValue("to",m_oa);

    // Now write the polarities.  If all polarities are equal then only write one polarity
    if(polarity(MIDDLE) == polarity(LEFT) && polarity(MIDDLE) == polarity(RIGHT) && 
	    polarity(MIDDLE) == polarity(HEAD) && polarity(MIDDLE) == polarity(TAIL)) {
	folder.setDoubleValue("polarity",polarity(LEFT));
	// Should delete the left/right polarities from the registry
    }
    else {
	folder.setDoubleValue("leftPolarity",polarity(LEFT));
	folder.setDoubleValue("rightPolarity",polarity(RIGHT));
	folder.setDoubleValue("endPolarity",polarity(HEAD));
    }

    // Finally, read the end cap value
    folder.setDoubleValue("endRC", endRC());
}

// Get data into a vector that is useable for optimization
Vector MAtom::getData(bool endRC) {
    int dim = endRC ? 6 : 5;
    Vector v(dim);

    v(0) = m_ba[MIDDLE].x.x;
    v(1) = m_ba[MIDDLE].x.y;
    v(2) = log2(m_r);
    v(3) = m_ba[MIDDLE].slope;
    v(4) = m_oa;
    if(endRC)
	v(5) = log2(m_endRC);

    return v;
}

// Read data from a vector useable for optimization
int MAtom::setData(const Vector &v,bool endRC,int idx) {
    dassert((endRC && v.size()>=6) || v.size()>=5);
    m_ba[MIDDLE].x.x = v(idx++);
    m_ba[MIDDLE].x.y = v(idx++);
    m_r = exp2(v(idx++));
    m_ba[MIDDLE].slope = v(idx++);
    m_oa = v(idx++);
    if(endRC)
	m_endRC = exp2(v(idx++));
    m_tsPrimary++;
    return idx;
}

/**********************************************************************************
 * Medial Node Class Implementation                                               *
 *********************************************************************************/

// Get a group of all of the neighbors of this atom
MNodeGroup MNode::getNeighbors() {
    MNodeGroup group(m_graph);
    for(int i=0;i<m_links.size();i++) {
	MNode *node = (m_links[i]->head()==this) ? m_links[i]->tail() : m_links[i]->head();
	group.append(node);
    }
    return group;
}

void MNode::clearTShape() {
    for(int i=LEFT;i<=TAIL;i++) {
	m_ba[i].tShape = BAtom::NOT_IN_SHAPE;
	m_ba[i].bnd = NULL;
    }
}

/**********************************************************************************
 * Medial Link Class Implementation                                               *
 *********************************************************************************/
const int MLink::CURVE = 0;
const int MLink::BRANCH = 1;
const int MLink::CTF = 2;
const int MLink::WRONG_TYPE = -1;


// Constructor
BranchLink::BranchLink(MGraph *graph,MNode *atom1,MNode *atom2,
	int bIndex1,int bIndex2,bool flip) :
MLink(graph) 
{
    dassert(bIndex1 > MAtom::MIDDLE && bIndex1 <= MAtom::TAIL);
    dassert(bIndex2 > MAtom::MIDDLE && bIndex2 <= MAtom::TAIL);

    m_bRef[0].node = atom1;
    m_bRef[1].node = atom2;
    m_bRef[0].bIndex = bIndex1;
    m_bRef[1].bIndex = bIndex2;
    m_flip = flip;
}

// Constructor
BranchLink::BranchLink(MGraph *graph,const BRef &tail,
	const BRef &head,bool flip) :
MLink(graph) 
{
    m_bRef[0] = tail;
    m_bRef[1] = head;
    m_flip = flip;
}

MLink *CurveLink::makeCopyWithSub(MGraph *newGraph,MNode *newTail,MNode *newHead) {
    return new CurveLink(newGraph,newTail,newHead);
}

MLink *BranchLink::makeCopyWithSub(MGraph *newGraph,MNode *newTail,MNode *newHead) {
    return new BranchLink(newGraph,newTail,newHead,m_bRef[0].bIndex,m_bRef[1].bIndex,m_flip);
}

/**********************************************************************************
 * Medial Node Group Class Implementation                                         *
 *********************************************************************************/
// Add an atom to the group

MNodeGroup MNodeGroupIF::ggInBox(const Vector2D &min,const Vector2D &max) {
    MNodeGroup ng(graph());
    for(MNode *n = first();n!=NULL;n = next()) {
	if(n->inBoundingBox(min, max)) {
	    ng.append(n);
	}
    }
    return ng;
}

MNodeGroup MNodeGroupIF::ggSelected() {
    MNodeGroup ng(graph());
    for(MNode *n = first();n!=NULL;n = next()) {
	if(n->selected()) {
	    ng.append(n);
	}
    }
    return ng;
}



MNode *MNodeGroupIF::closestAtom(const Vector2D &x) {
    MNode *node = first();
    if(node==NULL)
	return NULL;

    MNode *bestNode = node;
    double bestDist = node->x().distanceTo(x);

    while(node = next()) {
	int dist = (int) node->x().distanceTo(x);
	if (dist < bestDist) {
	    bestDist = dist;
	    bestNode = node;
	}		
    }

    return node;
}


void MNodeGroupIF::getExtents(Vector2D &extMin,Vector2D &extMax) 
{
    dassert(first()!=NULL);

    extMin = extMax = first()->x();
    for (MNode *n = first();n!=NULL;n = next()) {
	for(int bi=MAtom::MIDDLE;bi<=MAtom::TAIL;bi++) {
	    if(bi==MAtom::HEAD && n->getSuccessor())
		continue;
	    if(bi==MAtom::TAIL && n->getPredecessor())
		continue;

	    extMin = Vector2D::minimal(extMin,n->x(bi));
	    extMax = Vector2D::maximal(extMax,n->x(bi));
	}
    }
}

void MNodeGroupIF::scaleToFit(double x,double y,double w,double h) {
    dassert(first()!=NULL);
    dassert(w > 0 && h > 0);
    Vector2D ext1,ext2;

    // Find out how much to scale by
    getExtents(ext1,ext2);
    double maxdim = (ext2-ext1).infinityNorm();

    // Scale, keeping the top extent in place
    scaleBy(w/maxdim,ext1);

    // Translate top extent to x,y
    translateBy(Vector2D(x,y)-ext1);
}

Vector2D MNodeGroupIF::getCOG() 
{
    dassert(first()!=NULL);

    Vector2D cog(0,0,1);
    int count = 0;

    for (MNode *n = first();n!=NULL;n = next()) {
	cog += n->x();
	count++;
    }

    cog *= 1.0 / count;
    return cog;
}

void MNodeGroupIF::scaleBy(double scalefact, const Vector2D &cog)
{
    for (MNode *n = first();n!=NULL;n = next()) {
	n->x((n->x() - cog) * scalefact + cog);
	n->r(n->r()*scalefact);
    }
}

void MNodeGroupIF::elongate(double scalefact, const Vector2D &cog)
{
    for (MNode *n = first();n!=NULL;n = next()) {
	n->x((n->x() - cog) * scalefact + cog);
    }
}

void MNodeGroupIF::thin(double scalefact, const Vector2D &cog)
{
    for (MNode *n = first();n!=NULL;n = next()) {
	n->r(n->r()*scalefact);
    }
}

void MNodeGroupIF::rotateBy(double angledeg, const Vector2D &cog) {
    Vector2D x,xc;
    double angle = M_PI*angledeg/180.0;

    for (MNode *n = first();n!=NULL;n = next()) {
	xc = n->x() - cog;
	x.x = xc.x*cos(angle)-xc.y*sin(angle); 
	x.y = xc.x*sin(angle)+xc.y*cos(angle); 
	n->x(x+cog);
	n->fa(n->fa()-angle);
    }
}

void MNodeGroupIF::changeObjectAngleBy(double angledeg) {
    for (MNode *n = first();n!=NULL;n = next()) {
	n->oaDeg(n->oaDeg()-angledeg);
    }
};						

void MNodeGroupIF::translateBy(const Vector2D &dx) {
    for (MNode *n = first();n!=NULL;n = next()) {
	n->x(n->x()+dx);
    }   
}

vector<MNode *> MNodeGroupIF::getArray() {
    vector<MNode *> out;
    for(MNode *n = first(); n!=NULL; n = next()) {
	out.push_back(n);
    }
    return out;
}

vector<MAtom> MNodeGroupIF::getAtomArray() {
    vector<MAtom> out;
    for(MNode *n = first(); n!=NULL; n = next()) {
	out.push_back(*n);
    }
    return out;
}

void MNodeGroupIF::rho(double rho) {
    for(MNode *n = first(); n!=NULL; n = next()) {
	n->rho(rho);
    }
}

Matrix MNodeGroupIF::getBoundaryPoints() {
    vector<Vector2D> xList;

    for(MNode *n = first(); n!=NULL; n = next()) {
	for(int a=MNode::LEFT;a<=n->MNode::TAIL;a++) {
	    if(n->bAtom(a).bnd != NULL) {
		xList.push_back(n->bAtom(a).x);
	    }
	}
    }

    Matrix m(xList.size(),2);
    for(int i=0;i<xList.size();i++) {
	m(i,0) = xList[i].x;
	m(i,1) = xList[i].y;
    }

    return m;
}

void MNodeGroupIF::fitEllipse(Vector2D &center,Vector2D &major,Vector2D &minor) {
    Matrix m = getBoundaryPoints();
    Vector mean,eval;
    Matrix evec;

    m.computePCA(mean,eval,evec);
    double sMajor = sqrt(eval(1));
    double sMinor = sqrt(eval(0));

    center.x = mean(0);
    center.y = mean(1);

    major.x = sMajor * evec(0,1);
    major.y = sMajor * evec(1,1);

    minor.x = sMinor * evec(0,0);
    minor.y = sMinor * evec(1,0);
}


/**********************************************************************************
 * Medial Figure Class Implementation                                             *
 *********************************************************************************/


/**********************************************************************************
 * BRef Class Implementation                                                      *
 *********************************************************************************/
BRef BRef::next(bool cw) {
    BRef br;

    if(cw) {
	if(bIndex==MAtom::LEFT) {
	    if(node->hasSuccessor()) {
		br.bIndex = MAtom::LEFT;
		br.node = node->getSuccessor();
	    }
	    else {
		br.bIndex = MAtom::HEAD;
		br.node = node;
	    }
	}
	else if(bIndex==MAtom::RIGHT) {
	    if(node->hasPredecessor()) {
		br.bIndex = MAtom::RIGHT;
		br.node = node->getPredecessor();
	    }
	    else {
		br.bIndex = MAtom::TAIL;
		br.node = node;
	    }
	}
	else if(bIndex==MAtom::HEAD) {
	    br.bIndex = MAtom::RIGHT;
	    br.node = node;
	}
	else if(bIndex==MAtom::TAIL) {
	    br.bIndex = MAtom::LEFT;
	    br.node = node;
	}
	else {
	    dassert(0);
	}
    }
    else {
	if(bIndex==MAtom::RIGHT) {
	    if(node->hasSuccessor()) {
		br.bIndex = MAtom::RIGHT;
		br.node = node->getSuccessor();
	    }
	    else {
		br.bIndex = MAtom::HEAD;
		br.node = node;
	    }
	}
	else if(bIndex==MAtom::LEFT) {
	    if(node->hasPredecessor()) {
		br.bIndex = MAtom::LEFT;
		br.node = node->getPredecessor();
	    }
	    else {
		br.bIndex = MAtom::TAIL;
		br.node = node;
	    }
	}
	else if(bIndex==MAtom::HEAD) {
	    br.bIndex = MAtom::LEFT;
	    br.node = node;
	}
	else if(bIndex==MAtom::TAIL) {
	    br.bIndex = MAtom::RIGHT;
	    br.node = node;
	}
	else {
	    dassert(0);
	}
    }

    return br;
}

// Check if it's same as another bref
inline bool BRef::equals(const BRef &bref) const {
    return ((node==bref.node) && (bIndex == bref.bIndex));
}

Timer tBoundlet;

/**********************************************************************************
 * Boundlet (customized for MReps) Class Implementation                           *
 *********************************************************************************/
Vector *MRepBoundlet::tnCache = NULL;

inline void MRepBoundlet::compute(const BAtom &b0,const BAtom &b1) 
{
    // Compute distance between ends
    double d = b0.x.distanceTo(b1.x);

    // Compute scaled normals
    double tx0 = -d * b0.n.y;
    double ty0 = d * b0.n.x;
    double tx1 = -d * b1.n.y;
    double ty1 = d * b1.n.x;

    tx0 =  (b0.bndFlipped) ? -tx0 : tx0;
    ty0 =  (b0.bndFlipped) ? -ty0 : ty0;
    tx1 =  (b1.bndFlipped) ? -tx1 : tx1;
    ty1 =  (b1.bndFlipped) ? -ty1 : ty1;

    // Compute all the spline parameters
    ax(3) = 2.0 * (b0.x.x - b1.x.x) + tx0 + tx1;
    ax(2) = 3.0 * (b1.x.x - b0.x.x) - 2*tx0 - tx1;
    ax(1) = tx0;
    ax(0) = b0.x.x;

    ay(3) = 2.0 * (b0.x.y - b1.x.y) + ty0 + ty1;
    ay(2) = 3.0 * (b1.x.y - b0.x.y) - 2*ty0 - ty1;
    ay(1) = ty0;
    ay(0) = b0.x.y;

    bx(2) = 3.0 * ax(3);
    bx(1) = 2.0 * ax(2);
    bx(0) = ax(1);

    by(2) = 3.0 * ay(3);
    by(1) = 2.0 * ay(2);
    by(0) = ay(1);

    // Compute the scale parameters
    bScale = b0.scale;
    kScale = b1.scale - b0.scale;

    // Compute the polarity
    polarity = (fabs(b0.polarity) < fabs(b1.polarity)) ? b0.polarity : b1.polarity;

    // Clear the atom cache
    for(int i=0;i<=boundletCacheSize;i++) 
	cache[i].scale = 0.0;

    // Initialize the tCache if necessary
    if(tnCache == NULL) {
	tnCache = new Vector[boundletCacheSize + 1];
	for(int i=0;i <= boundletCacheSize;i++) {
	    double t = ((double)i) / boundletCacheSize;
	    tnCache[i].setSize(4);
	    tnCache[i](0) = 1;
	    tnCache[i](1) = t;
	    tnCache[i](2) = t * t;
	    tnCache[i](3) = t * tnCache[i](2);
	}
    }
};


inline MRepBoundlet::MRepBoundlet() :ax(4),ay(4),bx(4),by(4) {
    m_tsStart = m_tsEnd = -1;
}

/****************************************************************************
 * Automatic medial axis interpolation
 ***************************************************************************/
MAtom MAxelet::computeAtom(const BAtom &L,const BAtom &R) {
    double x1 = L.x.x;
    double y1 = L.x.y;
    double x2 = R.x.x;
    double y2 = R.x.y;

    double u1 = L.n.x;
    double v1 = L.n.y;
    double u2 = R.n.x;
    double v2 = R.n.y;

    double r1,r2;
    double den = v1 * u2 - u1 * v2;
    if(fabs(den) < rootThresh) {
	r1 = r2 = L.x.distanceTo(R.x) / 2;
    }
    else {
	r1 = (u2 * (y1-y2) - v2 * (x1-x2)) / den;
	r2 = (u1 * (y1-y2) - v1 * (x1-x2)) / den;
    }

    MAtom ma;
    ma.x((L.x - L.n * r1 + R.x - R.n * r2) * 0.5);
    ma.r((r1 + r2 ) * 0.5);

    double slope1 = (L.x - ma.x()).getSlope();
    double slope2 = (R.x - ma.x()).getSlope();
    ma.fa(-(slope1 + slope2) / 2.0);
    ma.oa(-(slope1 - slope2) / 2.0);

    if(ma.x(MAtom::LEFT).distanceTo(L.x) > ma.x(MAtom::LEFT).distanceTo(R.x)) {
	ma.fa(M_PI + ma.fa());
	ma.oa(-ma.oa());
    }

    return ma;
}

Vector2D MAxelet::computeGradient(double s,double t) {
    // S and T vectors
    Vector S(4,1.0,s,s*s,s*s*s),T(4,1.0,t,t*t,t*t*t);
    Vector2D x1,x2,x1_p,x2_p;

    // Compute up to the first derivative
    x1.x = bLeft->ax.dotProduct(S);
    x1.y = bLeft->ay.dotProduct(S);
    x1_p.x = bLeft->bx.dotProduct(S);
    x1_p.y = bLeft->by.dotProduct(S);
    x2.x = bRight->ax.dotProduct(T);
    x2.y = bRight->ay.dotProduct(T);
    x2_p.x = bRight->bx.dotProduct(T);
    x2_p.y = bRight->by.dotProduct(T);

    // Lengths of tangents
    double l1 = x1_p.twoNorm();
    double l2 = x2_p.twoNorm();

    // Second derivatives
    Vector2D x1_pp(bLeft->bx(2)*2*s + bLeft->bx(1),bLeft->by(2)*2*s + bLeft->by(1));
    Vector2D x2_pp(bRight->bx(2)*2*t + bRight->bx(1),bRight->by(2)*2*t + bRight->by(1));

    // Curvatures
    double k1 = (x1_p.x * x1_pp.y - x1_pp.x * x1_p.y) / (l1*l1*l1);
    double k2 = (x2_p.x * x2_pp.y - x2_pp.x * x2_p.y) / (l2*l2*l2);

    // Normalized variables
    Vector2D u1 = x1_p / l1, u2 = x2_p / l2;
    Vector2D u1_p = x1_p.getNormal() * k1;
    Vector2D u2_p = x2_p.getNormal() * k2;

    // Directional derivatives
    double Ft = - x2_p.dotProduct(u2-u1) - u2_p.dotProduct(x2-x1);
    double Fs = x1_p.dotProduct(u2-u1) + u1_p.dotProduct(x2-x1);	

    return Vector2D(Fs,Ft);
}

// Perform actual computation at a position t
void MAxelet::compute(CurveLink *link) {
    const BAtom &b1 = link->tail()->bAtom(MAtom::LEFT);
    const BAtom &b2 = link->head()->bAtom(MAtom::RIGHT);

    int t1 = (int)b1.tShape;
    int t2 = (int)b2.tShape;
    bLeft = (t1 >= 0) ? &b1.bnd->getBoundlet(t1) : NULL;
    bRight = (t2 >= 0) ? &b2.bnd->getBoundlet(t2) : NULL;

    // Nothing to do if bLeft or bRight are NULL
    if(bLeft==NULL || bRight==NULL) {
	traceable = false;
    }
    else {
	// Compute the gradient in s,t at 0 and 1.  This is important
	// because the gradient must be in the right direction for	
	// a valid medial axis to exist
	Vector2D g1 = computeGradient(0,1);
	Vector2D g2 = computeGradient(1,0);

	// Determine validity
	traceable = !(g1.x > 0 || g1.y > 0 || g2.x > 0 || g2.y > 0);
    }

    // Update time stamps
    m_tsStart = link->tail()->timestamp();
    m_tsEnd = link->head()->timestamp();

    // Create or clear the cache
    if(cache == NULL) {
	cache = new MAtom*[boundletCacheSize + 1];
	memset(cache,0,(boundletCacheSize+1)*sizeof(MAtom *));
    } else {
	clearCache();
    }
}

void MAxelet::clearCache() {
    for(int i= 0;i<=boundletCacheSize;i++) {
	if(cache[i] != NULL) {
	    delete cache[i];
	    cache[i] = NULL;
	}
    }
}

MAxelet::MAxelet() {
    // Clear the cache
    cache = NULL;
    m_tsStart = m_tsEnd = -1;
    brf = new BrentRootFinder(*this);
}


const double MAxelet::rootThresh = 1.0e-10;

// Get an interpolation of the medial atom between 0 and 1
MAtom MAxelet::getInterpolation(double t) {
    static MAtom maError(0,0,0,0,0);

    // Boundaries must exist
    if(bLeft == NULL || bRight == NULL)
	return maError;

    // Sanity check
    if(t == 0 || t == 1) {
	return computeAtom(bLeft->getInterpolation(1.0-t),bRight->getInterpolation(t));
    }

    // Cache the t position
    ba_t = bRight->getInterpolation(t);

    // Use root finder
    double s = brf->findRoot(0,1,rootThresh);

    return (s == BrentRootFinder::ROOT_ERROR) ? 
	maError : computeAtom(ba_s,ba_t);
}

MAxelet::~MAxelet() {
    if(cache != NULL) {
	clearCache();
	delete cache;
    }
    delete brf;
}

class MAxeletIntersectFn : public Function1D {
    private:
	MAxelet *axelet;
	Vector2D t;
	double xt;
    public:
	MAxeletIntersectFn(MAxelet *a,const Vector2D &pos, const Vector2D &nrm) {
	    t = nrm.getNormal();
	    xt = t.dotProduct(pos);
	    axelet = a;
	}

	double evaluate(double s) {
	    // Find the atom at t
	    MAtom atm = axelet->getInterpolation(s);

	    // Find the distance to the line
	    return xt - t.dotProduct(atm.x()); 
	}
};

/* 
   This method finds the intersection of the medial axis with
   the specified vector.  
 */
MAtom MAxelet::interpolateAcross(const Vector2D &x,const Vector2D &n) {
    // Set up optimizer
    MAxeletIntersectFn f(this,x,n);
    BrentRootFinder brf2(f);

    // Search
    double t = brf2.findRoot(0,1,rootThresh);

    // Return the atom
    return getInterpolation(t);
}

// Get a boundary atom on the boundary.  Parameter t ranges from 0 to the size
// of the boundary.  Boundary atoms are interpolated for noninteger t and looked
// up for integer t.
BAtom Boundary::interpolate(double t) {
    dassert(t >= 0);
    int i;

    // Find the boundary positions
    t = fmod(t,m_brefs.size());
    i = (int) floor(t);
    t -= i;

    // Get the interpolation
    return getBoundlet(i).getInterpolation(t);
}

const BAtom &Boundary::estimate(double t) {
    int i;
    int n = m_brefs.size();

    // Find the boundary positions
    t = fmod(t,m_brefs.size());
    i = (int) floor(t);
    t -= i;

    // Get the interpolation
    return getBoundlet(i).getEstimation(t);
}

inline Boundary::Boundary(const list<BRef> &bList) {
    m_brefs.reserve(bList.size());
    m_blets.resize(bList.size());

    list<BRef>::const_iterator p;
    int t = 0;
    for(p=bList.begin();p!=bList.end();p++) {
	m_brefs.push_back(*p);
	p->node->setTShape(p->bIndex,this,t++);
    }
}

inline Boundary::Boundary() {
}

/**********************************************************************************
 * Medial Figure Class Implementation                                             *
 *********************************************************************************/
vector<double> MFigure::computeMedialParametrization() {
    int n = size();
    vector<double>p(n);

    p[0] = 0;
    for(int i=1;i<n;i++) {
	p[i] = p[i-1] + node(i)->x().distanceTo(node(i-1)->x());
    }

    return p;
}

MAxelet &MFigure::getMedialAxis(int link) {
    dassert(link >= 0 && link < size() - 1);

    // Get the right axis
    MAxelet &axis = node(link)->getSuccessorLink()->getAxelet();

    // Make sure it is valid
    if(node(link)->timestamp() > axis.m_tsStart || 
	    node(link+1)->timestamp() > axis.m_tsEnd) 
    {
	axis.compute(node(link)->getSuccessorLink());
    }

    // Call the estimate routine
    return axis;
}

MAtom MFigure::findCrossing(const Vector2D &x,const Vector2D &n) {
    // For each link, check for crossing with this line segment
    Vector2D t = n.getNormal();
    double dmin = 0;
    int imin = -1;

    Vector2D z0(n.dotProduct(x),t.dotProduct(x));

    for(int i=0;i<size()-1;i++) {
	Vector2D xa = node(i)->x();
	Vector2D xb = node(i+1)->x();
	Vector2D a = Vector2D(n.dotProduct(xa),t.dotProduct(xa)) - z0;
	Vector2D b = Vector2D(n.dotProduct(xb),t.dotProduct(xb)) - z0;

	if((a.y < 0 && b.y > 0) || (a.y > 0 && b.y < 0)) {
	    double d = fabs(((b.y-a.y)*b.x - (b.x-a.x)*b.y) / (b.y - a.y));
	    if(dmin == 0 || d < dmin) {
		dmin = d;
		imin = i;
	    }
	}
    }

    // Found a crossing
    if(imin >= 0) 
	return getMedialAxis(imin).interpolateAcross(x,n);

    // Error - no crossing
    return MAtom(0,0,0,0,0);
}

/**********************************************************************************
 * Medial Shape Class Implementation                                              *
 *********************************************************************************/
void MShape::addBoundary(const list<BRef> &bList) {
    Boundary *bnd = new Boundary(bList);
    m_bounds.push_back(bnd);
}

MShape::~MShape() {
    destroyList(m_bounds);
}

/**********************************************************************************
 * Medial Graph Class Implementation                                              *
 *********************************************************************************/
// Remove a curve link
void MGraph::removeCurveLink(CurveLink *link) {
    dassert(link);
    dassert(link->tail()->m_links[MNode::SUCC] == link);
    dassert(link->head()->m_links[MNode::PRED] == link);
    dassert(topologyEdit = true);

    link->tail()->m_links[MNode::SUCC] = NULL;
    link->head()->m_links[MNode::PRED] = NULL;

    m_links.remove(link);
    delete link;
}

// Add a curve link
void MGraph::addCurveLink(MNode *tail,MNode *head) {
    dassert(topologyEdit = true);
    dassert(tail && head);

    if(tail->m_links[MNode::SUCC])
	removeCurveLink((CurveLink *)tail->m_links[MNode::SUCC]);
    if(head->m_links[MNode::PRED])
	removeCurveLink((CurveLink *)head->m_links[MNode::PRED]);

    CurveLink *link = new CurveLink(this,tail,head);
    m_links.push_back(link);

    // Add link to its nodes
    tail->m_links[MNode::SUCC] = link;
    head->m_links[MNode::PRED] = link;
}

/**
 * Place an atom inside an existing figure
 */
MNode *MGraph::insertAtomIntoFigure(MFigure *f,int n,const MAtom &atom) {
    // This is a topology edit mode function
    dassert(topologyEdit == true);
    dassert(f->size() >= n && n >= 0);

    // Create a new node
    MNode *node = new MNode(atom,this,f);
    m_nodes.push_back(node);

    // We need to determine if there is a next and previous node
    if(n > 0) 
	addCurveLink(f->node(n-1),node);
    if(n < f->size())
	addCurveLink(f->node(n),node);

    return node;
}

/**
 * This method adds a figure from another graph to this graph
 */
void MGraph::addFigure(const vector<MAtom> &atoms,MNode **outNodeList) {
    dassert(topologyEdit == true);

    // Add each atom
    MNode *lastNode = NULL;
    for(int i=0;i<atoms.size();i++) {
	MNode *node = new MNode(atoms[i],this);
	m_nodes.push_back(node);

	if(i > 0)
	    addCurveLink(lastNode,node);
	lastNode = node;

	if(outNodeList) 
	    outNodeList[i] = node;
    }
}

/**
 * This method adds a figure from another graph to this graph
 */
void MGraph::addGroup(MNodeGroupIF &source) {
    dassert(topologyEdit == true);
    MNode *srcNode;

    // Go through all the source nodes; copy them into the new model
    for(srcNode=source.first();srcNode!=NULL;srcNode=source.next()) {
	// Create a copy of the node
	MNode *myNode = new MNode(*srcNode,this);

	// Add the node to the list of nodes
	m_nodes.push_back(myNode);

	// Mark the node to preserve the links
	srcNode->m_mark = (int)myNode;
    }

    // Go through all the source nodes again; This time, copy the links as well
    for(srcNode=source.first();srcNode!=NULL;srcNode=source.next()) {
	// Go through the list of links; skip CTF for now...
	for(int iLink=0;iLink<MNode::CTF_START;iLink++) {
	    // Get the next link
	    MLink *link = srcNode->m_links[iLink];
	    MLink *newLink = NULL;

	    // Only deal with non-null links that have marks on both ends
	    if(link && link->head()->m_mark && link->tail()->m_mark) {
		// Get to the copy of these nodes
		MNode *tail = (MNode *)(link->tail()->m_mark);
		MNode *head = (MNode *)(link->head()->m_mark);

		// Create a new link from this link
		newLink = link->makeCopyWithSub(this,tail,head);

		// Add the link
		m_links.push_back(newLink);
	    }

	    // Set the link for the current node
	    MNode *myNode = (MNode *)srcNode->m_mark;
	    myNode->m_links[iLink] = newLink;
	}
    }

    // Finally, clear all the marks.  This is critical!
    for(srcNode=source.first();srcNode!=NULL;srcNode=source.next()) {
	// Mark the node to preserve the links
	srcNode->m_mark = 0;
    }
}



/**
 * Remove all atoms
 */
void MGraph::removeAll() {
    destroyList(m_nodes);
    destroyList(m_links);
    destroyList(m_figures);
    destroyList(m_shapes);
}

void MGraph::beginTopologyEdit() {
    dassert(topologyEdit == false);
    topologyEdit = true;
}

void MGraph::endTopologyEdit() {
    dassert(topologyEdit = true);

    // Clear shapes and figures
    destroyList(m_figures);
    destroyList(m_shapes);

    // As we are building figures, we will also be building a list of BRefs that are on the 
    // boundary of each figure
    list <BRef> bndList;

    // Build figures by following links in the graph
    for(MNodeItr p=m_nodes.begin();p!=m_nodes.end();p++) {      
	MNode *node = *p;

	// For each node clear it's boundary's tShape value
	node->clearTShape();

	// Add all successors to the figure
	if(!node->hasPredecessor()) {
	    // Found a starting node, create a figure
	    MFigure *f = new MFigure(this);

	    bndList.push_back(BRef(node,MAtom::TAIL));
	    while(node) {            
		f->append(node);

		bndList.push_back(BRef(node,MAtom::LEFT));
		bndList.push_back(BRef(node,MAtom::RIGHT));
		if(!node->hasSuccessor())
		    bndList.push_back(BRef(node,MAtom::HEAD));

		node->m_figure = f;
		node = node->getSuccessor();
	    }			

	    // Add figure to the object
	    m_figures.push_back(f);
	}
    }

    // Start by searching at branches.  Branches are guaranteed to connect pieces of 
    // boundary 
    BRefItr bref;
    for(bref=bndList.begin();bref!=bndList.end();bref++) {      
	// Create a list into which we will dump atoms until we make a complete loop
	list<BRef> chain;

	// Our bref iterator
	BRef b = *bref;

	// Start out with the shape being NULL, we might hit a shape later
	MShape *shape = b.node->figure()->shape();

	// Marching direction
	bool cw = true;

	while(!b.marked()) {

	    // Set the traversal direction in the atom
	    b.node->m_ba[b.bIndex].bndFlipped = !cw;

	    chain.push_back(b);
	    b.mark(true);

	    // See if there is a branch here
	    BranchLink *bl = b.node->getBranch(b.bIndex);
	    if(bl) {
		// We found a branch at the atom.  The branch is a semaphore, making us
		// go from tail to head, but if end up at the tail end of one, we must
		// drop the whole chain
		if(bl->head() == b.node) {
		    break;	
		}
		else {
		    dassert(bl->tail() == b.node);

		    // Move to the head node
		    b = bl->m_bRef[1];

		    // Check if this is a known shape
		    shape = (shape) ? shape : b.node->figure()->shape();

		    // See if we need to flip directions
		    if(bl->m_flip)
			cw = !cw;

		    // Set the traversal direction in the atom
		    b.node->m_ba[b.bIndex].bndFlipped = !cw;

		    chain.push_back(b);
		    b.mark(true);			
		}
	    }

	    // Follow the figure
	    b = b.next(cw);

	} // end while not marked

	// We either ran into an boundary atom that's already processed or into a backwards
	// branch.  If the chain is not a loop drop it
	if(chain.size() > 0 && chain.front().equals(b)) {
	    if(shape == NULL) {
		shape = new MShape(this);
		m_shapes.push_back(shape);
	    }

	    shape->addBoundary(chain);

	    // Each figure needs to be told who it's shape is...
	    for(BRefItr chnref = chain.begin();chnref!=chain.end();chnref++) {      
		chnref->node->figure()->m_shape = shape;
	    }
	}
    }

    // Clear all markings
    for(bref=bndList.begin();bref!=bndList.end();bref++) {      
	bref->mark(false);
    }

    topologyEdit = false;
    m_selectionDirty = true;
}

/**
 * This method loads the graph from an old-fashioned registry file
 */
void MGraph::loadFromModelFile(Registry &folder) {
    // Start editing topology
    beginTopologyEdit();

    // First of all, empty the model
    removeAll();

    // Read the number of figures
    int n = folder.getIntValue("figureCount",0);

    // We will need to maintain an index array
    MNode ***index = NULL;
    int *figSize = NULL;
    if(n > 0) {
	index = new MNode**[n];
	figSize = new int[n];
    }

    // Read each figure
    for (int f=0; f<n; f++) {
	// Let's work with a subfolder
	Registry &fFolder = folder.getSubFolder("figure[%d]",f);

	// Get the number of atoms
	figSize[f] = fFolder.getIntValue("primitiveCount",0);
	index[f] = new MNode*[figSize[f]];

	// Create a vector for the atoms we read in
	vector <MAtom> atoms;   
	atoms.reserve(figSize[f]);

	// Read rho (if present)
	static const double DEFAULT_RHO = 0.25;
	double rho = fFolder.getDoubleValue("rho",DEFAULT_RHO,f);

	// Read in each atom
	for(int j=0;j<figSize[f];j++) {
	    MAtom a;
	    a.read(fFolder.getSubFolder("primitive[%d]",j));

	    // Update the rho at the atom level
	    if(a.rho() == DEFAULT_RHO)
		a.rho(rho);

	    atoms.push_back(a);
	}

	// Add the figure to the graph
	addFigure(atoms,index[f]);
    }

    // Read in the number of branches
    int nBranch = folder.getIntValue("branchCount",0);

    // Read each branch
    for(int b=0;b<nBranch;b++) {
	// Let's work with a subfolder
	Registry &bFolder = folder.getSubFolder("branch[%d]",b);

	// Some fields we need to read
	BRef bref[2];

	// Read head/tail information about this branch
	bool valid = true;
	for(int j=0;j<2;j++) {
	    // Read all the indices
	    int fId = bFolder.getIntValue("figure[%d]",-1,j);
	    int pId = bFolder.getIntValue("primitive[%d]",-1,j);
	    int bId = bFolder.getIntValue("boundary[%d]",-1,j);

	    // Make sure this is a reference to a valid atom
	    if(fId < 0 || pId < 0 || bId < 0 || 
		    fId >= n || pId >= figSize[fId] || bId > 1) {  
		valid = false;
		break;
	    }

	    bref[j].node = index[fId][pId];
	    bref[j].bIndex = bId + MAtom::LEFT;
	}

	bool flip = bFolder.getBooleanValue("flip",false);

	// Proceed only if a valid set is read
	if(valid) {
	    // Add the branch
	    addBranch(bref[0],bref[1],flip);
	}
    }

    // Clear some stuff
    if(n > 0) {
	for(int i=0;i<n;i++)
	    if(figSize[i] > 0)
		delete[] index[i];
	delete[] index;
	delete[] figSize;
    }

    // Rebuild topology
    endTopologyEdit();
}

const MAtom &MGraph::removeNode(MNode *node) {
    int i;
    for(i=0;i<MNode::BRANCH;i++)
	if(node->m_links[i] != NULL) 
	    removeCurveLink((CurveLink*)node->m_links[i]);
    for(i=MNode::BRANCH;i<MNode::CTF_START;i++)
	if(node->m_links[i] != NULL) 
	    removeBranchLink(((BranchLink*)node->m_links[i])->tailBRef());
    m_nodes.remove(node);
    return *node;
}

// Remove group of atoms
void MGraph::removeGroup(MNodeGroupIF &group) {
    dassert(topologyEdit = true);
    for(MNode *node = group.first();node!=NULL;node=group.next()) {
	removeNode(node);
    }
}

// Remove a branch
void MGraph::removeBranchLink(const BRef &bref) {
    dassert(topologyEdit = true);

    int spot = MNode::BRANCH + bref.bIndex; 
    BranchLink *link = bref.node->getBranch(bref.bIndex);
    if(link) {
	m_links.remove(link);
	for(int i=0;i<2;i++) {
	    // Clear pointer to this link from both ends of the branch
	    link->m_bRef[i].node->m_links[MNode::BRANCH+link->m_bRef[i].bIndex] = NULL;
	}
	delete link;
    }
}

// Add a branch relationship to the object
void MGraph::addBranch(const BRef &tail,const BRef &head,bool flip) {   
    dassert(topologyEdit = true);

    // Adding a branch.  First create and add a link to the link array
    BranchLink *link = new BranchLink(this,tail,head,flip);
    m_links.push_back(link);

    // Remove existing links at these positions
    if(tail.node->hasBranch(tail.bIndex)) {
	removeBranchLink(tail);
    }
    if(head.node->hasBranch(head.bIndex)) {
	removeBranchLink(head);
    }

    // Add the link reference to the atoms
    tail.node->m_links[MNode::BRANCH + tail.bIndex] = link;
    head.node->m_links[MNode::BRANCH + head.bIndex] = link;
}

// Add a branch relationship to the object
void MGraph::removeBranch(const BRef &tail) {   
    dassert(topologyEdit = true);

    // Remove existing links at these positions
    if(tail.node->hasBranch(tail.bIndex)) {
	removeBranchLink(tail);
    }
}

// Copy constructor
MGraph& MGraph::operator =(const MGraph &graph) {
    removeAll();
    beginTopologyEdit();

    // Create new nodes
    for(MNodeItr n=graph.m_nodes.begin();n!=graph.m_nodes.end();n++) {
	// Create and store a copy node
	MNode *srcNode = *n;
	MNode *node = new MNode(*srcNode,this);
	m_nodes.push_back(node);

	// Use mark to point to copy of the source node
	srcNode->m_mark = (int)node;
    }

    // Copy links
    for(MLinkItr l=graph.m_links.begin();l!=graph.m_links.end();l++) {
	MLink *link = *l;
	if(link->type()==MLink::CURVE) {
	    addCurveLink((MNode *)link->tail()->m_mark,(MNode *)link->head()->m_mark);
	}
	else if(link->type()==MLink::BRANCH) {
	    BranchLink *bl = (BranchLink *)link;
	    BRef tail((MNode*)bl->m_bRef[0].node->m_mark,bl->m_bRef[0].bIndex);
	    BRef head((MNode*)bl->m_bRef[1].node->m_mark,bl->m_bRef[1].bIndex);
	    addBranch(tail,head,bl->m_flip);
	}
    }

    // Clear marks
    graph.clearMarks();

    // Compute shapes and figures
    endTopologyEdit();

    return *this;
}

void MGraph::clearMarks() const {
    for(MNodeItr n=m_nodes.begin();n!=m_nodes.end();n++) {
	(*n)->m_mark = 0;
    }
}

void MGraph::writeXML(const char *fname) {
    FILE *f = fopen(fname,"wt");

    fprintf(f,"<?xml version=\"1.0\"?>\n");
    fprintf(f,"<!--DOCTYPE MRep URI \"http://www.cs.unc.edu/~pauly/mreps.dtd\" -->\n");
    fprintf(f,"<MRep>\n");
    fprintf(f,"   <Nodes>\n");

    int idxNode = 0;
    for(MNodeItr iNode=m_nodes.begin();iNode!=m_nodes.end();iNode++) {
	// Get a node
	MNode *node = *iNode;
	node->m_mark = ++idxNode;

	// Write the node out 
	fprintf(f,"      <Node index=\"%d\">\n",node->m_mark);
	fprintf(f,"         <x>%g</x>\n",node->x().x);
	fprintf(f,"         <y>%g</y>\n",node->x().y);
	fprintf(f,"         <r>%g</r>\n",node->r());
	fprintf(f,"         <fa>%g</fa>\n",node->faDeg());
	fprintf(f,"         <oa>%g</oa>\n",node->oaDeg());
	fprintf(f,"         <rho>%g</rho>\n",node->rho());
	fprintf(f,"         <endRC>%g</endRC>\n",node->endRC());
	fprintf(f,"      </Node>\n");
    }

    fprintf(f,"   </Nodes>\n");

    fprintf(f,"   <Links>\n");

    for(MLinkItr iLink=m_links.begin();iLink!=m_links.end();iLink++) {
	// Get a link
	MLink *link = *iLink;

	// Need an array of strings for boundary indices
	char *bIdxName[] = {"middle","left","right","head","tail"};

	// Write the link out
	fprintf(f,"      <Link type=\"%s\">\n",(link->type()==MLink::CURVE) ? "curve" : "branch");
	fprintf(f,"         <head>%d</head>\n",link->head()->m_mark);
	fprintf(f,"         <tail>%d</tail>\n",link->tail()->m_mark);
	if(link->type()==MLink::BRANCH) {
	    BranchLink *bl = (BranchLink *)link;
	    fprintf(f,"         <hbi>%s</hbi>\n",bIdxName[bl->headBRef().bIndex]);
	    fprintf(f,"         <tbi>%s</tbi>\n",bIdxName[bl->tailBRef().bIndex]);
	    fprintf(f,"         <flipped>%d</flipped>\n",(int)bl->flipped());
	}
	fprintf(f,"      </Link>\n");
    }

    fprintf(f,"   </Links>\n");
    fprintf(f,"</MRep>\n");
    fclose(f);

    clearMarks();
}

void MGraph::writeToModelFile(Registry &folder) {
    // Write all figures
    folder.setIntValue("figureCount",m_figures.size());

    int fIdx = 0;
    for(MFigureItr f=m_figures.begin();f!=m_figures.end();f++) {
	MFigure *figure = *f;

	// Open a subfolder for the figure
	Registry &ff = folder.getSubFolder("figure[%d]",fIdx++);

	// Write figure size
	ff.setIntValue("primitiveCount",figure->size());

	// Check if all the rho's are the same
	bool writeEachRho = false;
	int aIdx;
	for(aIdx=1;aIdx<figure->size();aIdx++) {
	    if(figure->node(aIdx)->rho() != figure->node(0)->rho()) {
		writeEachRho = true;
		break;
	    }
	}

	// Write the figure's rho
	if(!writeEachRho && figure->size()) {
	    ff.setDoubleValue("rho",figure->node(0)->rho());
	}

	// Write each atom in the figure
	for(aIdx=0;aIdx<figure->size();aIdx++) {
	    Registry &af = ff.getSubFolder("primitive[%d]",aIdx);
	    figure->node(aIdx)->write(af);

	    if(writeEachRho) 
		ff.setDoubleValue("primitive[%d].rho",figure->node(aIdx)->rho(),aIdx);

	    // Mark the node with figure/atom IDX so that we can write out branch links later
	    int mark = ((fIdx-1) << 16) + aIdx;
	    figure->node(aIdx)->m_mark = mark;
	}
    }

    // Write all branch links
    int bIdx = 0;
    for(MLinkItr l=m_links.begin();l!=m_links.end();l++) {
	if((*l)->type()==MLink::BRANCH) {
	    BranchLink *bl = (BranchLink *)(*l);
	    Registry &bf = folder.getSubFolder("branch[%d]",bIdx++);
	    for(int i=0;i<2;i++) {
		BRef bref = bl->m_bRef[i];
		bf.setIntValue("figure[%d]",bref.node->m_mark >> 16,i);
		bf.setIntValue("primitive[%d]",bref.node->m_mark & 0xffff,i);
		bf.setIntValue("boundary[%d]",bref.bIndex - MAtom::LEFT,i);
	    }
	    bf.setBooleanValue("flip",bl->m_flip);
	}
    }

    // Clean up all the marks
    clearMarks();

    folder.setIntValue("branchCount",bIdx);
}

void MGraph::writeToMatlabFile(char *fname) {


}

MGraph::~MGraph() {
    // Delete all the nodes
    for(MNodeItr n=m_nodes.begin();n!=m_nodes.end();n++) {
	MNode *node = *n;
	delete node;
    }

    // Delete all the links
    for(MLinkItr l=m_links.begin();l!=m_links.end();l++) {
	MLink *link = *l;
	delete link;
    }

    // Delete all the figures 
    for(MFigureItr f=m_figures.begin();f!=m_figures.end();f++) {
	MFigure *figure = *f;
	delete figure;
    }

    // Delete all the shapes
    for(MShapeItr s=m_shapes.begin();s!=m_shapes.end();s++) {
	MShape *shape = *s;
	delete shape;
    }

    //cout << "Deleted Graph\n";
}

MGraph* MGraph::createSaucage(int nAtoms,double r,double span) {
    MGraph *mg = new MGraph();
    vector <MAtom> figure(nAtoms);
    for(int i=0;i<nAtoms;i++) {
	figure[i].x(Vector2D(0,i*span));
	figure[i].r(r);
	figure[i].faDeg(-90);
	figure[i].oaDeg(90);
	if(i==0 || i==nAtoms-1)
	    figure[i].endRC(1.2);
    }

    mg->beginTopologyEdit();
    mg->addFigure(figure);
    mg->endTopologyEdit();
    return mg;
}

// Print out contents of a model
//void MGraph::print(ostream &s) const {
//}
