#ifndef MREPS_PAULY_2D
#define MREPS_PAULY_2D

/*******************************************************************
* Medial Representations Data Structures
*******************************************************************

Underdanding Data Structures for Medial Representations
-------------------------------------------------------

An MRep is a collection of medial atoms.  Each medial atom has primary
attributes, such as position, orientation, radius, object angle and some
others.  Each medial atom implies two or more boundary atoms, each with a
position, orientation and aperture size.

Within an MRep, atoms are grouped logically.  Two main units of organization 
are Figure and Shape.  Figure is a chain of atoms, sampled from a continuous
medial curve.  Shape is a collection of figures that are connected by pieces of
boundary.  Figures within a shape are attached to each other by defining links
between pairs of boundary atoms.  

A good way to envision the organization of an MRep is by looking at neighbour 
relationships between atoms.  Atoms may have intrafigure relationship - an atom
may have successor, a predecessor or both.  Atoms between figures have neighbour
relationships if their boundary atoms are used to connect the figures.  Atoms between
shapes may also have neighbor relationships if one shape is a refinement of another.

There are two major aspects of an object.  One is geometrical - positions and orientations
of the atoms.  The second is topological - defining the organization of atoms into figures and
shapes, defining neighbour relationships between atoms.  Topological changes are considered 'rare'.
By rare we mean that they are not performed continuously during optimization, but that instead they
are always responses to user commands.  Geometrical changes, however, are common and should be very fast.

Atoms within an object can be accessed via the group mechanism.  A group is a collection of references to 
atoms within the object.

Topological changes are performed on an object as a whole.  Once a topological change has been made, the 
object enters a 'volatile' state.  Many things depending on the object may need to be recalculated. 

*/ 


#include <vector>
#include <list>
// #include <iostream>
#include <string>
#include <assert.h>
#include <math.h>
#include <registry.h>
#include <vector2d.h>
#include <optima.h>
#include <BrentLinearMethod.h>
#include <phquintic.h>

using namespace std;

class BAtom;
class MAtom;
class MFigure;
class MShape;
class Boundary;
class MNode;
class MNodeGroupIF;
class MNodeGroup;
class MLink;
class CurveLink;
class BranchLink;
class CtfLink;
class MGraph;
class MAxelet;
class MRepBoundlet;


/**
* Boundary Atom Class 
*
* A boundary atom.  Represents a boundary position, orientation, apertuire and weight tuple
* Atoms represent positions at sides and ends of medial atoms as well as interpolated boundary
* positions.
*/
class BAtom {
public:
	// Position, Normal to the boundary (expected to be of unit length)
	Vector2D x,n;

	// Slope of the vector n; Scale at which boundary is measured (r-proportional);
	// Polarity / weight combination.
	double slope,scale,polarity;

	// Boundary reference - this is very crappy, should not be part of BAtom, there 
	// should be a BNode, etc,etc,etc...
	Boundary *bnd;
	int tShape;

	// Direction in which the boundary is being tranversed at tis atom.  Default is true, false
	// is for atoms that are indentations
	bool bndFlipped;

	// Default value of the above
	static const int NOT_IN_SHAPE;

	// Return the tangent vector
	Vector2D tangent() const {
		return (bndFlipped) ? n.getNormal() : -n.getNormal();
	}

	// Constructors
	BAtom(const Vector2D &inX,const Vector2D &inN,double inScale,double inPolarity) :
	x(inX),n(inN),scale(inScale),polarity(inPolarity),tShape(NOT_IN_SHAPE),bnd(NULL),bndFlipped(false) {};

	BAtom() :
	x(0,0,1),n(0,1,0),scale(0.0),polarity(1.0),tShape(NOT_IN_SHAPE),bnd(NULL),bndFlipped(false) {};

	virtual ~BAtom() {
		// cout << "a";
	}
};

/**
* Medial Atom Class
*
* This class encapsulates information about a medial atom.  It has primary attributes (position,orientation,
* object angle, radius, rho (radius multiplier to get aperture size), endRC (radius multipler to get distance 
* to the end BAtom2D).  The primary attributes are persistant.
* There are also figure-set attributes.  These are flags that determine whether the primitive is selected, or is at one
* of the ends of the figure.  These attributes are maintained by the figure class.
* Finally there are secondary attributes.  These are computed/cahced on demand by the MAtom2D class.
*/
class MAtom {
protected:
	// A collection of BAtoms.  ba[MIDDLE] holds a couple of primary attributes
	BAtom m_ba[5];

	// Object angle; Radius of the atom; Rho; End radius of curvature
	double m_oa,m_r,m_rho,m_endRC;

	// Timestamp of the primary data in this atom.  Things that depend on this atom will check their
	// timestamp versus this timestamp and update itself if nessesary
	int m_tsPrimary,m_tsSecondary;

	// Selection mask of this atom.  Selection is used for grouping atoms
	bool m_selected;

	// This method is called whenever a secondary attribute is read.  This has to be inline.
	// We cast this to a non constant pointer to call the compute method within a const method.
	void prepareSecondary() const {
		if(m_tsSecondary < m_tsPrimary) {
			MAtom *cheat = (MAtom *)this;
			cheat->computeSecondary();
		}
	}

	// This actually computes the secondary attributes
	void computeSecondary();

public:
	// These are the labels for various BAtoms associated with this MAtom2D
	static const int MIDDLE,LEFT,RIGHT,HEAD,TAIL;

	// This constructor makes an atom given values for the five primary attributes
	MAtom(double x=0.5, double y=0.5, double r=1.0, double taDeg=90, double toDeg=90) 
	{
		m_ba[MIDDLE].x.x = x;
		m_ba[MIDDLE].x.y = y;
		m_ba[MIDDLE].slope = M_PI*taDeg/180;
		m_oa = M_PI*toDeg/180;
		m_r = r;
		m_rho = 0.25;
		m_endRC = 1.0;
		m_tsPrimary = m_tsSecondary = 0;
		m_selected = false;
	}

	virtual ~MAtom() {
		// cout << "A";
	};

	// Print atom for debugging
	// void print(ostream &s) const;

	// Set the position
	void x(const Vector2D &x) {
		m_ba[MIDDLE].x = x;
		m_tsPrimary++;
	}

	// Set the radius
	void r(double r) {
		m_r = r;
		m_tsPrimary++;
	}

	// Set the axis angle
	void faDeg(double faDeg) {
		fa(faDeg * M_PI / 180);
	}

	// Set the primitive object angle
	void oaDeg(double oaDeg) {
		oa(oaDeg*M_PI/180.0);
	}

	// Set the axis angle
	void fa(double faRad) {
		// Place axial angle in -180 to 180 range
		faRad = (faRad >= -M_PI) ? fmod(faRad+M_PI,2*M_PI)-M_PI : fmod(faRad+M_PI,2*M_PI)+M_PI;
		m_ba[MIDDLE].slope = faRad;
		m_tsPrimary++;
	}

	// Set the primitive object angle
	void oa(double oaRad) {
		m_oa = acos(cos(oaRad));
		m_tsPrimary++;
	}

	// Set the polarity.  The second parameter selects the boundary atom 
	// whose position is being affected.  (Allowed values are LEFT,RIGHT,MIDDLE and END, with
	// MIDDLE setting all three)
	void polarity(double value,int which) {
		m_ba[which].polarity = value;

		if(which==MIDDLE) {
			polarity(value,LEFT);
			polarity(value,RIGHT);
			polarity(value,HEAD);
			polarity(value,TAIL);
		}

		m_tsPrimary++;
	}

	// Set the end radius of curvature
	void endRC(double endRadiusOfCurvatureMultiplier) {
		m_endRC = endRadiusOfCurvatureMultiplier;
		m_tsPrimary++;
	}

	// Set rho.  Notice that when primitive is added to a figure, rho is overridden
	void rho(double rho) {
		m_rho = rho;
		m_tsPrimary++;
	}

	// Get position of primitive or boundary atoms
	const Vector2D& x(int which = MIDDLE) const {
		prepareSecondary();
		return m_ba[which].x;
	}

	// Return the radius
	double r() const {
		return m_r;
	}

	// Returns axis angle in degrees
	double faDeg() const {
		return m_ba[MIDDLE].slope*180.0/M_PI;
	}	

	// Returns object angle in degrees
	double oaDeg() const {
		return m_oa*180.0/M_PI;
	}	

	// Returns axis angle in degrees
	double fa() const {
		return m_ba[MIDDLE].slope;
	}	

	// Returns object angle in degrees
	double oa() const {
		return m_oa;
	}	

	// Return the polarity.  With default parameter (MIDDLE) it will return the polarity for the
	// whole primitive.  This polarity may or may not equal the polarities of the boundary atoms
	double polarity(int which) const {
		return m_ba[which].polarity;
	}

	// Get rho
	double rho() const {
		return m_rho;
	}

	// Get a reference to the entire boundary atom at giver index
	const BAtom &bAtom(int which) const {
		prepareSecondary();      
		return m_ba[which];
	}

	// Get the normal vector.  For MIDDLE it returns the frame vector and for all others it returns the
	// boundary normal
	const Vector2D &n(int which=MIDDLE) const {
		prepareSecondary();
		return m_ba[which].n;
	}

	// Returns the slope of the vector n.
	double slope(int which=MIDDLE) const {
		prepareSecondary();
		return m_ba[which].slope;
	}

	// Get the aperture scale associated with BAtom2D.  By default returns scale of the MAtom2D, (rho*r)
	double scale(int which) const {
		prepareSecondary();
		return m_ba[which].scale;
	}

	// Return the end-cap distance
	double endRC() const {
		return m_endRC;
	}

	// Select a atom using a selection mask.  Whatever bits are set in the mask will be set in the atom.
	virtual void select()	{
		m_selected = true;
	}

	// Unselect a atom using a selection mask.  Whatever bits are set in the mask will be cleared in the atom.
	virtual void deselect(short mask = 0x0001) {
		m_selected = false;
	}

	// Toggle selection of the atom using a mask.  The bit positions where mask is 1 will be flipped
	virtual void toggle(short mask = 0x0001)	{
		m_selected = !m_selected;
	}

	// Returns true if the ones in the mask are also ones in the atom's mask
	bool selected() const {
		return m_selected;
	}

	// Check if the primitive is inside a bounding box
	bool inBoundingBox(const Vector2D &p1, const Vector2D &p2) const;
	bool inBoundingBox(double x0, double y0, double x1, double y1) const;

	// Transformation functions. Scaling and rotation are about medial site location.
	void translateBy(const Vector2D &v)	{
		m_ba[MIDDLE].x += v;
		m_tsPrimary++;
	}

	void translateBy(double dx, double dy)	{
		m_ba[MIDDLE].x.x += dx;
		m_ba[MIDDLE].x.y += dy;
		m_tsPrimary++;
	}

	void rotateBy(double angleDeg) {
		faDeg(faDeg()+angleDeg);
	}

	void scaleBy(double mag) {
		m_r *= mag;
		m_tsPrimary++;
	}

	void changeAngleBy(double angleDeg) {
		oaDeg(oaDeg()+angleDeg);
	}

	int timestamp() const {
		return m_tsPrimary;
	}

	// This method backs up the primary characteristics of a medial atom.  Not all
	// characteristics are backed up: only position, orientation, angles and endRC.  
	// We do not backup rho or polarities
	Vector getData(bool endRC);

	// This method reads data from a vector, with end data if necessary and returns the
	// new position on the vector after it is done
	int setData(const Vector &v,bool endRC,int startIdx = 0);

	// Read and write to/from regisry folder
	void read(Registry &folder);
	void write(Registry &folder);
};

/**
* A Node is a wrapper around the atom that is aware of the graph to which it belongs and can return
* neighbour information.  MAtom is purely a geometrical construct while MNode is part of a large data structure
*
* MNodes can not be created other than by the MGraph object.  MAtoms however can exist and be created on their own
*/
class MNode : public MAtom {
private:
	// Collection of links in this atom.  The links are held in a uniform array.  The first two links are the 
	// pointers to CurveLinks (or NULLs) and the second two are pointers to branch links or NULLs.  The rest of
	// the links are arbitrary
	vector<MLink *> m_links;

	// Pointers to the figure in which the atom, the shape and the level of detail
	MFigure *m_figure;

	// Reference to the graph containing this node
	MGraph *m_graph;

	// Mark - this value is used when nodes are tranversed by methods such as save, or build.  At the end of 
	// such a method the marks are always guaranteed to be set to zero
	int m_mark;

	// Initialize a node
	MNode(const MAtom &atom,MGraph *graph,MFigure *figure=NULL) : MAtom(atom),m_links(CTF_START,(MLink*)NULL) {
		m_figure = figure;
		m_graph = graph;
		m_mark = 0;
	}

	// Dummy initializer
	MNode() : m_links(CTF_START,(MLink*)NULL) {
		m_graph = NULL;
		m_figure = NULL;
		m_mark = 0;
	}

	// Clear tShape in all boundary atoms.  This will have to be enhanced later
	void clearTShape();
	void setTShape(int which,Boundary *bnd,int t) {
		m_ba[which].tShape = t;
		m_ba[which].bnd = bnd;
	}

public:
	static const int SUCC,PRED,BRANCH,CTF_START;

	bool hasSuccessor(){
		return (m_links[SUCC] != NULL);
	}
	bool hasPredecessor(){
		return (m_links[PRED] != NULL);
	} 
	bool isEnd() {
		return !(hasSuccessor() && hasPredecessor());
	}

	MNode *getSuccessor();
	MNode *getPredecessor();

	bool hasBranch(int branchSide) {
		return (m_links[MNode::BRANCH + branchSide] != NULL);
	}

	// Return a pointer to the neighbor that is the left/right branch off this atom or NULL if one does not
	// exist.  Parameter branchSide takes value MAtom::RIGHT or MAtom::LEFT.  Return parameter nbrBranchSide 
	// also takes value RIGHT or LEFT and refers to the side of the neighbor that is linked with this atom
	BranchLink *getBranch(int branchSide) {
		return (BranchLink*)m_links[MNode::BRANCH + branchSide];
	}

	// Return a predecessor link
	CurveLink *getPredecessorLink() {
		return (CurveLink*)m_links[PRED];
	}

	// Return a successor link
	CurveLink *getSuccessorLink() {
		return (CurveLink*)m_links[SUCC];
	}

	// Get a group of all of the neighbors of this atom
	MNodeGroup getNeighbors();

	// Get the reference to the figure containing this atom
	MFigure *figure() {
		return m_figure;
	}

	virtual ~MNode() {
		// cout << "N";
	}

	// Selection methods are overriden to notify the parent graph
	virtual void select();
	virtual void deselect();
	virtual void toggle();

	friend class MGraph;
	friend class MLink;
	friend class BRef;
	friend class Boundary;
};
typedef list<MNode*>::const_iterator MNodeItr;

/**
* A reference to a Boundary Atom
*/
class BRef {
public:
	// Index of the medial atom in the graph
	MNode *node;

	// Index of the boundary atom in the medial atom
	int bIndex;

	// This method returns the next BRef in clockwise direction
	BRef next(bool clockWise);

	// This method sets a mark in the corresponding node
	void mark(bool onoff) {
		int m = 1 << bIndex;
		if(onoff)
			node->m_mark |= m;
		else
			node->m_mark &= (!m);
	}

	bool marked() {
		int m = 1 << bIndex;
		return (node->m_mark & m) != 0;
	}

	// Check if it's same as another bref
	bool equals(const BRef &bref) const;

	// Constructor
	BRef(MNode *node,int bIndex) {
		this->node = node;
		this->bIndex = bIndex;
	}

	// Default constructor
	BRef() {
		node = NULL;
		bIndex = -1;
	};

	virtual ~BRef() {
		// cout << "r";
	}

	// Return a pointer to the BAtom that we refer to or NULL of one is not found
	const BAtom *bAtom() const {
		return (node==NULL) ? NULL : &node->bAtom(bIndex);
	}
};


/**
* Atom Group Interface
*
* This is an ordered collection of atoms to which geometric operations can be applied.
* The group is a basic construct for dealing with atoms.  Groups do not contain atoms,
* instead they have references to the atoms contained within the source object.  Manipulating
* an atom in a group changes its value in the object.
*
* Deleting an atom from a group does not remove it from an object.  The topology of the object
* can not be affected by manipulating group object returned by the object.  To manipulate those
* use the topology interface of the object class.
*
* Notice the gg_ methods in this class.  The abbreviation 'gg' stands for 'Get Group'
*
* Selection is handled trough masks.  That means there may be up to 16 levels of selection.  This allows us
* to group things at multiple levels.  For example we could be working with an object that has mreps at different
* levels of detail.  We could isolate a level of detail by applying a mask to atoms that we are interested with
* and then limiting all operations to the group matching that mask.
*/
class MNodeGroupIF {
public:
	// These methods are used internally for parsing through the nodes.  These methods 
	// are implemented by subclasses in a way depending on type of containers used.  
	// This allows vector-based figures and shapes and list-based objects to share a
	// common interface.  The methods return NULL once the list has been parsed
	// CODE:
	// for(MNode *n = first();n != NULL;n = next())
	virtual MNode* first() = 0;
	virtual MNode* next() = 0;
	virtual MGraph* graph() = 0;

	// Select all the atoms in this group
	virtual void select() {
		for(MNode *n = first();n!=NULL;n = next())
			n->select();
	}

	// Deselect all the atoms in the group
	virtual void deselect() {
		for(MNode *n = first();n!=NULL;n = next())
			n->deselect();
	}

	// Toggle the selection status in the group
	virtual void toggle() {
		for(MNode *n = first();n!=NULL;n = next())
			n->toggle();
	}

	// Get the group of atoms that are selected.
	virtual MNodeGroup ggSelected();

	// Get the number of selected nodes
	virtual int nSelected() {
		int count = 0;
		for(MNode *n = first();n!=NULL;n = next())
			if(n->selected())
				count++;
		return count;
	}

	// Get a group of atoms in this group inside a rectangular region
	MNodeGroup ggInBox(const Vector2D &min,const Vector2D &max);

	// Get atom closest to a point
	MNode *closestAtom(const Vector2D &x);

	// Get the extents of the figure
	void getExtents(Vector2D &extMin,Vector2D &extMax);

	// Get the center of gravity of the figure
	Vector2D getCOG();

	// Elongate and thin the group (scale the group but keep the radii fixed)
	void elongate(double scalefact, const Vector2D &centroid);
	void elongate(double scalefact) {
		elongate(scalefact,getCOG());
	}

	void thin(double scalefact, const Vector2D &centroid);
	void thin(double scalefact) {
		thin(scalefact,getCOG());
	}

	// Uniformly scale the group
	void scaleBy(double scalefact, const Vector2D &centroid);	
	void scaleBy(double scalefact) {
		scaleBy(scalefact,getCOG());
	}

	// Rotate the group
	void rotateBy(double angledeg, const Vector2D &centroid);
	void rotateBy(double angledeg) {
		rotateBy(angledeg,getCOG());
	}

	// Translate the group
	void translateBy(const Vector2D &dx);

	// Set the object angle of the whole group
	void objectAngle(double angledeg);
	void changeObjectAngleBy(double angledeg);

	// Set the rho for the whole group
	void rho(double rho);

	// Set the polarity of the whole group
	void polarity(double value,int which=MAtom::MIDDLE);

	// Scale the group to fit inside a region
	void scaleToFit(double x,double y,double width,double height);

	// Fill a vector of MNode* with nodes from this group.  Useful when we need an array
	virtual vector<MNode *> getArray();

	// Fill a vector of MAtom* with nodes from this group.  Useful when we need an array of atoms
	virtual vector<MAtom> getAtomArray();

	// Get all boundary atoms in this group that are in an actual boundary.  Return in a form
	// of a nx2 matrix
	Matrix getBoundaryPoints();

	// This method fits an ellipse to the group
	void fitEllipse(Vector2D &center,Vector2D &major,Vector2D &minor);

	// Destructor
	virtual ~MNodeGroupIF() {
		//  "g";
	}
};


/**
* Medial atom group
*
* An implementation of the atom group interface
*/
class MNodeGroup : public MNodeGroupIF {
private:
	int itr;

protected:
	// Collection of atom indices
	vector <MNode *> m_nodes;

	// Reference to a graph to shich this group belongs
	MGraph *m_graph;

	// Constructor
	MNodeGroup(MGraph *graph = NULL) : m_graph(graph) {
	}


public:

	// See doc. in NodeGroupIF.  Uses list iterators
	virtual MNode* first() {
		itr = 0;
		return ((int)m_nodes.size() > itr) ? m_nodes[itr] : NULL;
	}

	virtual MNode* next() {
		itr++;
		return ((int)m_nodes.size() > itr) ? m_nodes[itr] : NULL;
	}

	virtual MGraph* graph() {
		return m_graph;
	}

	// Add an atom to the group
	void append(MNode *node) {
		dassert(node!=NULL);
		m_nodes.push_back(node);
	}

	// Access n-th node in the group
	MNode *node(int i) {
		dassert(i >= 0 && i < m_nodes.size());
		return m_nodes[i];
	}

	int size() {
		return m_nodes.size();
	}

	vector<MNode *> getArray() {
		return m_nodes;
	}

	friend class MNode;
	friend class MGraph;
	friend class MNodeGroupIF;
};

const static int boundletCacheSize = 120;
class MRepBoundlet : public CBoundlet {
private:
	// Timestamp of this boundlet, used to update it when nesessary
	int m_tsStart,m_tsEnd;

	// Coefficients for the spline
	Vector ax,ay,bx,by;

	// Distance between the ends of the spline
	double kScale, bScale, polarity;

	// Starting and ending BAtoms
	// BAtom x0,x1;

	// A collection of precomputed atoms (120) for now.  I use 120 because it's such a nice number
	// OK so this is incredibly cheesy.  Too bad.
	BAtom cache[boundletCacheSize+1];

	// A static collection of precomputed t vectors.  This saves us from computing t^n
	// This should make boundlet computation lightning fast
	static Vector *tnCache;

	// Get fast interpolation (at a cached t position)
	void fastInterpolation(int i) {
		BAtom &b = cache[i];
		const Vector &t = tnCache[i];

		// Compute the normal of the cached boundary atom
		b.n.y = - (bx(0) + bx(1)*t(1) + bx(2)*t(2));
		b.n.x = by(0) + by(1)*t(1) + by(2)*t(2);

		// Now comes the problem - it should be normalized, which needs a square root
		b.n.normalize();

		// While square root is churning, compute the boundary position
		b.x.y = ay(0) + ay(1)*t(1) + ay(2)*t(2) + ay(3)*t(3);
		b.x.x = ax(0) + ax(1)*t(1) + ax(2)*t(2) + ax(3)*t(3);

		// Scale interpolation.
		b.scale = bScale + t(1) * kScale;
		b.polarity = polarity;		
	}

public:
	void compute(const BAtom &start,const BAtom &end);

	// Get an interpolation of the atom between 0 and 1
	BAtom getInterpolation(double t) const {
		BAtom b;

		b.x.x = (((ax(3)*t)+ax(2))*t+ax(1))*t + ax(0);
		b.x.y = (((ay(3)*t)+ay(2))*t+ay(1))*t + ay(0);
		b.n.x = (((by(2)*t)+by(1))*t+by(0));
		b.n.y = - (((bx(2)*t)+bx(1))*t+bx(0));
		b.n.normalize();

		b.scale = bScale + t * kScale;
		b.polarity = polarity;

		return b;
	}

	// Get the precomputed atom nearest to the specified t
	const BAtom &getEstimation(double t) {
		dassert(0 <= t && t <= 1);
		int idx = (int)(t*boundletCacheSize);

		// Compute the atom if necessary
		if(cache[idx].scale == 0.0) 
			fastInterpolation(idx);
		//cache[idx] = getInterpolation(tnCache[idx](1));

		// Return the cached atom
		return cache[idx];
	}

	MRepBoundlet();

	friend class Boundary;
	friend class MAxelet;
};



/**
* A representation of the interpolated medial axis between two medial atoms
*/
class MAxelet : public Function1D {
private:
	// A cache of interpolated medial atoms and positions
	MAtom **cache;

	// Two boundlets that are used to perform computation
	const MRepBoundlet *bLeft,*bRight;

	// A pointer to the curve link to which this is attached
	// CurveLink *link;

	// This deletes the cache
	void clearCache();

	// Compute the medialness function for an s,t pair
	double computeMedFun(double s,double t) {
		ba_t = bRight->getInterpolation(t);
		ba_s = bLeft->getInterpolation(s);
		return (ba_t.x - ba_s.x).dotProduct(ba_t.tangent() - ba_s.tangent());
	}

	// Compute the value of the medialness function for a fixed t
	double evaluate(double s) {
		ba_s = bLeft->getInterpolation(s);
		return (ba_t.x - ba_s.x).dotProduct(ba_t.tangent() - ba_s.tangent());
	}

	// Cached t value for above method
	BAtom ba_s,ba_t;

	// Compute the gradient of the medialness function
	Vector2D computeGradient(double s,double t);

	// Compute the atom at the s, t pair
	static MAtom computeAtom(const BAtom &s,const BAtom &t);

	// The root finder
	BrentRootFinder *brf;

	// The root finder threshold
	static const double rootThresh;

	// Is this a traceable medial axis
	bool traceable;

	// Timestamps of the left and right atoms
	int m_tsStart,m_tsEnd;

public:
	// This method computes all the necessary matrices and clears the cache.
	// The method takes a link as a parameter
	void compute(CurveLink *link);

	// Get an interpolation of the medial atom between 0 and 1
	MAtom getInterpolation(double t);

	// Get the precomputed atom nearest to the specified t
	const MAtom &getEstimation(double t) {
		dassert(0 <= t && t <= 1);
		int idx = (int)(t*boundletCacheSize);

		// Compute the cached medial atom
		if(cache[idx] == NULL) {
			double tInt = idx / ((double)boundletCacheSize);
			cache[idx] = new MAtom();
			*cache[idx] = getInterpolation(tInt);
		}

		return *cache[idx];
	}

	// Find crossing of the interpolated axis with a line defined by x,n
	MAtom interpolateAcross(const Vector2D &x,const Vector2D &n);

	// Is this traceable
	bool isTraceable() {
		return traceable;
	}

	// Constructor, destructor
	MAxelet();
	~MAxelet();

	friend class MFigure;
};






/**
* Medial Link Class
*
* Medial Link is a parent class for more descriptive links. 
*/
class MLink {
protected:
	// Reference to the graph to which this link belongs
	MGraph *m_graph;

	// Constructor
	MLink(MGraph *graph=NULL) {
		m_graph = graph;
	}

	// This is a weird method.  It creates a copy of a link, but uses substitution
	// for nodes involved in the nodes.  This should be overridden by other links
	virtual MLink *makeCopyWithSub(MGraph *newGraph,MNode *newTail,MNode *newHead) {
		return NULL;
	}


public:
	// This method returns a value that identifies the type of the link that this is
	virtual int type() const {
		return MLink::WRONG_TYPE;
	}

	// These are the different recognized link types
	static const int WRONG_TYPE,CURVE,BRANCH,CTF;

	// Return the reference to the 'head' of the link
	virtual MNode *head() {
		dassert(0);
		return NULL;
	}

	// Return the reference to the tail of the link
	virtual MNode *tail() {
		dassert(0);
		return NULL;
	}

	// Always need a virtual destructor
	virtual ~MLink() {
		// cout << "l";
	}

	friend class MNode;
	friend class MGraph;
};
typedef list<MLink*>::const_iterator MLinkItr;

/**
* Curve Link Structure
*
* A curve link is a link from a medial atom to the next/previuos medial atom on the 
* continuous medial curve.  Curve links are also known an Intrafigural links.  This is 
* a directed link with a head and a tail.
*/
class CurveLink : public MLink {
private:
	// Indices of the two atoms involved in the link
	MNode *m_node[2];

	// A medial interpolant between the atoms
	MAxelet axelet;

	// Empty constructor
	CurveLink() : MLink(NULL){
		m_node[0] = m_node[1] = NULL;
	}	

	// Constructor
	CurveLink(MGraph *graph,MNode *tail,MNode *head) : MLink(graph) {
		m_node[0] = tail;
		m_node[1] = head;
	}

	MLink *makeCopyWithSub(MGraph *newGraph,MNode *newTail,MNode *newHead);

public:
	// Identify this type of link
	virtual int type() const {
		return MLink::CURVE;
	}

	// Return the reference to the 'head' of the link
	virtual MNode *head() {
		return m_node[1];
	}

	// Return the reference to the tail of the link
	virtual MNode *tail() {
		return m_node[0];
	}

	// Get the axelet
	MAxelet &getAxelet() {
		return axelet;
	}

	friend class MNode;
	friend class MGraph;
};

/**
* Branch Link Structure
*
* A branch link is a link between atoms in a connected figure.
* continuous medial curve.  Curve links are also known an Interfigural links.
*/
class BranchLink : public MLink {   
private:
	// Indices within the medial atom of the boundary atoms involved in the relationship (should only
	// be MAtom2D::LEFT or MAtom2D::RIGHT
	BRef m_bRef[2];

	// This flag specifies whether the sense in with we trace the boundary changes at this
	// branch (cw to ccw) and is used to distinguish between protrusions and indentations
	bool m_flip;

	// Empty constructor
	BranchLink() : MLink(NULL) {}

	// Constructor
	BranchLink(MGraph *graph,MNode *atom1,MNode *atom2,int bIndex1,int bIndex2,bool flip);

	// Constructor
	BranchLink(MGraph *graph,const BRef &tail,const BRef &head,bool flip);

	// Copy operation
	MLink *makeCopyWithSub(MGraph *newGraph,MNode *newTail,MNode *newHead);

public:
	// Identify this type of link
	virtual int type() const {
		return MLink::BRANCH;
	};

	// Return the reference to the 'head' of the link
	virtual MNode *head() {
		return m_bRef[1].node;
	}

	// Return the reference to the tail of the link
	virtual MNode *tail() {
		return m_bRef[0].node;
	}

	// Return the reference to the 'head' BRef
	virtual const BRef &headBRef() {
		return m_bRef[1];
	}

	// Return the reference to the tail Bref
	virtual const BRef &tailBRef() {
		return m_bRef[0];
	}

	// Return whether this link is flipped
	bool flipped() {
		return m_flip;
	}

	friend class MNode;
	friend class MGraph;
};

/**
* Coarse To Fine Link Structure
*
* A coarse to fine link defines a relationship between an atom and its refinement.  
* Each atom may be a parent or child of multiple atoms
*/
class CtfLink : public MLink {
public:
	// Identify this type of link
	virtual int type() const {
		return MLink::CTF;
	};

	friend class MNode;
	friend class MGraph;
};

/*
template <class T>
class SubSpline {
static const int width;

vector<T> x;
vector <bool> computed;

// End points of the spline
static void buildParentArray();

// Recursive get procedure
const T &get(int i);
public:
SubSpline() : x(width+1),computed(width+1,false) {
}

virtual ~SubSpline() {
}

void setEnds(const T &xStart,const T &xEnd){
x[0] = xStart;
x[width] = xEnd;
computed.assign(width+1,false);
computed[0] = computed[width] = 1;
}

const T &getEstimation(double t) {
dassert(0.0 <= t && t <= 1.0);
int i = floor(width * t);
return get(i);
}

static int *computeLogArray() {
int *logArray = new int[width+1];
for(int i=1;i<=width;i++) {
int j = 0;
while(((i >> j) & 1) == 0)
j++;
logArray[i] = 1 << j;
}
return logArray;
}
};
*/


class Boundary {
	vector <BRef> m_brefs;
	vector <MRepBoundlet> m_blets;

	Boundary();
	Boundary(const list<BRef> &);

public:
	// Get a boundary atom on the boundary.  Parameter t ranges from 0 to the size
	// of the boundary.  Boundary atoms are interpolated for noninteger t and looked
	// up for integer t.
	BAtom interpolate(double t);

	// Get an approximation of the atom at position t.  The boundary is precached at 120 positions 
	// and this returns the nearest neighbour of one og these positions
	const BAtom &estimate(double t);

	// Return a reference to the atom
	const BRef &bRef(int idx) {
		dassert(idx >= 0 && idx < m_brefs.size());
		return m_brefs[idx];
	}

	// Return a boundlet that starts at a position idx, initialize if necessary
	MRepBoundlet &getBoundlet(int i) {
		int j = i+1;
		int n = m_brefs.size();
		i = (i < 0) ? n + (i % n) : i % n;
		j = (j < 0) ? n + (j % n) : j % n;

		// Check the timestamps
		if(m_brefs[i].node->timestamp() > m_blets[i].m_tsStart ||
			m_brefs[j].node->timestamp() > m_blets[i].m_tsEnd) 
		{
			// Recomputethe boundlets on timestamp failure
			m_blets[i].compute(m_brefs[i].node->bAtom(m_brefs[i].bIndex),
				m_brefs[j].node->bAtom(m_brefs[j].bIndex));
			m_blets[i].m_tsStart = m_brefs[i].node->timestamp();
			m_blets[i].m_tsEnd = m_brefs[j].node->timestamp();
		}

		return m_blets[i];
	}

	// Number of boundary atoms that form this boundary
	int size() {
		return m_brefs.size();
	}

	// Later we may provide read only access to internal lists...
	friend class MShape;
};
typedef list<Boundary *>::const_iterator  BoundaryItr;

/**
* A shape is a figure with a boundary description
*/
class MShape : public MNodeGroup {
private:
	// A list of boundary vertices in order
	list <Boundary*> m_bounds;

	// The index of the level of detail containing this shape
	int m_lod;

	// Constructor
	MShape(MGraph *graph=NULL) : MNodeGroup(graph) {
	}

	// Add a collection of boundary references to this shape
	void addBoundary(const list <BRef> &bList);

public:
	// Get the shape that contains this figure
	// MLevel &lod() {
	//    return m_graph.lod(m_lod);
	// }

	// Access the boundary list
	const list<Boundary *> &bounds() {
		return m_bounds;
	}

	virtual ~MShape();

	friend class MNode;
	friend class MGraph;
};
typedef list<MShape*>::const_iterator MShapeItr;

/**
* A figure is just a group...
*/
class MFigure : public MNodeGroup {
private:
	// The index of the shape containing this figure
	MShape *m_shape;

	// Constructor
	MFigure(MGraph *graph=NULL,MShape *shape=NULL) : MNodeGroup(graph) {
		m_shape = shape;
	}

public:
	// Get the index of this figure in the parent
	int index();

	// Get the shape that contains this figure
	MShape *shape() {
		return m_shape;
	}

	// Get a handle on the medial axis for a link on this figure
	MAxelet &getMedialAxis(int link);

	// Compute the length of the figure
	vector<double> computeMedialParametrization();

	// Find a crossing of the figure with a line, nearest to
	// the specified point.  Return empty atom if crossing not found.
	MAtom findCrossing(const Vector2D &x,const Vector2D &n);

	friend class MNode;
	friend class MGraph;
};
typedef list<MFigure*>::const_iterator MFigureItr;

/**
* Graph of medial atoms and links
*
The topology of the graph can be modified by adding atoms and links to the graph.  Once all the
atoms and links have been added


*
* At
*/
class MGraph : public MNodeGroupIF {
private:
	// Medial Atoms 
	list <MNode *> m_nodes;

	// Links between different medial atoms
	list <MLink *> m_links;

	// Distinct figures contained in this graph
	list <MFigure *> m_figures;

	// Distinct shapes contained in this graph
	list <MShape *> m_shapes;

	// Current selection cache
	MNodeGroup m_selection;
	bool m_selectionDirty;

	// Flag indicating whether topology is being edited
	bool topologyEdit;

	void removeCurveLink(CurveLink *);
	void removeBranchLink(const BRef &bRef);
	void addCurveLink(MNode *head,MNode *tail);
	void clearMarks() const;	

	// Node list iterator for first/next
	MNodeItr firstNextItr;

public:
	virtual MNode *first() {
		firstNextItr = this->m_nodes.begin();
		return (firstNextItr == m_nodes.end()) ? NULL : *firstNextItr;
	}

	virtual MNode *next() {
		firstNextItr++;
		return (firstNextItr == m_nodes.end()) ? NULL : *firstNextItr;
	}

	virtual MGraph*graph() {
		return this;
	}

	// === Access to Nodes and Links =====================================
	const list<MNode*> &nodes() {
		return m_nodes;
	}

	const list<MLink*> &links() {
		return m_links;
	}

	const list<MFigure*> &figures() {
		return m_figures;
	}

	const list<MShape*> &shapes() {
		return m_shapes;
	}

	/**
	* Selection Retrieval
	*
	* The graph maintains a special group of nodes that are currently selected.
	* When nodes are selected or deselected, the selection is marked dirty and 
	* recomputed on next retrieval.  Selection is only maintained for the default 
	* selection mark, 0x0001.  In topology edit mode, selection is recomputed each
	* time it is asked for.
	*/
	MNodeGroup *selection() {
		if(m_selectionDirty || topologyEdit) {
			m_selection = MNodeGroupIF::ggSelected();
			if(!topologyEdit)
				m_selectionDirty = false;
		}
		return &m_selection;
	}

	void dirtySelection() {
		m_selectionDirty = true;
	}

	int nSelected() {
		return selection()->size();
	}

	MNodeGroup ggSelected() {
		return *selection();
	}

	// === Figure topology manipulation methods ==========================
	void beginTopologyEdit();
	void endTopologyEdit();

	// Add atom, place it into a figure
	MNode *insertAtomIntoFigure(MFigure *figure,int position,const MAtom &atom);

	// Create a new figure and add it to the object
	void addFigure(const vector<MAtom> &atoms,MNode **outNodeList = NULL);

	// Remove atom from a figure
	const MAtom &removeNode(MNode *node);

	// Remove figure
	void removeGroup(MNodeGroupIF &group);

	// This method adds a group of atoms (and connecting links) from one graph to another
	// graph.  The source group should be from a different graph that the target group!
	void addGroup(MNodeGroupIF &source);

	// Remove all figures
	void removeAll();

	// === Interfigural Relationships ====================================
	void addBranch(const BRef &tail,const BRef &head,bool flip);

	// Remove a branch at this boundary atom
	void removeBranch(const BRef &tail);

	// === IO ============================================================
	// void print(ostream &s) const;
	void loadFromModelFile(Registry &folder);
	void writeToModelFile(Registry &folder);

	// Keeping up with technology...
	void writeXML(const char *fname);

	// Writes the model to a  matlab file
	void writeToMatlabFile(char *fname);

	MGraph() {
		topologyEdit = false;
		m_selectionDirty = true;
	}

	MGraph(const MGraph &src) {
		topologyEdit = false;
		m_selectionDirty = true;
		*this = src;
	}

	virtual ~MGraph();

	// Assignment operator
	MGraph& operator =(const MGraph &graph);

	// Quickly append a node to one of the figures in the graph.  This add method is 
	// different from the other topology editors because it allows us to add a node 
	// without rebuilding the graph.
	// MNode *quickAppendNode(MFigure *figure,const MAtom &atom);


	// A static method to create a simple model
	static MGraph* createSaucage(int nAtoms,double r=0.1,double span=0.1);
};

class MException {
private:
	string m_text;
public:
	MException(char *text) : m_text(text) {      
	}

	MException() {
	}
	/*
	void print(ostream &s) const {
	s << m_text;
	}
	*/
};

/*

// Overloaded display operators.  
ostream& operator << (ostream &s,const MAtom &atom) {
atom.print(s);
return s;
}

ostream& operator << (ostream &s,const MGraph &graph) {
graph.print(s);
return s;
}

ostream& operator << (ostream &s,const MException &exception) {
exception.print(s);
return s;
}
*/

/******************
Some inlined methods that had to be defined after the delarations
****************************************************************/
inline MNode *MNode::getSuccessor() {
	return m_links[SUCC] ? m_links[SUCC]->head() : NULL;
}

inline MNode *MNode::getPredecessor() {
	return m_links[PRED] ? m_links[PRED]->tail() : NULL;
}

// Selection methods are overriden to notify the parent graph
inline void MNode::select()	{
	MAtom::select();
	m_graph->dirtySelection();
}

// Unselect a atom using a selection mask.  Whatever bits are set in the mask will be cleared in the atom.
inline void MNode::deselect() {
	MAtom::deselect();
	m_graph->dirtySelection();
}

// Toggle selection of the atom using a mask.  The bit positions where mask is 1 will be flipped
inline void MNode::toggle()	{
	MAtom::toggle();
	m_graph->dirtySelection();
}



#endif 
