#include "mreps2D.h"
#include "interp2D.h"
#include <optima.h>

MIProblem::MIProblem(MNode *start,MNode *end,double t,double help) {
	// Get the boundary and set start and end locations on the boundary
	b = start->bAtom(MAtom::LEFT).bnd;	

	sLeft = start->bAtom(MAtom::LEFT).tShape;
	sRight = end->bAtom(MAtom::RIGHT).tShape;
	bLenLeft = end->bAtom(MAtom::LEFT).tShape - sLeft;
	bLenRight = start->bAtom(MAtom::RIGHT).tShape - sRight;

	this->t = t;
   this->help = help;

	// Compute the parameters of the line on which the atom should lie
	rx = start->x() * (1-t) + end->x() * t;
	ru = (end->x() - start->x()).getNormal();
	ru.normalize();  
}

// This method returns the intersection of two rays if they intersect, false otherwise
inline bool rayIntersect(Vector2D x1,Vector2D n1,Vector2D x2,Vector2D n2,double &d1,double &d2) {
	// Make sure vectors are not parallel
	double den = n2.dotProduct(n1.getNormal());
	if(den == 0)	return false;

	Vector2D dxBar = (x2-x1).getNormal() / den;
	
	d1 = n2.dotProduct(dxBar);
	d2 = n1.dotProduct(dxBar);
	
	return true;
}

double MIProblem::evaluate(const Vector &x) {
   // If x is out of limits, return a huge value
	for(int i=0;i<2;i++) 
      if(x(i) < 0.0 || x(i) > 1.0)
         return 1000;

	// Get the positions and normals to the boundary
	BAtom L = b->interpolate(sLeft + x(0)*bLenLeft);
	BAtom R = b->interpolate(sRight + x(1)*bLenRight);

   // Tweak the normals a tidsy bit
	// double help = x(2);
	if(fabs(help) > 0.0) {
		L.n = Transform2D::rotation(help * (M_PI/360) * sin(2*M_PI*x(0))) * L.n;
		R.n = Transform2D::rotation(help * (M_PI/360) * cos(3*M_PI*x(0))) * R.n;
	}

   // K,L are signed distances from L,R to middle ray
	// M,N are signed distances from midpoint of the ray to intersection points with L,R
	double k,l,m,n;
	if(!rayIntersect(L.x,L.n,rx,ru,k,m) || !rayIntersect(R.x,R.n,rx,ru,l,n)) {
		return 1000;
	}

   // Initial penalty is the sum of distance along ray and difference in k,l
	double dist1 = fabs(m-n);
   double dist2 = fabs(fabs(k)-fabs(l));

   // Additional penalty is accrued if k+l > sqrt(2)*dist(pt1,pt2) (object angle under 45 deg)
   double d2 = R.x.distanceToSqr(L.x);
   double p1 = k*k / d2;
	double p2 = l*l / d2;
	double penalty = (p1<2 ? 0 : p1-2) + (p2<2 ? 0 : p2-2);

   // Done - return score
	double p = dist1 + dist2 + penalty;

   return dist1 + dist2 + penalty;
}

MAtom MIProblem::getAtom(const Vector &x) {
	// Target atom
	MAtom ma;

	// Get the positions and normals to the boundary
	BAtom L = b->interpolate(sLeft + x(0)*bLenLeft);
	BAtom R = b->interpolate(sRight + x(1)*bLenRight);

   // Tweak the normals a tidsy bit
	//double help = x(2);
	if(fabs(help) > 0.0) {
		L.n = Transform2D::rotation(help * (M_PI/360) * sin(2*M_PI*x(0))) * L.n;
		R.n = Transform2D::rotation(help * (M_PI/360) * cos(3*M_PI*x(0))) * R.n;
	}

   // K,L are signed distances from L,R to middle ray
	// M,N are signed distances from midpoint of the ray to intersection points with L,R
	double k,l,m,n;
	if(!rayIntersect(L.x,L.n,rx,ru,k,m) || !rayIntersect(R.x,R.n,rx,ru,l,n)) {
		return 1000;
	}

   // Take the average of m and n as the distance from rx,rt
   double dr = (m+n) / 2;

   // Coordinates of the primitive
   ma.x(rx + ru*dr);

   // Compare distance from both ends
   // double ds = sqrt((mi.start.x()-xi)*(mi.start.x()-xi)+(mi.start.y()-yi)*(mi.start.y()-yi));
   // double de = sqrt((mi.end.x()-xi)*(mi.end.x()-xi)+(mi.end.y()-yi)*(mi.end.y()-yi));

   // Average the radii
   ma.r((fabs(l) + fabs(k)) / 2.0);

   // Axial angle and object angle
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


// This method interpolates a selected figure
void interpolateFigure(MFigure *figure,int n) {
	// Compute the step required
	double tStep = 1.0 / n;

	// The list to which we will add atoms
	vector<MAtom> aList;

	// Go through all nodes
	for(int i=0;i<figure->size() - 1;i++) {
		// Push the existing node
		aList.push_back(*figure->node(i));

		// Get a handle of the axis
		MAxelet &axelet = figure->getMedialAxis(i);

		// A couple points
		MAtom x0 = *figure->node(i);
		MAtom x1 = *figure->node(i+1);
		Vector2D d = x1.x()-x0.x();
		Vector2D n = d.getNormal();

		// Compute each interpolant
		for(double t=tStep;t<1.0 - 0.5*tStep;t+=tStep) { 
			Vector2D x = x0.x() + d * t;
			MAtom a = axelet.interpolateAcross(x,n);
			aList.push_back(a);
		}

		if(i==figure->size()-2)
			aList.push_back(*figure->node(i+1));
	}

	// Add the new list of atoms as a figure
	MGraph *graph = figure->graph();
	graph->beginTopologyEdit();
	graph->removeGroup(*figure);
	graph->addFigure(aList);
	graph->endTopologyEdit();

	// Select the figure
	graph->deselect();
	graph->nodes().back()->figure()->select();
}

/*

// This method interpolates a selected figure
void interpolateFigure(MFigure *figure,int n) {
	// Compute the step required
	double tStep = 1.0 / n;

	// The list to which we will add atoms
	vector<MAtom> aList;

	// Create a graph to do this on.  We need this because the figure involved might be a part of a larger shape
	MGraph g;
	g.beginTopologyEdit();
	g.addFigure(figure->getAtomArray());
	g.endTopologyEdit();

	// Get the figure from the new graph
	MFigure *f = g.figures().front();

	// Go through all nodes
	for(int i=0;i<f->size()-1;i++) {
		MNode *s = f->node(i);
		MNode *e = f->node(i+1);

		aList.push_back(*s);
		for(int j=1;j<n;j++) {
			MAtom a = getMedialInterpolation(s,e,tStep*j);
			aList.push_back(a);
		}
		if(i==f->size()-2)
			aList.push_back(*e);
	}

	// Add the new list of atoms as a figure
	model.beginTopologyEdit();
	model.removeGroup(*figure);
	model.addFigure(aList);
	model.endTopologyEdit();

	// Select the figure
	model.deselect();
	model.nodes().back()->figure()->select();
}
*/

// This method will resample a figure at even intervals
void resampleFigureUniform(MFigure *figure) {
	// Compute the position of each node along the medial axis
	vector<double> tArray = figure->computeMedialParametrization();

	// The list to which we will add atoms
	vector<MAtom> aList;

	double tStep = tArray[figure->size()-1] / (figure->size()-1);
	int curAtom = 0;
	
	// Go through all nodes
	aList.push_back(*figure->node(0));	
	for(int i=1;i<figure->size()-1;i++) {
		double t = i*tStep;
		while(tArray[curAtom] < t)
			curAtom++;
		
		MNode *s = figure->node(curAtom-1);
		MNode *e = figure->node(curAtom);
		double x = (t-tArray[curAtom-1]) / (tArray[curAtom]-tArray[curAtom-1]);

		Vector2D lineX = (e->x() - s->x()) * x + s->x();
		Vector2D lineN = (e->x() - s->x()).getNormal();
		MAxelet &ax = figure->getMedialAxis(curAtom-1);
		MAtom a = ax.interpolateAcross(lineX,lineN);

		a.select();
		aList.push_back(a);
	}

	aList.push_back(*figure->node(figure->size()-1));

	// Add the new list of atoms as a figure
	MGraph *graph = figure->graph();
	graph->beginTopologyEdit();
	graph->removeGroup(*figure);
	graph->addFigure(aList);
	graph->endTopologyEdit();
}


