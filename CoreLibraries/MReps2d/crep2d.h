#ifndef _CREP_H_
#define _CREP_H_

#include "phquintic.h"
#include <matrix/src/matrix.h>
#include <DSL2DLib/src/DSL2D.h>
/*
class CBoundlet {
private:
   // This describes the curve
   PHQuintic curve;
   
   // Vectors that describe normals of the ends of the boundlet when it
   // is mapped onto a unit interval.
   Vector2D vStart,vEnd;

   // A transform that maps the boundlet into real space
   Transform2D transform,inverseTransform;

   // This multiplyer is applied to the normals
   double nMult;

public:
   // Create a dummy boundlet
   CBoundlet();

   // Get the interpolation of the boundary.  The vectors that are input as references will
   // be returned as 3 element vectors with 1 in the third position for vector x and
   // 0 in the third position for vector n.
   void getInterpolation(double t,Vector2D &p,Vector2D &n);

   // Compute the boundlet given absolute position of the two boundary
   // sites.  The normal vector direction is the desired normal to the
   // boundary and the length is the desired radius of curvature
   void compute(const Vector2D &xStart,const Vector2D &xEnd,const Vector2D &nStart,const Vector2D &nEnd);
};
*/

// A set of boundary chunks for a given link
class DSLLinkBoundary2D {
private:
   
public:
   // Is this a dirty link
   bool dirty;

   // Tho boundlets
	CBoundlet bLeft,bRight;

   DSLLinkBoundary2D() {
      dirty = false;
   };

   DSLLinkBoundary2D(const DSLPrimitive2D &s,const DSLPrimitive2D &e) {
      setEnds(s,e);
   }

	void setEnds(const DSLPrimitive2D &s,const DSLPrimitive2D &e);
};

// A set of boundary chunks for a given link
class DSLEndBoundary2D {
private:
   
public:
   // Is this a dirty link
   bool dirty;

   // Tho boundlets
	CBoundlet bLeft,bRight;

   // Cached lengths of the boundary
   double lengthLeft,lengthRight;

   DSLEndBoundary2D() {
      dirty = false;
   };

   DSLEndBoundary2D(const DSLPrimitive2D &p,double direction) {
      setEnd(p,direction);
   }

	void setEnd(const DSLPrimitive2D &p,double direction);
};

class DSLBoundary2D : public DSLFigure2DData {
private:
	// A collection of 'links'
	vector <DSLLinkBoundary2D> links;

	// A pair of end interpolants
	DSLEndBoundary2D first,last;

	// A dirty bit for the whole structure
	bool dirty;

	// Rebuild the links in the structure to match figure
	void build();

	// Rebuild the links in the structure to match figure
	void buildLink(int i);
   void buildEnd(bool firstLink);

   // The figure
   DSLFigure2D *figure;

   // The total length of the boundary
   bool arcLengthDirty;
   double totalArcLength;

public:
   // Constructor
   DSLBoundary2D();
   ~DSLBoundary2D();

   // Set figure pointer
   void setFigure(DSLFigure2D *figure);

   // This method will be called whenever one of the primitives in the object gets updated.
   // The data object can then recompute (or flag as dirty) the attributes derived from 
   // primitives in the owner figure.
   void primitiveChanged(int i);

   // This method will be called whenever a primitive is inserted into a figure, or removed from
   // a figure, or a major change occurs to the whole figure.
   void figureChanged();

   // The data object must be able to make a copy of itself and return a pointer.  
   DSLFigure2DData *makeCopy(DSLFigure2D *newFigure);

   // Compute boundary position and normal at a certain point between two primitives
   void interpolateBoundary(int link,bool right,double t,Vector2D &p,Vector2D &n);

   // Compute boundary position and normal at point on the endcap
   void interpolateEndCap(bool first,bool right,double t,Vector2D &p,Vector2D &n);

   // Interpolate the whole boundary.  Note that the boundary is a set of pieces that are not nessesarily in order.
   // t runs from 0 to tmax();
   void interpolateBoundary(double t,Vector2D &p,Vector2D &n);

   // Max for t.  Equal to number of boundlets
   double tmax() {
      return figure->nprimitives()*2+2;
   }

   // Get the total arc-length
   double getTotalArcLength() {
      if(arcLengthDirty) {
         arcLengthDirty = false;
         totalArcLength = first.bLeft.getArcLength() + first.bRight.getArcLength() + 
            last.bLeft.getArcLength() + last.bRight.getArcLength();
         for(int i=0;i<links.size();i++) {
            totalArcLength += links[i].bLeft.getArcLength() + links[i].bRight.getArcLength();
         }         
      }
      return totalArcLength;
   }

   // Access a boundlet
   CBoundlet &getEndBoundlet(bool isFirst,bool isRight) {
      if(isFirst) {
         if(dirty)
            build();
         else if(first.dirty) 
            first.setEnd(figure->getPrimitive(0),-1);
         return isRight ? first.bRight : first.bLeft;
      }
      else {
         if(dirty)
            build();
         else if(last.dirty) 
            last.setEnd(figure->getPrimitive(figure->nprimitives()-1),1);
         return isRight ? last.bRight : last.bLeft;
      }
   }

   CBoundlet &getLinkBoundlet(int link,bool isRight) {
      if(dirty)
         build();
      else if(links[link].dirty) 
         links[link].setEnds(figure->getPrimitive(link),figure->getPrimitive(link+1));
      return isRight ? links[link].bRight : links[link].bLeft;
   }
};



#endif