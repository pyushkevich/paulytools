#include "crep2d.h"

/*
CBoundlet::CBoundlet() : vStart(0,0,true), vEnd(0,0,true) {
}

void CBoundlet::compute(const Vector2D &xStart,const Vector2D &xEnd,
                        const Vector2D &nStart,const Vector2D &nEnd) 
{
   

   // We will compute the transform that takes a vector in real space and maps it
   // onto the normalized space.  This vector maps (xEnd-xStart) to (1,0)   
   transform = Transform2D::frameChange(xStart,xEnd-xStart);
   inverseTransform = transform.inverse();

   // Test the inverse
   Transform2D T = transform * inverseTransform;
      
   // Now let's map all the normal vectors to our unit space
   Vector2D vs = inverseTransform * nStart;
   Vector2D ve = inverseTransform * nEnd;

   // Now let's compare the old vectors to the new vectors.  If the difference in one-norm
   // is less than epsilon, the quintic stands
   double dv = (vStart-vs).oneNorm() + (vEnd-ve).oneNorm();

   // OK, store the new vectors and go home
   vStart = vs;
   vEnd = ve;
   
   if(dv > 0.00001) {
      // Now, remember that the user passes in normal vectors, but what we need is tangent vectors,
      // and in the right direction.
      
      vs = Vector2D(vs.y,-vs.x,true);
      ve = Vector2D(ve.y,-ve.x,true);

      if(vs.x < 0) {
         nMult = -1;
         vs = -vs;
         ve = -ve;
      }
      else {
         nMult = 1;
      }

      curve.setEndVectors(vs.x,vs.y,ve.x,ve.y);
   }
}

// Get the interpolation of the boundary
void CBoundlet::getInterpolation(double t,Vector2D &p,Vector2D &n) {
   Vector2D P(0,0,true),N(0,0,false);
   
   // Get interpolation from ph curve
   curve.getInterpolation(t,P.x,P.y,N.x,N.y);

   // Lets transform these vectors into the real world
   p = transform * P;
   n = (transform * N)*nMult;


}
*/

/***********************************************************************************
 * This method returns a list of boundary primitives for a figure
 ***********************************************************************************/
void DSLLinkBoundary2D::setEnds(const DSLPrimitive2D &s,const DSLPrimitive2D &e) {
   for(int which=1;which<=2;which++) {
      // Positions of the end primitives
      Vector2D xs = s.getPosition(which);
      Vector2D xe = e.getPosition(which);

      // A cheat that makes the boundary look nice and smooth.  The length of the
      // tangent vector is proportional to the length of the segment
      double d = 0.8 * xs.distanceTo(xe);   
      Vector2D ns = s.getNormal(which) * d;
      Vector2D ne = e.getNormal(which) * d;

      // Compute left or right boundary
      if(which==1)
         bLeft.compute(xs,xe,ns,ne);
      else
         // Right side is flipped, t going backward
         bRight.compute(xs,xe,ns,ne);
   }
 
   // A cheat that makes the boundary look nice and smooth.  The length of the
   // tangent vector is proportional to the length of the segment
   dirty = false;
}

void DSLEndBoundary2D::setEnd(const DSLPrimitive2D &p,double direction) {
   Vector2D xs,ns;

   // Distance to the endcap
   double de = p.getEndRC() * p.getR() * direction;
  
   // Compute position and normal at end cap
   xs = p.getPosition() + p.getNormal(0) * de;

   // Length of the normal vector depends on the distance between curve ends
   double nLen = 0.8 * p.getPosition(1).distanceTo(xs);
   ns = p.getNormal(0) * nLen;
   
   // bLeft.compute(xs,p.getPosition(1),ns,p.getNormal(1) * nLen);
   bLeft.compute(p.getPosition(1),xs,p.getNormal(1) * nLen,ns * direction);
   bRight.compute(p.getPosition(2),xs,p.getNormal(2) * nLen,ns * direction);

   dirty = false;
}

DSLBoundary2D::DSLBoundary2D() {
   dirty = false;
   arcLengthDirty = false;
   figure = NULL;

   //printf("DSLBoundary2D constructed\n");
}

DSLBoundary2D::~DSLBoundary2D() {
   // printf("DSLBoundary2D destroyed\n");
}

void DSLBoundary2D::setFigure(DSLFigure2D *figure) {
   this->figure = figure;
   dirty = true;
   arcLengthDirty = true;
}

void DSLBoundary2D::figureChanged() {
   dirty = true;
   arcLengthDirty = true;
}

void DSLBoundary2D::primitiveChanged(int p) {
   if(dirty)
      return;

   if(p==0) 
      first.dirty = true;
   
   if(p==figure->nprimitives()-1) 
      last.dirty = true;
   
   if(p > 0) 
      links[p-1].dirty = true;

   if(p < links.size())
      links[p].dirty = true;

   arcLengthDirty = true;
}

DSLFigure2DData *DSLBoundary2D::makeCopy(DSLFigure2D *newFigure) {
   DSLBoundary2D *rtn = new DSLBoundary2D();
   *rtn = *this;
   rtn->figure = newFigure;
   return rtn;
}


void DSLBoundary2D::build() {
   dassert(figure);
   
   first.setEnd(figure->getPrimitive(0),-1);
   last.setEnd(figure->getPrimitive(figure->nprimitives()-1),1);

   links.clear();
   for(int i=1;i<figure->nprimitives();i++) {
      links.push_back(DSLLinkBoundary2D(figure->getPrimitive(i-1),figure->getPrimitive(i)));
   }

   dirty = false;
}

// Compute boundary position and normal at a certain point between two primitives
void DSLBoundary2D::interpolateBoundary(int link,bool right,double t,Vector2D &p,Vector2D &n) {
   dassert(figure);

   if(dirty)
      build();
   else if(links[link].dirty) {
      links[link].setEnds(figure->getPrimitive(link),figure->getPrimitive(link+1));
   }

   if(right)
      links[link].bRight.getInterpolation(t,p,n);
   else
      links[link].bLeft.getInterpolation(t,p,n);

   n = n.getNormal();
   n.normalize();
}

// Compute boundary position and normal at point on the endcap
void DSLBoundary2D::interpolateEndCap(bool firstLink,bool right,double t,Vector2D &p,Vector2D &n) {
   dassert(figure);

   // Pick an end
   DSLEndBoundary2D &end = (firstLink) ? first : last;
   int primNo = (firstLink) ? 0 : figure->nprimitives()-1;

   if(dirty)
      build();
   else if(end.dirty) {
      end.setEnd(figure->getPrimitive(primNo),(firstLink) ? -1 : 1);
   }
   
   if(right)
      end.bRight.getInterpolation(t,p,n);
   else
      end.bLeft.getInterpolation(t,p,n);

   n = n.getNormal();
   n.normalize();
}

// Interpolate the whole boundary.  t ranges from 0 to tmax
void DSLBoundary2D::interpolateBoundary(double t,Vector2D &p,Vector2D &n) {
   int segment = fmod(t,tmax());
   double tseg = fmod(t,1.0);

   bool right = (segment > figure->nprimitives());
   int link = (right) ? 2*figure->nprimitives()+1-segment : segment;

   if(link == 0 || link==figure->nprimitives()) 
      interpolateEndCap(link==0,right,tseg,p,n);
   else
      interpolateBoundary(link-1,right,tseg,p,n);
}













