#include "optimization.h"
#include <xmmintrin.h>
#include <iostream>
#include <string>
#include <vnl/vnl_rotation_matrix.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

inline float SplinePrior::priorProbabilityFull(SplineDataCache *sdc) {
  return priorProbability(sdc,0,0,spline->dim(0),spline->dim(1)); 
}

float SplinePrior::priorProbability(SplineDataCache *sdc,int cuFirst,int cvFirst,int cuLast,int cvLast) {
  return 0;
}

// Prior probability for a region of control points
float SplinePatchPrior::priorProbabilityFull(SplineDataCache *sdc) {
  float ppTotal = 0;
  float wgtTotal = 0;
  float wgtDefault = 1.0f / paTimeStamp.size();
  float wgt;

  for (int iPatch=0;iPatch<sdc->patch.width();iPatch++)
    {
    for (int jPatch=0;jPatch<sdc->patch.height();jPatch++)
      {
      if (spline->getPatchTS(iPatch,jPatch) > paTimeStamp(iPatch,jPatch))
        {
        paTimeStamp(iPatch,jPatch) = spline->getPatchTS(iPatch,jPatch);

        wgt = wgtDefault;
        paMeasure(iPatch,jPatch) = priorOverPatch(sdc,iPatch,jPatch,wgt);
        paWeight(iPatch,jPatch) = wgt;
        }
      ppTotal += paMeasure(iPatch,jPatch);
      wgtTotal += paWeight(iPatch,jPatch);
      }
    }

  return ppTotal / wgtTotal;
}

// Constructor
SplinePatchPrior::SplinePatchPrior(MSpline *spline)
: SplinePrior(spline) 
{
  int w = spline->dim(0)-2;
  int h = spline->dim(1)-2;
  paTimeStamp.resize(w,h);
  paWeight.resize(w,h);
  paMeasure.resize(w,h);
  paTimeStamp.setAll(-2);
}

float distanceVariance(MSpline *spline,int cuFirst,int cvFirst,int cuLast,int cvLast) {
  // Compute the range of the application of the prior
  int cu0 = cuFirst == 0 ? 0 : cuFirst - 1;
  int cv0 = cvFirst == 0 ? 0 : cvFirst - 1;
  int cu1 = cuLast == spline->dim(0) ? cuLast : cuLast+1;
  int cv1 = cvLast == spline->dim(1) ? cvLast : cvLast+1;
  int u,v,i;

  // Prior probability
  float prior = 0;

  // Create a bag of distances between neighboring atoms
  vector<float> dvec;

  // Compute the horizontal and vertical
  for (v=cv0;v<=cv1;v++)
    {
    for (u=cu0;u<cu1;u++)
      {
      SMLVec3f delta = spline->getControl(u,v) - spline->getControl(u+1,v);
      dvec.push_back(delta.magnitude());
      }
    }
  for (v=cv0;v<cv1;v++)
    {
    for (u=cu0;u<=cu1;u++)
      {
      SMLVec3f delta = spline->getControl(u,v) - spline->getControl(u,v+1);
      dvec.push_back(delta.magnitude());
      }
    }

  // Compute the mean of the links
  float mean = 0;
  for (i=0;i<dvec.size();i++)
    {
    mean += dvec[i];
    }
  mean /= dvec.size();

  // Compute the variance
  float var = 0;
  for (i=0;i<dvec.size();i++)
    {
    var += (dvec[i] - mean)*(dvec[i] - mean);
    }
  var /= (dvec.size()-1);

  // Return the variance
  return sqrt(var);
}


float ControlPointConstraint::priorProbability(SplineDataCache *sdc,int cuFirst,int cvFirst,int cuLast,int cvLast) {
  // The penalty term
  float penalty = 0;

  // Make sure that all control points are inside the control box
  for (int j=cvFirst;j<=cvLast;j++)
    {
    for (int i=cuFirst;i<=cuLast;i++)
      {
      SMLVec3f X = spline->getControl(i,j);
      for (int d=0;d<3;d++)
        {
        float x = X[d];
        if (x > 1.0f)
          {
          penalty += X[0];
          }
        else if (x < 0.0f)
          {
          penalty += 1.0f - x;
          }
        }
      }
    }

  return penalty;
}




float smoothness(MSpline *spline,int cuFirst,int cvFirst,int cuLast,int cvLast) {

  // Compute the range of the application of the prior
  int cu0 = cuFirst == 0 ? 0 : cuFirst - 1;
  int cv0 = cvFirst == 0 ? 0 : cvFirst - 1;
  int cu1 = cuLast == spline->dim(0) ? cuLast : cuLast+1;
  int cv1 = cvLast == spline->dim(1) ? cvLast : cvLast+1;
  int u,v;

  // Prior probability
  float prior = 0;

  // Compute the horizontal links
  for (v=cv0;v<=cv1;v++)
    {
    // Displacement at the left of the link
    SMLVec4f dfLeft = spline->getControl4f(cu0,v);

    for (u=cu0;u<cu1;u++)
      {
      // Displacement at the right of the link
      SMLVec4f dfRight = spline->getControl4f(u+1,v);

      // Compute the prior
      prior += (dfRight-dfLeft).squared_magnitude();

      // Update the left link
      dfLeft = dfRight;
      }
    }

  // Compute the vertical links
  for (u=cu0;u<=cu1;u++)
    {
    // Displacement at the left of the link
    SMLVec4f dfTop = spline->getControl4f(u,cv0);

    for (v=cv0;v<cv1;v++)
      {
      // Displacement at the right of the link
      SMLVec4f dfBottom = spline->getControl4f(u,v+1);

      // Compute the prior
      prior += (dfBottom-dfTop).squared_magnitude();

      // Update the left link
      dfTop = dfBottom;
      }
    }

  // Return the prior information
  return prior;
}

float MarkovControlPrior::markovPairPenalty(int i1,int j1,int i2,int j2) {
  // Distance minimization
  SMLVec3f p1 = spline->getControl(i1,j1);
  SMLVec3f p2 = spline->getControl(i2,j2);
  return (p1-p2).squared_magnitude();




  /*
  SMLVec3f p1 = spline->getControl(i1,j1);
  SMLVec3f p2 = spline->getControl(i2,j2);
  SMLVec3f o1(ctlBase(i1,j1).x,ctlBase(i1,j1).y,ctlBase(i1,j1).z);
  SMLVec3f o2(ctlBase(i2,j2).x,ctlBase(i2,j2).y,ctlBase(i2,j2).z);

  return (p2+o1-(p1+o2)).LengthSquared() / (o1-o2).LengthSquared();
  */
  //(p1-o1).DistanceSquared(p2-o2) / o1.DistanceSquared(o2);
  //return p1.DistanceSquared(p2);

}

float MarkovControlPrior::priorProbabilityFull(SplineDataCache *sdc) {
  int m = spline->dim(0);
  int n = spline->dim(1);

  float prior = 0;

  // Loop over all control points in the lattice
  for (int ic=0;ic<=m;ic++)
    {
    for (int jc=0;jc<=n;jc++)
      {
      if (ic > 0)
        prior += markovPairPenalty(ic-1,jc,ic,jc);
      if (jc > 0)
        prior += markovPairPenalty(ic,jc-1,ic,jc);
      if (ic < m)
        prior += markovPairPenalty(ic,jc,ic+1,jc);
      if (jc < n)
        prior += markovPairPenalty(ic,jc,ic,jc+1);
      }
    }

  // Return the prior
  return prior;
}

float MarkovCurvaturePrior::priorProbabilityFull(SplineDataCache *sdc) {
  int m = spline->dim(0);
  int n = spline->dim(1);
  int ic,jc;

  float prior = 0;

  for (ic=0;ic<=m;ic++)
    {
    for (jc=0;jc<=n;jc++)
      {
      if (ic > 0 && ic < m)
        {
        const SMLVec3f &X0 = spline->getControl(ic-1,jc);
        const SMLVec3f &X1 = spline->getControl(ic,jc);
        const SMLVec3f &X2 = spline->getControl(ic+1,jc);

        SMLVec3f Xii = (X2+X0) - (X1+X1);
        prior += Xii.squared_magnitude();
        }
      if (jc > 0 && jc < n)
        {
        const SMLVec3f &X0 = spline->getControl(ic,jc-1);
        const SMLVec3f &X1 = spline->getControl(ic,jc);
        const SMLVec3f &X2 = spline->getControl(ic,jc+1);

        SMLVec3f Xjj = (X2+X0) - (X1+X1);
        prior += Xjj.squared_magnitude();
        }
      if (ic > 0 && jc > 0)
        {
        const SMLVec3f &X00 = spline->getControl(ic,jc);
        const SMLVec3f &X01 = spline->getControl(ic,jc-1);
        const SMLVec3f &X10 = spline->getControl(ic-1,jc);
        const SMLVec3f &X11 = spline->getControl(ic-1,jc-1);

        SMLVec3f Xij = (X00+X11) - (X01+X10);
        prior += Xij.squared_magnitude();
        }
      if (ic < m && jc < n)
        {
        const SMLVec3f &X00 = spline->getControl(ic,jc);
        const SMLVec3f &X01 = spline->getControl(ic,jc+1);
        const SMLVec3f &X10 = spline->getControl(ic+1,jc);
        const SMLVec3f &X11 = spline->getControl(ic+1,jc+1);

        SMLVec3f Xij = (X00+X11) - (X01+X10);
        prior += Xij.squared_magnitude();
        }
      }
    }

/*
    for(ic=0;ic<=m;ic++) {
        for(jc=1;jc<n;jc++) {
            SMLVec3f X01 = spline->getControl(ic-1,jc);
            SMLVec3f X10 = spline->getControl(ic,jc-1);
            SMLVec3f X11 = spline->getControl(ic,jc);
            SMLVec3f X12 = spline->getControl(ic,jc+1);
            SMLVec3f X21 = spline->getControl(ic+1,jc);
        
            SMLVec3f D01 = X11 - X01;
            SMLVec3f D10 = X11 - X10;           
            SMLVec3f D12 = X12 - X11;
            SMLVec3f D21 = X21 - X11;

            float d01 = D01.LengthSquared();
            float d10 = D10.LengthSquared();
            float d21 = D21.LengthSquared();
            float d12 = D12.LengthSquared();
            
            float alpha = D21.Dot(D01);
            float  beta = D12.Dot(D10);
            float gamma = D21.Dot(D12);
            float delta = D10.Dot(D01);

            prior += 2.0f;
            prior -= alpha*alpha / (d21*d01);
            prior -= beta * beta / (d12*d10);
            prior += gamma*gamma / (d21*d12);
            prior += delta*delta / (d01*d10);
        }
    }
*/
/*
    // Loop over all control points in the lattice
    for(ic=0;ic<=m;ic++) {
        for(jc=1;jc<n;jc++) {
            SMLVec3f X0 = spline->getControl(ic,jc-1);
            SMLVec3f X1 = spline->getControl(ic,jc);
            SMLVec3f X2 = spline->getControl(ic,jc+1);

            SMLVec3f Xjj = (X2+X0) - (X1+X1);
            prior += Xjj.LengthSquared();
        }
    }

    // Loop over all control points in the lattice
    for(jc=0;jc<=n;jc++) {
        for(ic=1;ic<m;ic++) {
            SMLVec3f X0 = spline->getControl(ic-1,jc);
            SMLVec3f X1 = spline->getControl(ic,jc);
            SMLVec3f X2 = spline->getControl(ic+1,jc);

            SMLVec3f Xii = (X2+X0) - (X1+X1);
            prior += Xii.LengthSquared();
        }
    }
*/
  // Return the prior
  return prior;
}




MarkovControlPrior::MarkovControlPrior(MSpline *spline) 
: SplinePrior(spline)
{
  ctlBase.resize(spline->dim(0)+1,spline->dim(1)+1);
  spline->getControlArray(ctlBase,0);
}


/**
 * Compute the prior probability on a range of values
 */
float CurvaturePrior::priorOverPatch(SplineDataCache *sdc,int iPatch,int jPatch,float &weight) {
  return sdc->patch(iPatch,jPatch)->integrateMedialMeasure(*this,weight);
}

float CurvaturePrior::computeMedialMeasure(const MedialPoint &mp) {
  return(mp.sinTheta2 < 0) ? sqrtf((mp.kappa1*mp.kappa1) + (mp.kappa2*mp.kappa2)) : 0;
}

MRepConstraint::MRepConstraint(Registry *settings, MSpline *spline) : SplinePatchPrior(spline) {
  epsRadius = settings->getDoubleValue("epsilonRadius",0.0001);
  epsBlum = settings->getDoubleValue("epsilonBlum",0.0001);
}

MRepConstraint::MRepConstraint(MSpline *spline) : SplinePatchPrior(spline) {
  epsRadius = 0.0001;
  epsBlum = 0.0001;
}

float MRepConstraint::priorOverPatch(SplineDataCache *sdc,int iPatch,int jPatch,float &weight) {
  return sdc->patch(iPatch,jPatch)->integrateBoundaryMeasure(*this,weight);
}

// Computes value > 0 ? 0 : value;


#ifdef CMREPS_SIMD
inline float maskNegative(float value) {
  __m128 r0,r1;

  r0 = _mm_set_ps1(0.0f);
  r1 = _mm_set_ps1(value);

  r0 = _mm_cmpgt_ss(r0,r1);
  r0 = _mm_and_ps(r0,r1);

  return *((float *)(&r0));
}
#else
inline float maskNegative(float value) {
  return value > 0 ? 0 : value;
}
#endif

// Check if any of the normals when dotted with the area vector for a triangle give a negative
// area 
inline float triangleTest(float *A,float *N1,float *N2,float *N3) {
  register __m128 r0,r1,r2,r3,r4;

  // Load the data
  //r0 = _mm_loadu_ps(A);
  //r1 = _mm_loadu_ps(N1);
  //r2 = _mm_loadu_ps(N2);
  //r3 = _mm_loadu_ps(N3);
  r0 = _mm_set_ps(0,A[2],A[1],A[0]);
  r1 = _mm_set_ps(0,N1[2],N1[1],N1[0]);
  r2 = _mm_set_ps(0,N2[2],N2[1],N2[0]);
  r3 = _mm_set_ps(0,N3[2],N3[1],N3[0]);

  // Multiply
  r1 = _mm_mul_ps(r1,r0);                     //  x       13      12      11
  r2 = _mm_mul_ps(r2,r0);                     //  x       23      22      21
  r3 = _mm_mul_ps(r3,r0);                     //  x       33      32      31

  // Should be using a 3x3 transpose
  r4 = _mm_unpacklo_ps(r1,r2);                //  22      12      21      11
  r0 = _mm_unpackhi_ps(r1,r2);                //  x       x       23      13
  r1 = _mm_shuffle_ps(r4,r3,0xC4);            //  x       31      21      11
  r2 = _mm_shuffle_ps(r4,r3,0xDE);            //  x       32      22      12 
  r3 = _mm_shuffle_ps(r0,r3,0xE4);            //  x       33      23      13

  // _MM_TRANSPOSE4_PS(r1,r2,r3,r0);

  // Add the elements
  r1 = _mm_add_ps(r1,r2);

  // Sneak in a set to zero instruction
  r0 = _mm_xor_ps(r0,r0);

  // Add the elements
  r1 = _mm_add_ps(r1,r3);

  // Compare the elements to 0
  r0 = _mm_cmpgt_ps(r0,r1);
  r0 = _mm_and_ps(r0,r1);

  // Add the elements of the vector
  r1 = _mm_shuffle_ps(r0,r0,0x01);
  r2 = _mm_shuffle_ps(r0,r0,0x02);
  r1 = _mm_add_ps(r1,r0);
  r1 = _mm_add_ps(r1,r2);

  // Return the last element in r1
  return *((float *)(&r1));
}


float MRepFastConstraint::priorOverPatch(SplineDataCache *sdc,int iPatch,int jPatch,float &weight) 
{
  // The error term
  float penalty = 0;

  // Go through all patches
  PatchDataCache *pdc = sdc->patch(iPatch,jPatch);

  // Count all the bad triangles
  for (int t=0;t<pdc->nTriangles;t++)
    {
    PatchTriangle &T = pdc->PT[t];
    if (T.in)
      {
      MedialPoint& MP0 = pdc->MP(T.v[0]);
      MedialPoint& MP1 = pdc->MP(T.v[1]);
      MedialPoint& MP2 = pdc->MP(T.v[2]);

      SMLVec3f NN = -T.NBnd[0];
      float p11 = triangleTest(
        NN.data_block(),
        MP0.bp[0].N.data_block(),MP1.bp[0].N.data_block(),MP2.bp[0].N.data_block());
      float p12 = triangleTest(
        T.NBnd[1].data_block(),
        MP0.bp[1].N.data_block(),MP1.bp[1].N.data_block(),MP2.bp[1].N.data_block());

      penalty += p11 + p12;
      }
    }

  return -10000.0f * penalty;
}

float MRepRadiusConstraint::priorOverPatch(SplineDataCache *sdc,int iPatch,int jPatch,float &weight) 
{
  // The error term
  float penalty = 0;

  // Go through all patches
  PatchDataCache *pdc = sdc->patch(iPatch,jPatch);

  // Count all the bad triangles
  for (int t=0;t<pdc->idxTrimmed.size();t++)
    {
    int idx = pdc->idxTrimmed[t];
    if(pdc->MP(idx).R() < 0) 
      cout << "Rad: " << pdc->MP(idx).R() << endl;
    penalty += maskNegative(pdc->MP(idx).R());
    }

  return -10000.0f * penalty;
}

float applyPenalty(float value, float eps) {


  //return (value < 0) ? 1.0f+value*value : 0;

  if (value >= eps)
    {
    return 0;
    }
  else if (value >= 0)
    {
    return 1000 * (1.0 - value / eps);
    //return (1.0 - value / eps) * (1.0 - value / eps);
    }
  else
    {
    return 1000 - value;
    }

  // This is bound to be a very large value...
  // return (1 - value / eps)*(1 - value / eps); */
}

float MRepConstraint::computeBoundaryMeasure(const MedialPoint &mp,int side) {
  // Compute the applicable penalties at a point
  float penalty = 0;

  // Negative radius is not allowed
  penalty += applyPenalty(mp.R(),epsRadius);

  // Minimum curvature may not be less than -1/r
  penalty += applyPenalty(1.0f + mp.bp[side].kappa2 * mp.R(),-epsBlum);   

  // Return the penalty
  return penalty;
}

float MRepConstraint::computeCrestBoundaryMeasure(const MedialPoint &mp) {
  // Compute the applicable penalties at a point
  float penalty = 0;

  // Negative radius is not allowed
  penalty += applyPenalty(mp.R(),epsRadius);

  // Return the penalty
  return penalty;
}






void SplineOptProblem::makeSDVector(Vector &vec) {

  // Compute the length of the control point spread
  double sx = 0.05;
  double sy = 0.05;
  double sz = 0.05;
  double sr = 0.05;

  int off = 0;
  for (int v=cpFirstV;v<=cpLastV;v++)
    {
    for (int u=cpFirstU;u<=cpLastU;u++)
      {
      vec(off++) = sx;
      vec(off++) = sy;
      vec(off++) = sz;
      vec(off++) = sr;
      }
    }   
}

inline int SplineOptProblem::getVectorSize() {
  return vecSize;
}

void SplineOptProblem::makeVector(Vector &vec) {
  int off = 0;
  for (int v=cpFirstV;v<=cpLastV;v++)
    {
    for (int u=cpFirstU;u<=cpLastU;u++)
      {
      vec(off++) = spline->getControl(u,v,0);
      vec(off++) = spline->getControl(u,v,1);
      vec(off++) = spline->getControl(u,v,2);
      vec(off++) = spline->getControl(u,v,3);
      }
    }
}

void SplineOptProblem::applyVector(const Vector &vec) {
  int off = 0;
  for (int v=cpFirstV;v<=cpLastV;v++)
    {
    for (int u=cpFirstU;u<=cpLastU;u++)
      {
      spline->setControl(u,v,0,vec(off++));
      spline->setControl(u,v,1,vec(off++));
      spline->setControl(u,v,2,vec(off++));
      spline->setControl(u,v,3,vec(off++));
      }
    }
}

double SplineOptProblem::evaluate(const Vector &vec) {
  // Apply the vector
  applyVector(vec);

  // Compute the image match over the affected patches in a nice timestamped fashion
  float lkhd = match->getMatch();

  // Apply prior terms to the deformation
  // float pp0 = constrControl.priorProbability(spline,splineData,cpFirstU,cpFirstV,cpLastU,cpLastV);
  // cout << "pp0 = " << pp0 << endl;
  /*
  if(pp0 > 0) {
      applyVector(backup);
      return 1000.0f * pp0;
  }*/

  // float pp = priorMCP.priorProbability(spline,splineData,cpFirstU,cpFirstV,cpLastU,cpLastV);
  // float pp1 = priorCrv.priorProbabilityFull(spline,splineData);
  //float pp2 = constrMRep.priorProbability(spline,splineData,cpFirstU,cpFirstV,cpLastU,cpLastV);

  // float pp1 = priorCrv.priorProbabilityFull(spline,splineData);
  // float pp2 = constrMRep.priorProbabilityFull(spline,splineData);

  // cout << "lkhd = " << lkhd << "\t prior1 = " << pp1 << "\t prior2 = " << pp2 << endl;

  //float pp1 = priorFast.priorProbability(spline,splineData,cpFirstU,cpFirstV,cpLastU,cpLastV);
  //float pp2 = priorRad.priorProbability(spline,splineData,cpFirstU,cpFirstV,cpLastU,cpLastV);

  float pp1 = priorFast.priorProbabilityFull(splineData);
  float pp2 = priorRad.priorProbabilityFull(splineData);
  // float pp3 = priorCtl.priorProbabilityFull(splineData);
  float pp4 = priorCtl.priorProbabilityFull(splineData);

  // Unapply the vector
  applyVector(backup);

  // If the count is all zero, return (?) zero
  // float total = (lkhd + wgtPriorCrv*pp1 + wgtConstrMRep*pp2);
  // float total = (lkhd + wgtConstrMRep*pp2 + wgtPriorCrv*pp1);
  float total = (lkhd + wgtFast*pp1 + wgtRadius*pp2 + wgtPriorCtl*pp4);
  //float total = lkhd;

  // printf("%10g\t%10g\t%10g\t%10g\n",lkhd,pp1,pp2,total);

  // return total * 0.001;
  evaluationCost++;

  return total;
}

SplineOptProblem::SplineOptProblem(MSpline *spline, SplineDataCache *cache, SplineImageMatcher *match, 
                                   int cpFirstU, int cpFirstV, int cpRangeU, int cpRangeV,
                                   Registry *stg)
: priorFast(spline), priorRad(spline), priorCtl(spline), priorCrv(spline)
{
  this->spline = spline;
  this->splineData = cache;
  this->match = match;

  this->cpFirstU = cpFirstU;
  this->cpFirstV = cpFirstV;
  this->cpRangeU = cpRangeU;
  this->cpRangeV = cpRangeV;
  this->cpLastU = cpFirstU + cpRangeU - 1;
  this->cpLastV = cpFirstV + cpRangeV - 1;

  // Size of the vector
  vecSize = cpRangeU*cpRangeV*4;

  // Create a backup vector for optimization
  backup.setSize(vecSize);
  makeVector(backup);

  // Figure out the patches affected by this optimization
  uPatchFirst = cpFirstU < 4 ? 0 : cpFirstU - 3;
  vPatchFirst = cpFirstV < 4 ? 0 : cpFirstV - 3;
  uPatchLast = cpLastU > spline->dim(0)-4 ? spline->dim(0) - 3 : cpLastU;
  vPatchLast = cpLastV > spline->dim(1)-4 ? spline->dim(1) - 3 : cpLastV;

  // Configure prior weights
  // wgtPriorCrv = 0.1f;
  wgtPriorCtl = (float) stg->getDoubleValue("wgt.cpPrior",0.2);
  wgtConstrMRep = (float) stg->getDoubleValue("wgt.mrepConstraint",50.0);
  wgtFast = (float) stg->getDoubleValue("wgt.selfCrossConstraint",1.0);
  wgtRadius = (float) stg->getDoubleValue("wgt.radiusConstraint",1.0);;
  wgtPriorCrv = (float) stg->getDoubleValue("wgt.crvPrior",0.2);;
}


// Constructor.  Takes a range of control points to optimize over.
SplineRigidProblem::SplineRigidProblem(MSpline *spline, SplineDataCache *cache, SplineImageMatcher *match) {
  // Save the parameters
  this->spline = spline;
  this->splineData = cache;
  this->match = match;

  // Create the backup vector
  B.resize(spline->dim(0)+1,spline->dim(1)+1);

  // Compute the spline's center point    
  C.fill(0);
  for (int i=0;i<=spline->dim(0);i++)
    {
    for (int j=0;j<=spline->dim(1);j++)
      {
      B(i,j) = spline->getControl4f(i,j);
      C += SMLVec3f(B(i,j)[0],B(i,j)[1],B(i,j)[2]);
      }
    }
  C /= (spline->dim(0)+1) * (spline->dim(1)+1);
}

// Evaluate a solution
double SplineRigidProblem::evaluate(const Vector &v) {
  // Apply the vector
  applyVector(v);

  // Compute the image match
  float lkhd = match->getMatch();
  float totalArea = 0;

  // Compute the area of the spline
  for (int ip=0;ip<splineData->patch.size();ip++)
    {
    PatchDataCache *pdc = splineData->patch.getData()[ip];
    totalArea -= pdc->awBndTotal;
    }

  // Unapply the vector
  revertSpline();

  // Compute the prior term
  // float pp = -log(penaltyFunction(v(3),0,0.2));

  // Print out pp
  // cout << lkhd << "\t" << pp << endl;

  // The result is the ratio of match to area
  // return lkhd + pp;
  cout << lkhd << "\t" << totalArea << endl;
  return lkhd + totalArea * 10;
}

// Apply a vector to the spline
void SplineRigidProblem::applyVector(const Vector &v) {

  float scale = pow(2.0,v(3));

  // Create an appropriate matrix
  SMLMatrix3f S,M1,M2,M;
  SMLMatrix3f Rx = vnl_matops::d2f(vnl_rotation_matrix(SMLVec3d(v(0),0,0)));
  SMLMatrix3f Ry = vnl_matops::d2f(vnl_rotation_matrix(SMLVec3d(0,v(1),0)));
  SMLMatrix3f Rz = vnl_matops::d2f(vnl_rotation_matrix(SMLVec3d(0,0,v(2))));

  S.fill(0);
  S.fill_diagonal(scale);

  // Compute the transformation matrix
  M = (Rx * Ry) * (Rz * S);

  // The translation
  SMLVec3f T(v(4),v(5),v(6));

  // Transform the points
  for (int i=0;i<=spline->dim(0);i++)
    {
    for (int j=0;j<=spline->dim(1);j++)
      {
      // Get the new X coordinate of the spline control point
      SMLVec3f X(B(i,j)[0],B(i,j)[1],B(i,j)[2]);
      SMLVec3f Y = M * (X-C);
      spline->setControl(i,j,Y + C + T);

      // Scale the radius of the spline control point
      spline->setControl(i,j,3,B(i,j)[3] * scale);
      }
    }
}

void SplineRigidProblem::revertSpline() {
  // Restore the points
  for (int i=0;i<=spline->dim(0);i++)
    {
    for (int j=0;j<=spline->dim(1);j++)
      {
      spline->setControl(i,j,B(i,j));
      }
    }
}

// Get the size of the vector used to drive the problem
int SplineRigidProblem::getVectorSize() {
  return 7;
}

// Get the vector for current spline state.  Vector should be properly sized
void SplineRigidProblem::makeVector(Vector &v) {
  v = Vector(7,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
}

// Get the standard deviation vector
void SplineRigidProblem::makeSDVector(Vector &v) {
  v = Vector(7,M_PI/6.0,M_PI/6.0,M_PI/6.0,0.2,0.2,0.2,0.2);
}

void SplineOptDriver::applyBestSolution() {
  problem->applyVector(method->getBestEverX());
}

bool SplineOptDriver::isFinished() {
  return method->isFinished();
}

bool SplineOptDriver::optimize(int ms) {
  // Clock the starting time
  clock_t tStart = clock();

  // Go for a while
  while (clock()-tStart < ms && !method->isFinished())
    {
    method->performIteration();
    // cout << method->getBestEverValue() << endl;
    }

  return method->isFinished();
}

void SplineOptDriver::initConjGrad(Registry *settings) {
  numFunc = new NumericalFunction(*problem,0.001);

  Vector mean(problem->getVectorSize());
  problem->makeVector(mean);

  ConjugateGradientMethod *cgm = new ConjugateGradientMethod(*numFunc,mean);
  method = cgm;
}

void SplineOptDriver::initPowell(Registry *settings) {
  Vector mean(problem->getVectorSize());
  problem->makeVector(mean);

  PowellMethod *pm = new PowellMethod(problem,mean);

  // Set the step matrix to be small
  Matrix M(mean.size(),mean.size());
  M.setAll(0.001);
  pm->setDirectionMatrix(M);

  // Set the tolerance to be low
  pm->setFTolerance(0.01);

  method = pm;
}

void SplineOptDriver::initEvolution(Registry *settings) {
  int mu = settings->getIntValue("es.mu",2);
  int lambda = settings->getIntValue("es.lambda",4);

  // Generate mean and standard deviation vectors
  int n = problem->getVectorSize();
  Vector mean(n),sd(n),delta(n),sigma0(n),sigma1(n);
  problem->makeVector(mean);
  problem->makeSDVector(sd);

  // Create a new solution space
  ss = new GaussianSS(mean,sd);

  // Create a new evolutionary strategy
  EvolutionaryStrategy *es = new EvolutionaryStrategy(*problem,*ss,mu,lambda,SELECTION_MuPlusLambda);

  // Make sure that with 90% probability the solutions will not fall farther than 0.05 away from the source
  es->computeSigmaFactorForConfidence(0.2);

  // Set the delta sigmas - the amount by which each sigma can change.  This will make it very unlikely that any
  // sigma will change by more than a factor of 2
  delta.setAll(log(2.0) / 3.0);
  es->setDeltaSigma(delta);

  // Bound the sigmas
  sigma0.setAll(0.0001);
  sigma1.setAll(1.0);
  es->setSigmaBounds(sigma0,sigma1);

  // Disable the crossover
  es->setMutationProbability(n < 4 ? 1.0 : 1.0 / n);
  //es->setGranularity(4);
  //es->useCustomCrossover(2,noCrossover); 

  // Set the first element
  es->setNth(0,mean);

  // Use the evolutionary method
  method = es;
}

SplineOptDriver::SplineOptDriver(MSpline *spline, SplineDataCache *cache, SplineImageMatcher *match, 
                                 int cpFirstU, int cpFirstV, int cpRangeU, int cpRangeV, Registry *rSettings) 
{
  this->spline = spline;
  this->splineData = cache;
  this->match = match;

  problem = new SplineOptProblem(spline,cache,match,cpFirstU,cpFirstV,cpRangeU,cpRangeV,rSettings);
  numFunc = NULL;
  ss = NULL;

  printf("Control Point %d, %d, range %d, %d\n",cpFirstU, cpFirstV, cpRangeU, cpRangeV);

  // Initialize appropriate optimizer
  //char *sMethod = rSettings->getStringValue("method","cg");
  const char *sMethod = rSettings->getStringValue("method","es");
  if (!strcmp(sMethod,"es") || !strcmp(sMethod,"evolution"))
    initEvolution(rSettings);
  else if (!strcmp(sMethod,"cg") || !strcmp(sMethod,"conjugate") || !strcmp(sMethod,"conjgrad"))
    initConjGrad(rSettings);
  else if (!strcmp(sMethod,"powell"))
    initPowell(rSettings);
  else
    throw "Unrecognized optimization method";

}

SplineOptDriver::~SplineOptDriver() {
  delete problem;
  delete method;
  if (numFunc)
    delete numFunc;
  if (ss)
    delete ss;
}


double SplineOptDriver::getBestValue() {
  return method->getBestEverValue();  
}

double SplineOptDriver::getEvalCost() {
  return problem->getEvaluationCost();
}

/**
 * This method finds a transformation that best aligns an mspline to an image
 *
void alignSplineToImage(MSpline *spline,SplineDataCache *cache, DistanceTransform *dt) {
    // Find the center of weight of all boundary points
    vector<SMLVec3f> IX;
    int i;
    
    // Go through all the voxels and add to a point list
    for(int cz=0;cz<dt->cube.size(2);cz++) {
        for(int cy=0;cy<dt->cube.size(1);cy++) {
            for(int cx=0;cx<dt->cube.size(0);cx++) {
                DataCube<short> &cube = dt->cube(cx,cy,cz);
                int bx = cx * cube.size(2);
                int by = cy * cube.size(1);
                int bz = cz * cube.size(2);

                for(int z=0;z<cube.size(2);z+=2) {
                    for(int y=0;y<cube.size(1);y+=2) {
                        for(int x=0;x<cube.size(0);x+=2) {
                            if(cube(x,y,z) < 2) {
                                SMLVec3f X(bx+x,by+y,bz+z),Y;
                                dt->IS.TransformPoint(X,Y);
                                IX.push_back(Y);
                            }
                        }
                    }
                }
            }
        }
    }

    // Compute the center of mass
    SMLVec3f IC;
    for(i=0;i<IX.size();i++) {
        IC += IX[i];
    }
    IC /= IX.size();

    // Compute the inertia matrix
    SMLMatrix3f II;
    for(i=0;i<IX.size();i++) {
        for(int j=0;j<3;j++) {
            for(int k=0;k<3;k++) {
                
            }
        }
        IC += IX[i];
    }



}*/
