#ifndef _IMATCH_H_
#define _IMATCH_H_

class PatchDataCache;
class DistanceTransform;

#include <iostream>

#include "mspline.h"
#include "Array2D.h"

/**
 * Generic spline/image matcher
 */
class SplineImageMatcher {
public:
    SplineImageMatcher();
    virtual ~SplineImageMatcher() {};

    virtual float getMatch() = 0;
    virtual void printMatchInfo(std::ostream &out) = 0;

    // Set the spline
    virtual void setSpline(MSpline *spline, SplineDataCache *cache);

protected:
    // Spline and its cache
    MSpline *spline;
    SplineDataCache *cache;

    // Time stamp array for all computed overlap measurements
    Array2D<long> tsPatch;

};

/**
 * Area normalized boundary image match measure
 */
class BoundaryImageMatcher : public SplineImageMatcher, public BoundaryMeasure {
public:

    BoundaryImageMatcher() {};
    virtual ~BoundaryImageMatcher() {};


    // Set the spline
    void setSpline(MSpline *spline, SplineDataCache *cache);

    // Compute the match
    float getMatch();

protected:
    // Stored array of computed areas and distances per patch
    Array2D<float> patchArea, patchMeasure;

    // Total area and distance
    float totalArea,totalMeasure;

    // Match one patch in the image
    float matchPatch(PatchDataCache *pdc,float &area);

    // Compute the area and mean distance
    float matchPatches(float &area);
};

/**
 * Distance transform matcher
 */
class SplineDistanceMatcher : public BoundaryImageMatcher {
public:
    SplineDistanceMatcher(DistanceTransform *distMap);
    virtual ~SplineDistanceMatcher() {};

    // Boundary and crest measurements
    float computeBoundaryMeasure(const MedialPoint &mp,int side);
    float computeCrestBoundaryMeasure(const MedialPoint &mp);

    // Print out match information
  void printMatchInfo(std::ostream &out);

private:
    DistanceTransform *distMap;
    int res,mask,sz0,sz1,sz2;

    // Compute distance at just one point
    float computeDistanceAtPoint(const MedialPoint &mp,int side);
};


/**
 * Grayscale profile matcher
 */
class ProfileMatcher : public BoundaryImageMatcher {
public:
    ProfileMatcher(GrayImage *image);
    virtual ~ProfileMatcher() {};

    // The rho constant for profile matching
    float rho;

    // Boundary and crest measurements
    float computeBoundaryMeasure(const MedialPoint &mp,int side);
    float computeCrestBoundaryMeasure(const MedialPoint &mp);

    // Print out match information
    void printMatchInfo(std::ostream &out);

private:
    // Number of samples on each side of the zero
    enum {NSAMPLES=3};

    // Gaussian derivative weight array
    float wgt[NSAMPLES];
    
    // The image
    GrayImage *image;

    // Profile computation method
    float computeProfile(const SMLVec3f &X, const SMLVec3f &N, float r);
};

/** 
 * Grayscale/binary image gradient image match
 */
class GradientImageMatcher : public BoundaryImageMatcher {
public:
    // Constructor
    GradientImageMatcher(GrayImage *image)
      {
      this->image = image;
      }

    virtual ~GradientImageMatcher() {};

    // Boundary and crest measurements
    float computeBoundaryMeasure(const MedialPoint &mp,int side)
      {
        // Compute the image gradient at the point
        SMLVec3f G;
        image->interpolateVoxelGradient(
          mp.bp[side].X[0],
          mp.bp[side].X[1],
          mp.bp[side].X[2],
          G.data_block());

        // Dot the image gradient with the m-rep boundary normal
        return dot_product(mp.bp[side].N,G);
      }

    float computeCrestBoundaryMeasure(const MedialPoint &mp)
      {
        // Compute the image gradient at the point
        SMLVec3f G;
        image->interpolateVoxelGradient(
          mp.bp[0].X[0],
          mp.bp[0].X[1],
          mp.bp[0].X[2],
          G.data_block());

        // Dot the image gradient with the m-rep boundary normal
        return dot_product(mp.bp[0].N,G);
      }

    // Print out match information
    void printMatchInfo(std::ostream &out)
      {
      float match = getMatch();
      out << "Gradient Image Match" << std::endl;
      out << "   total advection force  : " << match << std::endl;
      out << "   surface area           : " << totalArea << std::endl;
      }

private:

    // The image
    GrayImage *image;

};

/**
 *
 */
class SplineBinaryVolumeMatcher : public SplineImageMatcher {
public:
    // Volume match object
    SplineBinaryVolumeMatcher(BinaryImage *bim);
    virtual ~SplineBinaryVolumeMatcher() {};

    // Set the spline
    void setSpline(MSpline *spline, SplineDataCache *cache);

    // Compute the total match measure
    float getMatch();

    // Get a printout
    void printMatchInfo(std::ostream &out);

private:
    // Binary image
    BinaryImage *bim;

    // Array of computed volumes
    Array2D<float> volOverlap, volObject;

    // Volume of the image
    float volImage,volTotalObject,volTotalOverlap,volWhole;

    // Compute the match over one patch of the spline
    void getPatchMatch(int i,int j);

    // Compute the match on just one triangle in the image
    // Last two parameters are 6 x Voverlap, 6 x Vtotal
    void getTriangleMatch(MedialPoint &M1,MedialPoint &M2,MedialPoint &M3,float &vo6,float &vt6);
};


/*
 This method parses a 3D line in a 3D volume.  It reports all voxel values and voxel boundaries
 that it passes through.

 tVoxelStart array - every value of t at which a new voxel is entered
 Intensity array - intensity value at each of the voxels.

 The u,v vectors are in image coordinates
 */
void parseLine3D(AbstractImage3D *img,float u[3],float v[3],float tEnd,int *Iout,float *tOut,int &outSize);


/**
 * Fast trilinear interpolation for images
 */
float triLerp(float V000,float V001,float V010,float V100,
             float V011,float V101,float V110,float V111,
             float x, float y, float z);


#endif