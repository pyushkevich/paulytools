#include "mspline.h"
#include "imaging.h"
#include "imatch.h"
#include <fvec.h>

SplineImageMatcher::SplineImageMatcher() {
    spline = NULL;
    cache = NULL;
}

void SplineImageMatcher::setSpline(MSpline *spline, SplineDataCache *cache) {
    this->spline = spline;
    this->cache = cache;
    tsPatch.resize(cache->patch.width(),cache->patch.height());
    tsPatch.setAll(-3);
}

SplineDistanceMatcher::SplineDistanceMatcher(DistanceTransform *distMap) {
    this->distMap = distMap;
    
    // Resolution and mask
    res = 1 << distMap->LDCUBE;
    mask = res-1;

    // Cube dimensions
    sz0 = distMap->size(0);
    sz1 = distMap->size(1);
    sz2 = distMap->size(2);
}

void BoundaryImageMatcher::setSpline(MSpline *spline, SplineDataCache *cache) {
    SplineImageMatcher::setSpline(spline,cache);
    patchArea.resize(cache->patch.width(),cache->patch.height());
    patchMeasure.resize(cache->patch.width(),cache->patch.height());
}

extern int round(float f);

/**
 * Compute match for a vector
 */
float SplineDistanceMatcher::computeDistanceAtPoint(const MedialPoint &mp,int side) {
    // Get the distance value   
    return fabs(distMap->interpolateVoxel(mp.bp[side].X[0],mp.bp[side].X[1],mp.bp[side].X[2]));
}

/***********************************
    Equivalent to the following but using only 4 multiply operations.
    (1-x)*(1-y)*(1-z)*V000 + (1-x)*(1-y)*(  z)*V001 +
    (1-x)*(  y)*(1-z)*V010 + (  x)*(1-y)*(1-z)*V100 +
    (1-x)*(  y)*(  z)*V011 + (  x)*(1-y)*(  z)*V101 +
    (  x)*(  y)*(1-z)*V110 + (  x)*(  y)*(  z)*V111;
 ************************************/
float triLerp(float V000,float V001,float V010,float V100,
              float V011,float V101,float V110,float    V111,
              float x,   float  y,  float z)    
{
    F32vec4 r0,r1,r2,r3,r4,r5,r6,r7,r8;
    
                                                            //  3       2       1       0

    r0 = _mm_loadu_ps(&x);                                  //  ?       z       y       x
    r1 = _mm_shuffle_ps(r0,r0,0x19);                        //  x       y       z       y
    r2 = _mm_shuffle_ps(r0,r0,0x62);                        //  y       z       x       z
    r3 = _mm_mul_ps(r1,r2);                                 //  xy      yz      zx      yz
    r4 = _mm_shuffle_ps(r0,r0,0x84);                        //  z       x       y       x           
    r5 = _mm_mul_ps(r4,r3);                                 //  xyz     xyz     xyz     xyz
    r4 = _mm_sub_ps(r1,r3);                                 //  x-xy    y-yz    z-zx    y-yz
    r1 = _mm_shuffle_ps(r3,r3,0x71);                        //  zx      yx      yz      xz
    r3 = _mm_sub_ps(r1,r5);                                 //  zx-xyz  yx-xyz  yz-xyz  xz-xyz
    r1 = _mm_shuffle_ps(r4,r4,0xc7);                        //  x-xy    y-yz    z-zx    x-xy
    r6 = _mm_sub_ps(r1,r3);                                 //  x-xy-xz y-yz-yx z-zx-yz x-xy-xz
                                                            //   +xyz    +xyz    +xyz    +xyz
    r1 = _mm_add_ps(r6,r4);                                 //  ...     ...     ...     x+y-yz-xy-xz+xyz
    r8 = _mm_set_ss(1.0f);                                  //  0       0       0       1
    r8 = _mm_sub_ss(r8,r2);                                 //  0       0       0       1-z 
    r8 = _mm_sub_ss(r8,r1);                                 //  0       0       0       1-z-x-y+xy+yz+xz-xyz
    r5 = _mm_shuffle_ps(r5,r3,0xa0);                        //  yx-xyz  yx-xyz  xyz     xyz
    r3 = _mm_shuffle_ps(r3,r5,0x2d);                        //  xyz     yx-xyz  zx-xyz  yz-xyz
    r8 = _mm_shuffle_ps(r8,r6,0x54);                        //  z-zx-.. z-zx-.. 0       1-z-x-y+xy+yz+xz-xyz
    r8 = _mm_shuffle_ps(r8,r6,0xe8);                        //  x-xy..  y-yz..  z-zx..  1-z-x-y+xy+yz+xz-xyz

    r1 = _mm_loadu_ps(&V000);                               //  V100    V010    V001    V000
    r2 = _mm_loadu_ps(&V011);                               //  V111    V110    V101    V011
    r1 = _mm_mul_ps(r1,r8);                                 
    r2 = _mm_mul_ps(r2,r3);                                 

    // Now add up all elements
    r1 = _mm_add_ps(r1,r2);
    r2 = _mm_shuffle_ps(r1,r1,0xb1);
    r1 = _mm_add_ps(r1,r2);
    r2 = _mm_shuffle_ps(r1,r1,0x4e);
    r1 = _mm_add_ps(r1,r2);

    // The result is in r1
    return r1[0];
}

float BoundaryImageMatcher::getMatch() {
    
    // Initialize totals
    totalArea = 0;
    totalMeasure = 0;

    // Make sure that all the patches have been freshly computed
    cache->refreshMedialAllPatches(cache->getMaxLevel());
    cache->refreshBoundaryAllPatches(cache->getMaxLevel());
    
    // Check timestamps over the image
    for(int i=0;i<cache->patch.width();i++) {
        for(int j=0;j<cache->patch.height();j++) {

            // Get the timestamp if the patch and compare it to the volume timestamps
            long ts = cache->getPatchTimeStamp(i,j,cache->getMaxLevel(),SplineDataCache::BOUNDARY);
            if(ts > tsPatch(i,j)) {

                // Compute the volume match over the patch
                patchMeasure(i,j) = matchPatch(cache->patch(i,j),patchArea(i,j));

                // Update the time stamp
                tsPatch(i,j) = ts;
            }

            // Add the volumes
            totalArea += patchArea(i,j);
            totalMeasure += patchMeasure(i,j);
        }
    }

    // Return the mean square distance
    return totalMeasure / totalArea;
}

void SplineDistanceMatcher::printMatchInfo(std::ostream &out) 
{
  float match = getMatch();
  out << "Distance Match" << std::endl;
  out << "   mean square distance  : " << match << std::endl;
  out << "   surface area          : " << totalArea << std::endl;
}

/**
 * Integrate spline to distance transform match over a single patch
 */
float BoundaryImageMatcher::matchPatch(PatchDataCache *pdc,float &area) {
    return pdc->integrateBoundaryMeasure(*this,area);
}

float SplineDistanceMatcher::computeBoundaryMeasure(const MedialPoint &mp,int side) {
    return computeDistanceAtPoint(mp,side);
}

float SplineDistanceMatcher::computeCrestBoundaryMeasure(const MedialPoint &mp) {
    return computeDistanceAtPoint(mp,0);
}



ProfileMatcher::ProfileMatcher(GrayImage *image) {
    this->image = image;

    // Set the default rho value
    rho = 0.25;
    
    // Gaussian derivative weight function
    wgt[2] = 0.194276f;      // -1.5 
    wgt[1] = 0.24197f;       // -1.0
    wgt[0] = 0.176033f;      // -0.5
}

inline float ProfileMatcher::computeProfile(const SMLVec3f &X, const SMLVec3f &N, float r) {
    float total = 0;
    SMLVec3f X1 = X, X2 = X;

    // Compute the step vector
    SMLVec3f step = N * (0.5f * rho * r);

    for(int i=0;i<NSAMPLES;i++) {
        X1 += step;
        X2 -= step;
        total -= wgt[i] * image->interpolateVoxel(X1[0],X1[1],X1[2]);
        total += wgt[i] * image->interpolateVoxel(X2[0],X2[1],X2[2]);
    }

    return  -total;
}

// Boundary and crest measurements
float ProfileMatcher::computeBoundaryMeasure(const MedialPoint &mp,int side) {
    return computeProfile(mp.bp[side].X,mp.bp[side].N,mp.R());
}

float ProfileMatcher::computeCrestBoundaryMeasure(const MedialPoint &mp) {
    return computeProfile(mp.bp[0].X,mp.bp[0].N,mp.R());
}

// Print out match information
void ProfileMatcher::printMatchInfo(std::ostream &out) {
    out << "Profile Match" << endl;
    out << "   integrated image grad : " << getMatch() << endl;
    out << "   surface area          : " << totalArea << endl;
}












SplineBinaryVolumeMatcher::SplineBinaryVolumeMatcher(BinaryImage *bim) {
    // Store the image
    this->bim = bim;

    // Compute the overall volume of the image
    int nvox = 0;

    for(int i=0;i<bim->cube.sizeOfCube();i++) {
        DataCube<unsigned char> &ic = bim->cube(i);
        for(int z=0;z<ic.size(2)-1;z++)
            for(int y=0;y<ic.size(1)-1;y++)
                for(int x=0;x<ic.size(0)-1;x++)
                    if(ic(x,y,z))
                        nvox ++;
    }

    // Compute the voxel size based on the transform
    SMLVec3f unit(1.0f,1.0f,1.0f);
    SMLVec3f voxSize = TransformVector(bim->IS,unit);
    float voxVol = voxSize[0] * voxSize[1] * voxSize[2];
    
    // Compute the volume of the image
    volImage = voxVol * nvox;
    volWhole = voxVol * bim->size(0) * bim->size(1) * bim->size(2);
}

void SplineBinaryVolumeMatcher::setSpline(MSpline *spline, SplineDataCache *cache) {
    SplineImageMatcher::setSpline(spline,cache);
    volObject.resize(cache->patch.width(),cache->patch.height());
    volOverlap.resize(cache->patch.width(),cache->patch.height());
}

void SplineBinaryVolumeMatcher::printMatchInfo(std::ostream &out) {
    float ovl = -getMatch();

    out << "Volume Match" << endl;
    out << "   volume of image box      : " << volWhole << endl;
    out << "   volume of image data     : " << volImage << endl;
    out << "   volume of spline object  : " << volTotalObject << endl;
    out << "   volume of intersection   : " << volTotalOverlap << endl;
    out << "   relative overlap         : " << ovl << endl;
}

float SplineBinaryVolumeMatcher::getMatch() {
    // Total volumes
    volTotalObject = 0.0f;
    volTotalOverlap = 0.0f;

    // Make sure that all the patches have been freshly computed
    cache->refreshMedialAllPatches(cache->getMaxLevel());
    cache->refreshBoundaryAllPatches(cache->getMaxLevel());
    
    // Check timestamps over the image
    for(int i=0;i<cache->patch.width();i++) {
        for(int j=0;j<cache->patch.height();j++) {
            // Get the timestamp if the patch and compare it to the volume timestamps
            long ts = cache->getPatchTimeStamp(i,j,cache->getMaxLevel(),SplineDataCache::BOUNDARY);
            if(ts > tsPatch(i,j)) {
                // Compute the volume match over the patch
                getPatchMatch(i,j);

                // Update the time stamp
                tsPatch(i,j) = ts;
            }

            // Add the volumes
            volTotalObject += volObject(i,j);
            volTotalOverlap += volOverlap(i,j);
        }
    }

    // Return the intersection of volumes divided by the union
    return - volTotalOverlap / (volTotalObject + volImage - volTotalOverlap);
}

void SplineBinaryVolumeMatcher::getPatchMatch(int i,int j) {
    // Get the patch
    PatchDataCache *pdc = cache->patch(i,j);
    SplineGridDefinition *grid = pdc->getGrid();
    
    // An array analogous to MP of volume overlaps
    volOverlap(i,j) = 0;
    volObject(i,j) = 0;

    // Patch resolution
    int w = pdc->MP.width();
    int h = pdc->MP.height();

    // Parse all interior squares in the patch
    for(int v=0;v<h-1;v++) {
        for(int u=0;u<w-1;u++) {
            // Check if the square is inside the object
            if(pdc->in(u,v) & pdc->in(u+1,v) & pdc->in(u,v+1) & pdc->in(u+1,v+1) & 0x01) {
                getTriangleMatch(pdc->MP(u,v),pdc->MP(u,v+1),pdc->MP(u+1,v),
                    volOverlap(i,j),volObject(i,j));
                getTriangleMatch(pdc->MP(u,v+1),pdc->MP(u+1,v+1),pdc->MP(u+1,v),
                    volOverlap(i,j),volObject(i,j));
            }
        }
    }

    // Parse all other triangles
    for(int k=0;k<pdc->idxTrimStripIdx;k+=3) {
        getTriangleMatch(pdc->MP(pdc->idxTrimStrip[k]),pdc->MP(pdc->idxTrimStrip[k+1]),
            pdc->MP(pdc->idxTrimStrip[k+2]),volOverlap(i,j),volObject(i,j));
    }

    // Scale the volumes
    volOverlap(i,j) *= 0.166667f;
    volObject(i,j) *= 0.166667f;
}

/*
 This method parses a 3D line in a 3D volume.  It reports all voxel values and voxel boundaries
 that it passes through.

 tVoxelStart array - every value of t at which a new voxel is entered
 Intensity array - intensity value at each of the voxels.

 The u,v vectors are in image coordinates
 */
void parseLine3D(AbstractImage3D *img,float u[3],float v[3],float tEnd,int *Iout,float *tOut,int &outSize) {
    
    // The starting t position
    float t0 = 0;
    int i;
    int X[3], S[3];
    float tMax[3],tDelta[3],w[3];

    // Output array index
    outSize = 0;

    // Initialize the step sizes
    for(i=0;i<3;i++) {
        tDelta[i] = fabs(1.0f / v[i]);
        w[i] = u[i] + tEnd * v[i];
    }

    // Figure out the starting position
    for(i=0;i<3;i++) {
        if(u[i] < 0) {
            if(v[i] > 0) {
                t0 = vnl_math_max(t0,-u[i] * tDelta[i]);
            }
            else {
                return;
            }
        }
        else if(u[i] > img->size(i)) {
            if(v[i] < 0) {
                t0 = vnl_math_max(t0,(u[i] - img->size(i)) * tDelta[i]);
            }
            else {
                return;
            }
        }
    }

    // Figure out the ending position
    for(i=0;i<3;i++) {
        if(w[i] < 0) {
            if(v[i] < 0) {
                tEnd = vnl_math_min(tEnd,w[i] * tDelta[i]);
            }
            else {
                return;
            }
        }
        else if(w[i] > img->size(i)) {
            if(v[i] > 0) {
                tEnd = vnl_math_min(tEnd,-(w[i] - img->size(i)) * tDelta[i]);
            }
            else {
                return;
            }
        }
    }

    // Make sure there is anything to trace
    if(t0 >= tEnd)
        return;

    // Adjust the starting position to the u,v vector.  We have to be careful here not to get 
    // rounding error and hence negative values
    if(t0 > 0) {
        for(i=0;i<3;i++) {
            u[i] += t0 * v[i];
            if(u[i] < 0.0)
                u[i] = 0.0;
            if(u[i] >= img->size(i)) 
                u[i] = img->size(i)-0.00001f;
        }
    }

    // Initialize the data
    for(i=0;i<3;i++) {
        X[i] = (int)u[i];
        if(v[i] == 0) { 
            S[i] = 0;
            tMax[i] = 1.0e32f;
        }
        else if(v[i] < 0) {
            S[i] = -1;
            tMax[i] = t0 + (u[i] - X[i]) * tDelta[i];
        }
        else {
            S[i] = 1;
            tMax[i] = t0 + (X[i]+1-u[i]) * tDelta[i];       
        }
    }

    // Get the starting intensity
    int iLast = img->getVoxelNBC(X[0],X[1],X[2]);
    
    // Store the initial t and intensity
    Iout[outSize] = iLast;
    tOut[outSize++] = t0;

    // Which max is the smallest?
    int idx = 1;
    while(1) {
        int iMax = tMax[0] < tMax[1] ? (tMax[0] < tMax[2] ? 0 : 2) : (tMax[1] < tMax[2] ? 1 : 2);
        if(tMax[iMax] < tEnd) {

            // Go to the next voxel
            X[iMax] += S[iMax];

            // Get the intensity at the new voxel
            int iNew = img->getVoxelNBC(X[0],X[1],X[2]);

            // Save the position if intensity changed
            if(iNew != iLast) {
                // Store the initial t and intensity
                Iout[outSize] = iNew;
                tOut[outSize++] = tMax[iMax];
                iLast = iNew;
            }

            // Update the t value
            tMax[iMax] += tDelta[iMax];
        }
        else {
            tOut[outSize++] = tEnd;
            return;
        }
    }
}

void SplineBinaryVolumeMatcher::getTriangleMatch(MedialPoint &M0,MedialPoint &M1,MedialPoint &M2,
                                            float &vOverlap6,float &vTotal6) 
{
    // Medial triangle center point
    SMLVec3f CM = (M0.X() + M1.X() + M2.X()) * 0.333333f;
    
    // Same point, in image coordinates
    SMLVec3f ICM = TransformPoint(bim->SI,CM);

    // Distances between vectors adj points
    SMLVec3f d1 = M1.X() - M0.X();
    SMLVec3f d2 = M2.X() - M0.X();

    // Zero volume overlap
    // vOverlap6 = vTotal6 = 0;

    float sign = 1.0;
    for(int d=0;d<2;d++) {
        // Boundary points
        BoundaryPoint &bp0 = M0.bp[d];
        BoundaryPoint &bp1 = M1.bp[d];
        BoundaryPoint &bp2 = M2.bp[d];

        // Compute the dipstick vector
        SMLVec3f CB = (bp0.X + bp1.X + bp2.X) * 0.33333333f;
        SMLVec3f D = CB - CM;

        // Map vector D to image coordinates
        SMLVec3f ID = TransformVector(bim->SI,D);

        // 3D-Rasterize the vector
        int iVal[100],nVal;
        float tVal[100];
        parseLine3D(bim,ICM.data_block(),ID.data_block(),1.0,iVal,tVal,nVal);
        
        // The volume of a triangular prism is given by
        // 6V = t[(d1 x d2) . (u0+u1+u2)] + t^2[d1.(u2xu1)+d2.(u0xu1)+d1.(u2xu0)] + t^3[u0.u1xu2]
        SMLVec3f u0 = bp0.N * M0.R();
        SMLVec3f u1 = bp1.N * M1.R();
        SMLVec3f u2 = bp2.N * M2.R();
        SMLVec3f u0u1 = vnl_cross_3d(u0,u1);
        SMLVec3f u1u2 = vnl_cross_3d(u1,u2);
        SMLVec3f u2u0 = vnl_cross_3d(u2,u0);
        
        // Now, 6V = At + Bt^2 + Ct^3
        float A = sign * dot_product(vnl_cross_3d(d1,d2),u0+u1+u2);
        float B = sign * dot_product(d1,u2u0 - u1u2) + dot_product(d2,u0u1);
        float C = sign * dot_product(u0,u1u2);

        // Starting point on the profile depends on the first voxel
        int i0 = iVal[0] ? 0 : 1;       
        for(int i=i0;i<nVal-1;i+=2) {
            vOverlap6 += tVal[i+1] * (A + tVal[i+1] * (B + tVal[i+1] * C));
            vOverlap6 -= tVal[i] * (A + tVal[i] * (B + tVal[i] * C));
        }

        // The total volume
        vTotal6 += (A + B + C);

        // Change the sign
        sign -= 2.0f;
    }
}


