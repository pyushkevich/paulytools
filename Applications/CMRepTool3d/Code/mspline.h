#ifndef M_SPLINE_H
#define M_SPLINE_H

// Includes
#include <cmath>
#include <smlmath.h>
#include <vector>
#include "array2d.h"
#include "align.h"

class DynamicBSpline2D;


/**
 * A control point in the grid
 */
class BSplineControl {
public:
    // The value at this point
    float x;

    // The values in the 4x4 neighborhood of control points
    SMLMatrix4f nbr;

    // The equivalent of nbr, but transposed and stored for fast access
    float *mmNbrArray,*mmNbrArrayT;

    // This MGM matrix is used for uniform B-splines.  It is M * nbr * Mt
    // SMLMatrix4f MGM;

    BSplineControl() {
        nbr.fill(0);
        mmNbrArray = (float *)_aligned_malloc(sizeof(float)*16,32);
        mmNbrArrayT = (float *)_aligned_malloc(sizeof(float)*16,32);
    }

    ~BSplineControl() {
        _aligned_free(mmNbrArray);
        _aligned_free(mmNbrArrayT);
    }
};

typedef Array2D<BSplineControl> BSplineControlArray;


/**
 * Spline Grid Definition
 * A grid definition is basically an array of u and v values.  
 */
class SplineGridDefinition {
    // U and V values at grid positions
    float *gridValues[2];

    // Knot value at each grid position
    int *knotValues[2];

    // Grid index for each patch
    int *patchIndex[2];

    // Dimensions of the grid
    int dim[2], size[2];

    // For each grid position, there is a vector of precomputed basis function weights
    SMLVec4f **W[2];

    // There are also weight vectors for the patches
    SMLVec4f ***patchWeight[2];

    // Time stamp for the weights
    int tsWeights;

    // Hidden constructor
    SplineGridDefinition(); 

    // Update the weights for a given spline
    void updateWeights(DynamicBSpline2D *spline);

    // Allow the spline access to this class
    friend class DynamicBSpline2D;
    friend class MSpline;
public:
    // Destructor
    ~SplineGridDefinition();

    // Get the UV coordinate at an index
    float uv(int dim,int index) const {
        return gridValues[dim][index];
    }

    // Get the UV coordinate at an index in a patch
    float uv(int dim,int iPatch,int iGrid) const {
        return gridValues[dim][patchIndex[dim][iPatch]+iGrid];
    }


    // Get the knot value of a value
    int knot(int dim,int index) const {
        return knotValues[dim][index];
    }

    // Get the patch index for a grid value
    int patch(int dim,int index) const {
        return knotValues[dim][index]-3;
    }

    // Find the index at which i-th patch begins.  (Can pass 0..m-2, m-2 for end of last patch)
    int patchStart(int dim,int iPatch) const {
        return patchIndex[dim][iPatch];
    }

    // Find the index of a u or v value
    int getIndexForUV(int dim,float uv);

    // Get the size of the grid
    int getSize(int dim) {
        return size[dim];
    }

    // Create a new grid definition where each patch has a 2^k resolution
    // and knot values are interpolated precisely.
    static SplineGridDefinition *newSGDForCache(DynamicBSpline2D *spline,int resLevel);
};

/**
 * Boundary point data
 */
struct BoundaryPoint {
    // The boundary position and its normal
    SMLVec3f X,N,Xu,Xv,Nu,Nv;

    // Two principal curvatures
    float kappa1,kappa2;
};

// This is a point on a 3D surface 
/*
struct SurfacePoint {
    // The partial derivatives
    SMLVec3f X,Xu,Xv,Xuu,Xuv,Xvv,N;
    float E,F,G,L,M,N;
};

/**
 * This structure represents a point sampled from the medial spline
 */
ALIGN32_PRE struct MedialPoint {
    // The function (xyzr) and its partial derivatives
    SMLVec4f F,Fu,Fv,Fuu,Fuv,Fvv;

    // The normal to the spline (normalized and unit)
    SMLVec4f F3,F3raw;

    // Grad R vector and other assorted geometrical information
    SMLVec3f GradR;

    // The length of the normal vector
    float normN;

    // The elements of the first fundamental form on the surface
    float IE,IF,IG;

    // The principal curvatures
    float kappa1,kappa2;

    // The term sin(theta)^2, that should be greater than zero
    float sinTheta2;

    // Boundary point data
    BoundaryPoint bp[2];

    // The u,v position
    float u,v;

    // References that access the data as 3-vectors
    SMLVec3f &X(), &Xu(), &Xv(), &Xuu(), &Xuv(), &Xvv(), &N(), &NRaw();
    float &R(), &Ru(), &Rv(), &Ruu(), &Ruv(), &Rvv();

    const SMLVec3f &X() const;
    const SMLVec3f &Xu() const;
    const SMLVec3f &Xv() const;
    const SMLVec3f &Xuu() const;
    const SMLVec3f &Xuv() const;
    const SMLVec3f &Xvv() const;
    const SMLVec3f &N() const;
    const SMLVec3f &NRaw() const;
    float R() const;
    float Ru() const;
    float Rv() const;
    float Ruu() const;
    float Ruv() const;
    float Rvv() const;

    // This is needed for callback functions
    void *data;
} ALIGN32_POST;

/** 
 * Two dimensional rectangular uniform B-spline
 */
class DynamicBSpline2D {
public:
    // Create a new spline with default control points.  The spline will have m+1 and n+1 control points.  
    // Each control point has d dimensions (3 for a 3D spline)
    DynamicBSpline2D(int m,int n,int d,int zeroKnotMargin=2);
    ~DynamicBSpline2D();

    // This method sets one of the control values
    void setControl(int i,int j,int k,float value);
    void setControl(int i,int j,const SMLVec3f X) {
        setControl(i,j,0,X[0]);
        setControl(i,j,1,X[1]);
        setControl(i,j,2,X[2]);
    }
    void setControl(int i,int j,const SMLVec4f X) {
        setControl(i,j,0,X[0]);
        setControl(i,j,1,X[1]);
        setControl(i,j,2,X[2]);
        setControl(i,j,3,X[3]);
    }

    // Get a control point value
    float getControl(int i,int j,int k) {
        return P[k](i,j).x;
    }

    // Get a control as a 3-vector
    SMLVec3f getControl(int i,int j) {
        SMLVec3f XC;
        XC[0] = P[0](i,j).x;
        XC[1] = P[1](i,j).x;
        XC[2] = P[2](i,j).x;
        return XC;
    }

    // Get a control as a 3-vector
    SMLVec4f getControl4f(int i,int j) {
        SMLVec4f XC;
        XC[0] = P[0](i,j).x;
        XC[1] = P[1](i,j).x;
        XC[2] = P[2](i,j).x;
        XC[3] = P[3](i,j).x;
        return XC;
    }

    // Get one of the dimensions of the spline
    int dim(int i) {
        return m[i];
    }

    // Return the u or v value for a given knot
    float getKnotUV(int coordinate,int knotIndex) {
        return knots[coordinate][knotIndex];
    }

    // This method samples the B-spline over the whole grid.
    //void compute(int uJet,int vJet,int iControl,Array2D<float> &output);
    //void computeBox(int uJet,int vJet,int iControl,Array2D<float> &output,int u0,int v0,int u1,int v1);

    // This method samples the B-spline over the whole grid, 3 control values at a time
    // iJet : index into the jet (0 = X, 1 = Xu, 2 = Xv, ... 5 = Xvv)
    // iControl : index into the control array (e.g. 0 = X, 1 = Y, 2 = Z)
    //void compute(int uJet,int vJet,int iControl,Array2D<SMLVec3f> &output);
    //void computeBox(int uJet,int vJet,int iControl,Array2D<SMLVec3f> &output,int u0,int v0,int u1,int v1);

    // Compute a patch of the spline
    //void computePatch(int uJet,int vJet,int iControl,int iPatch,int jPatch,Array2D<SMLVec3f> &output,int stepLog=0);

    // Compute a range of u,v values affected by a rectangle of control points
    // void findUVBox(int pu,int pv,int pw,int ph,int &u0,int &v0,int &u1,int &v1);

    // This method pushes the control points to move the spline
    void push(int uKnot,int vKnot,float u,float v,int iControlFirst,int iControlLast,float *vector);
    
    // This method pushes the control points to move the spline
    void push(const SplineGridDefinition &grid,int iGrid,int jGrid,int iControlFirst,int iControlLast,float *vector) {
        push(grid.knot(0,iGrid),grid.knot(1,jGrid),grid.uv(0,iGrid),grid.uv(1,jGrid),iControlFirst,iControlLast,vector);
    }

    // This method return the control point array for the first three dimensions.
    void getControlArray(Array2D<SMLVec3f> &out,int startIdx=0);
    void getControlArray(Array2D<SMLVec4f> &out,int startIdx=0);

    // Get the timestamp for a patch
    long getPatchTS(int i,int j) {
        return tsPatch(i,j);
    }

    // Get the timestamp for all patches, based on control point state
    long getControlTimeStamp() {
        return tsControl;
    }

    // This method computes the basis function Nip at a specific value.  Third parameter is an array of length 3, in it the 
    // 0th,1st and 2nd derivative of the basis function is stored
    void basisJet(int dim,int i,float u,SMLVec4f *W);

    // Interpolate the spline or its derivatives over a grid with variable spacing.
    // The out array must match the (dimensions of the grid)*outStride!!!
    void interpolateGrid(const SplineGridDefinition &grid,
                                       int uvFirst[2],int uvLast[2],int uvStep[2],int uJet,int vJet,
                                       int iControlFirst,int iControlLast,
                                       float *out);

    // Same as the above but just one point
    void interpolateGridPoint(const SplineGridDefinition &grid,int u,int v,int uJet,int vJet,int iControlFirst,int iControlLast,float *out);

    /**
     * Interpolate a point that is not part of a grid.  Pass in the patch index as well as the weights, computed by calling
     * basisJet.
     */
    void interpolatePoint(int iPatch,int jPatch,SMLVec4f *Wu,SMLVec4f *Wv,int uJet,int vJet,int iControlFirst,int iControlLast,float *out);

protected:
    // For uniform knot sequences this matrix is used
    SMLMatrix4f *UB;
    SMLVec4f *jetMul[2];

    // Dimensions of the spline
    int m[2];

    // Timestamps for each patch in the system
    Array2D<long> tsPatch;

    // Upper Bound timestamp: changed with each control point
    long tsControl;

    // Array of control points, sorted first by the coordinate and then by index
    BSplineControlArray *P;

    // Knot sequences in u and v
    float *knots[2];

    // Whether or not the knot sequences are 'regular'
    bool regKnots;
};

class MSpline : public DynamicBSpline2D {
private:
    // These four methods are used for crest search.  They interpolate just the values needed to compute
    // the crest and nothing else
    void interpolateMedialCrestD1(int iPatch,int jPatch,SMLVec4f *Wu,SMLVec4f *Wv,MedialPoint &mp);
    void interpolateMedialCrestD2(int iPatch,int jPatch,SMLVec4f *Wu,SMLVec4f *Wv,MedialPoint &mp);
    void interpolateMedialCrestD1Missing(int iPatch,int jPatch,SMLVec4f *Wu,SMLVec4f *Wv,MedialPoint &mp);
    void interpolateMedialCrestD2Missing(int iPatch,int jPatch,SMLVec4f *Wu,SMLVec4f *Wv,MedialPoint &mp);
    friend class PatchDataCache;

public:
    // Constructor for the MSpline
    MSpline(int m,int n);

    // This is the fastest way to interpolate the patch
    void interpolateMedialSurfacePatch(const SplineGridDefinition &grid,int iPatch,int jPatch,
                                       int step,Aux2D<MedialPoint> *MP);

    // This method computes the second derivatives for a patch
    void interpolateMedialSurfacePatch02(const SplineGridDefinition &grid,int iPatch,int jPatch,
                                       int step,Aux2D<MedialPoint> *MP);

    // This interpolates the 2-jet at a medial location, stores result in a MedialPoint structure
    void interpolateMedialPoint02(int iPatch,int jPatch,SMLVec4f *Wu,SMLVec4f *Wv,MedialPoint &mp);

    // This interpolates the 2-jet at a medial location, stores result in a MedialPoint structure
    void interpolateMedialPoint02(const SplineGridDefinition &grid,int u,int v,MedialPoint &mp);

    // Interpolate the boundary up to the second order
    void computeCurvatures(MedialPoint &MP,bool isCrest=false);

    // Boundary interpolation routines that don't compute anything but the boundary
    bool interpolateBoundary(MedialPoint &MP);
    void interpolateBoundaryCrest(MedialPoint &MP);
};

// This is an interface for measuring something on a medial boundary
class BoundaryMeasure {
public:
    virtual float computeBoundaryMeasure(const MedialPoint &mp,int side) = 0;
    virtual float computeCrestBoundaryMeasure(const MedialPoint &mp) = 0;
};

// This is an interface for measuring something on a medial sheet
class MedialMeasure {
public:
    virtual float computeMedialMeasure(const MedialPoint &mp) = 0;
};

// This object represents data associated with a triangle in a patch
class PatchTriangle {
public:
    // Area of the boundary triangle as a vector
    SMLVec3f NBnd[2],NMed;

    // Actual area of the triangle
    float areaBnd[2],areaMed;

    // The vertices that form this triangle
    int v[3];

    // Is the triangle contained by the trim curve?
    int in;
};

/**
 * The spline data cache is a dynamic storage for spline data at a 'pyramid' of levels
 * Each type of data is cached separately
 */
class PatchDataCache {
public:
    // Enumeration of types of data
    enum Types {MEDIAL=0, BOUNDARY=1, AW_MEDIAL=2, AW_BOUNDARY=3, CURVATURES = 4, NUM_TYPES=5};

    // Constructor for the data cache
    PatchDataCache(MSpline *spline,SplineGridDefinition *grid,int iPatch,int jPatch,int maxLevel=6);

    // Destructor
    virtual ~PatchDataCache();

    // Get the grid for this cache
    SplineGridDefinition *getGrid() {
        return grid;
    }

    // Get the maximum level used
    int getMaxLevel() {
        return maxLevel;
    }

    // Get the timestamp at a patch at a level
    long getTimeStamp(int level,int type) {
        return tsPatch[type][level];
    }

    // Patch index
    int iPatch,jPatch;


    // This is the combination of X and R derivatives, conveniently stored in an aligned fashion
    // Aux2D<SMLVec4f> F,Fu,Fv,Fuu,Fuv,Fvv,Nraw,N;

    // These are the curvatures
    // Aux2D<float> K[2],H[2];

    // The medial surface data - may not be fresh!
    // Aux2D<SMLVec3f> X,Xu,Xv,N,Nraw,Xuu,Xuv,Xvv;

    // The radius function data - may not be fresh!
    // Aux2D<float> R,Ru,Rv,Ruu,Ruv,Rvv;

    // The sine square of theta, this value must be positive or zero, and where it is negative, we have a problem
    // Aux2D<float> sinTheta2;

    // In/out indicator array.  Value 0 means out, 1 means in, 2 means edge has been found here.
    Aux2D<MedialPoint> MP;
    Aux2D<unsigned char> in;

    // A list of triangles formed by the patch, also the areas of these triangles
    PatchTriangle *PT;

    // Number of these triangles, number of regular triangles and number of trim triangles
    int nTriangles, nTrianglesTrim, nTrianglesReg;

    // The boundary surface patch data - may not be fresh!
    // Aux2D<SMLVec3f> XB[2],NB[2];

    // Array of trimmed vertices (same as checking in(i) for each point, but faster)
  std::vector<int> idxTrimmed;

    // The number of trimmed interior points
    int idxTrimmedInteriorSize;

    // Area weights for each point - used to integrate over the boundary
  std::vector<float> awBnd[2],awMed,awMedTrimmed;

    // Sum of area weights of the patch, boundary and medial.  Equal to 6*area.
    float awBndTotal,awMedTotal,awMedTrimmedTotal;

    // Trim curve segments' start and end indices (assuming not too many per patch)
    int idxTrimCurves[16],idxTrimCurvesIdx;

    // Trim curve strips - vertex indices for a bunch of triangles
    int *idxTrimStrip,idxTrimStripIdx;

    // Length of the crest curve
    int getCrestSize() {
        return iCrest;
    }

    // This method integrates any boundary measure along the patch.  The area returned is the total area
    float integrateBoundaryMeasure(BoundaryMeasure &measure,float &outArea);

    // This method integrates any boundary measure along the patch.  The measure is computed over the whole
    // patch, not just the trimmed region
    float integrateMedialMeasure(MedialMeasure &measure,float &outArea);
    float integrateTrimmedMedialMeasure(MedialMeasure &measure,float &outArea);

    // Refresh the second order boundary information
    void refreshCurvatures(int level);

private:
    // The maximum resolution level
    int maxLevel;

    // The spline that is our data source
    MSpline *spline;

    // A grid of points (with weights)
    SplineGridDefinition *grid;

    // Time stamps for each patch at all levels of resolution
    long *tsPatch[NUM_TYPES];

    // The index into the crest array of the current crest point
    int iCrest;

    // This method is used to find the crest along a line
    void findCrest(int dir,int iIn,int jIn,int iOut,int jOut);

    // This class is used in crest finding
    struct CrestSearch {
        // Current u, v position
        float u, v;

        // Search direction
        int dir;

        // Weights for the computation
        SMLVec4f Wu[3],Wv[3];
    };

    // This is the Newton-Raphson method from NRC
    void rtsafe(CrestSearch &crestSearch,float x1, float x2,float xacc);

    // Update the data at a patch at a given level
    void refreshMedial(int level);
    void refreshBoundary(int level);
    
    // Refresh the weights
    void refreshMedialAreaWeights(int level);
    void refreshBoundaryAreaWeights(int level);

    // Dimension and step for trim curve search
    int last, step;

    // Search for medial trim curve on the patch
    void computeTrimCurve(int level);

    // Marching methods for the trim curve finding
    void startMarch(int i,int j,int inSide);
    void doMarch(int i,int j,int inSide);

    // Some crest search computations
    void computeCrestFunction(CrestSearch &cs,float coordVal,float *f,float *df);
    void computeCrestMedialBoundary(CrestSearch &cs);

    friend class SplineDataCache;

    // Arrays used for boundary and medial computations
    Aux2D<float> MM,BM[2];
};

/**
 * The spline data cache is a dynamic storage for spline data at a 'pyramid' of levels
 * Each type of data is cached separately
 */
class SplineDataCache {
public:
    // Enumeration of types of data
    enum Types {MEDIAL=0, BOUNDARY=1, NUM_TYPES=2};

    // Constructor for the data cache
    SplineDataCache(MSpline *spline,int maxLevel=6);

    // Destructor
    virtual ~SplineDataCache();

    // Update the data at a patch at a given level
    bool refreshMedialPatch(int i,int j,int level);
    bool refreshBoundaryPatch(int i,int j,int level);

    // Update all data at a given level
    void refreshMedialAllPatches(int level);
    void refreshBoundaryAllPatches(int level);

    // Get the grid for this cache
    SplineGridDefinition *getGrid() {
        return grid;
    }

    // Get the maximum level used
    int getMaxLevel() {
        return maxLevel;
    }

    // Get the timestamp at a patch at a level
    long getPatchTimeStamp(int iPatch,int jPatch,int level,Types type) {
        return patch(iPatch,jPatch)->getTimeStamp(level,type);
    }

    // Get the timestamp for the whole level
    long getTimeStamp(int level,Types type) {
        return tsPatchUB[type][level];
    }

    // The patches
    Array2D<PatchDataCache*> patch;

private:
    // The maximum resolution level
    int maxLevel;

    // The size of the spline in patches
    int p[2];

    // The spline that is our data source
    MSpline *spline;

    // A grid of points (with weights)
    SplineGridDefinition *grid;

    // Overall time stamps
    long *tsPatchUB[NUM_TYPES],*tsPatchLB[NUM_TYPES];
};

#endif
