#ifndef _IMAGING_CPP_
#define _IMAGING_CPP_

#include <smlmath.h>
#include "array2d.h"
#include "align.h"

#include <iostream>

//#define FFTW_ENABLE_FLOAT 1
//#include <fftw.h>

/**
 * Data cube.  An arrangement of data into slices, rows and elements
 */ 
template<class T> class DataCube
{
private:    
  // Raw voxel data in z, y, x order
  T* voxels;

  // Indexing information 
  int *offSlice,*offRow;

  // Dimensions
  int dim[3];

  // Sizes
  int szCube,szSlice,szRow;

  // Allocation operator
  void alloc(int x,int y,int z) {
    // Initialize the dimensions
    dim[0] = x;
    dim[1] = y;
    dim[2] = z;

    // Initialize the sizes
    szRow = x;
    szSlice = y * szRow;
    szCube = z * szSlice;

    // Allocate the array (use 16 bit allocation)
    // voxels = (T*) _mm_malloc(szCube * sizeof(T),16);
    voxels = new T[szCube];
    offSlice = new int[z];
    offRow = new int[y];

    // Build the offset arrays
    offSlice[0] = 0;
    for (int i=1;i<z;i++)
      {
      offSlice[i] = offSlice[i-1] + szSlice;
      }

    offRow[0] = 0;
    for (int j=1;j<y;j++)
      {
      offRow[j] = offRow[j-1] + szRow;
      }
  }

  void dealloc() {
    delete[] offSlice;
    delete[] offRow;
    // _mm_free(voxels);
    delete[] voxels;
  }

public:

  // Access to slices and rows
  T *root() {
    return voxels;
  }

  T *slice(int z) {
    return voxels + offSlice[z];
  }

  T *row(int y,int z) {
    return voxels + offSlice[z] + offRow[y];
  }

  T *voxel(int x,int y,int z) {
    return voxels + offSlice[z] + offRow[y] + x;
  }

  T& operator()(int x,int y,int z) {
    return voxels[offSlice[z] + offRow[y] + x];
  }

  T& operator()(int offset) {
    return voxels[offset];
  }

  /** This method rapidly loads eight adjacent voxels at indices
   x,y,z,x+1,y+1,z+1 */
  void getEightVoxelCube(int x, int y, int z, float *voxels)
    {
    // Traverse the voxels in a pattern that minimizes the number of 
    // pointer addition operations
    T* root = voxel(x,y,z);                           // x:0  y:0  z:0
    voxels[0] = (float) *root; root += 1;             // x:1  y:0  z:0
    voxels[3] = (float) *root; root += szRow;         // x:1  y:1  z:0
    voxels[6] = (float) *root; root -= 1;             // x:0  y:1  z:0
    voxels[2] = (float) *root; root += szSlice;       // x:0  y:1  z:1

    voxels[4] = (float) *root; root += 1;             // x:1  y:1  z:1
    voxels[7] = (float) *root; root -= szRow;         // x:1  y:0  z:1
    voxels[5] = (float) *root; root -= 1;             // x:0  y:0  z:1
    voxels[1] = (float) *root;   
    }


  int size(int d) {
    return dim[d];
  }

  const int *size() {
    return dim;
  }

  int sizeOfCube() {
    return szCube;
  }

  int sizeOfSlice() {
    return szSlice;
  }

  int sizeOfRow() {
    return szRow;
  }

  void resize(int x,int y,int z) {
    if (voxels)
      dealloc();
    alloc(x,y,z);
  }

  void resize(int x,int y,int z,const T &value) {
    if (voxels)
      dealloc();
    alloc(x,y,z);

    for (int i=0;i<szCube;i++)
      voxels[i] = value;
  }

  // Simple Constructor
  DataCube(int x,int y,int z,const T &value) {
    alloc(x,y,z);
    for (int i=0;i<szCube;i++)
      voxels[i] = value;
  }

  // Simple Constructor
  DataCube(int x,int y,int z) {
    alloc(x,y,z);
  }

  // Dumb Constructor
  DataCube() {
    voxels = NULL;
    szCube = szSlice = szRow = 0;
    dim[0] = dim[1] = dim[2] = 0;
  }

  virtual ~DataCube() {
    if (voxels)
      dealloc();
  }
  };

/**
 * A generic type independant parent for image objects
 */
class AbstractImage3D
  {
public:
  AbstractImage3D() {
    mmData = (float*)_aligned_malloc(sizeof(float)*20,16);
  }

  virtual ~AbstractImage3D() {
    _aligned_free(mmData);
  }

  // Generic method to compute a slice
  virtual unsigned char *getSlice(int dim,int slice,int &res,float &trimS,float &trimT) = 0;

  // Are we loaded
  virtual bool isLoaded() = 0;

  // Image dimensions
  virtual int size(int d) = 0;

  // Get a pixel, no bounds checking
  virtual short getVoxelNBC(int x,int y,int z) = 0;
  virtual float getVoxelNBCFloat(int x,int y,int z) = 0;

  // Get a pixel value at a floating point position, fast, uses SIMD
  // The parameter vector should be in space coordinates
  virtual float interpolateVoxel(float x,float y,float z) = 0;
  virtual void interpolateVoxelGradient(float x,float y,float z, float *G) = 0;

  // Scale the image
  void setVoxelSize(float sx,float sy,float sz);

  // Blur the image
  virtual void gaussianBlur(float *scale) = 0;

  // Get the voxel size
  SMLVec3f getVoxelSize() {
    return vd;
  }

  // Image transforms into space
  SMLMatrix4f SI,IS;

protected:
  // Create a default space-image transform
  void makeDefaultTransform(SMLVec3f vd);

  // A method to perform trilinear interpolation
  void doInterpolation(const __m128 &r0,const __m128 &i0,const __m128 &i1, float *rtn);

  // Image dimensions (x,y,z)
  int dim[3];

  // Voxel size
  SMLVec3f vd;

  // Scale and translation vectors, aligned for speed
  float *mmData;
  };

/**
 * A generic cube of cubes used to access an image efficiently.  This can be used for images,
 * distance transforms,or any other 3D structure.
 */
template<class T> class ImageCube : public AbstractImage3D
  {
public:
  enum
    {
    LDCUBE = 4
    };

  // Generate the distance transform from an image file
  ImageCube() : AbstractImage3D()
  {
    // Resolution of each cube
    res = 1 << LDCUBE;
  }


  // Load from my own cube format file.  Catch exceptions here!
  virtual void loadFromFile(const char *file);
  virtual void saveToFile(const char *file);

  // Load from GIPL file
  virtual void loadFromITKReadableFile(const char *giplFile,std::ostream &out = cout);

  // Voxel interpolation method for cubed images
  virtual float interpolateVoxel(float x,float y,float z);
  virtual void interpolateVoxelGradient(float x,float y,float z, float *G);

  // Generic method to compute a slice
  unsigned char *getSlice(int dim,int slice,int &res,float &trimS,float &trimT);

  // Get a pixel, no bounds checking
  short getVoxelNBC(int x,int y,int z) {
    return (short) cube(x >> LDCUBE,y >> LDCUBE, z >> LDCUBE)(x % res, y % res, z % res);
  }

  float getVoxelNBCFloat(int x,int y,int z) {
    return cube(x >> LDCUBE,y >> LDCUBE, z >> LDCUBE)(x % res, y % res, z % res);
  }

  // Are we loaded
  bool isLoaded() {
    return cube.sizeOfCube() > 0;
  }

  // Image dimensions
  int size(int d) {
    return dim[d];
  }

  // Blur the image
  virtual void gaussianBlur(float *scale);

  // Crop to ROI
  virtual void autoCrop(float padding = 0.2,int resampling = 1,unsigned short background=0);

  // This very useful method takes an ZYX order array of data and places the data into cubes.
  // They array of data must be of the right size
  void populateCubes(const T* dataZYX);
  void unpackCubes(T* dataZYX);

  // A data cube of image cubes
  DataCube< DataCube<T> > cube;

protected:
  // Cached image slice for texture mapping or similar tasks
  Array2D<unsigned char> slice;

  // Resolution of the cube
  int res;

  // Number of cubes in directions (x,y,z)
  int cd[3];

  // Allocate the slice for texture mapping
  void allocateSlice();

  // Utility method
  // int roundup(int x);

  // Method to load from a padded cube
  virtual void loadFromPaddedCube(const T *pixels, const int *size);

  // Internal method to get four voxels around an xyz triple
  bool getEightVoxels(float x,float y,float z,__m128 &r0,__m128 &r1,__m128 &r2);

  // Method used to convert T to pixel for slice computation
  virtual unsigned char genSlicePixel(const T& pixel) {
    return (unsigned char) pixel;
  };

  // Image background : value for pixels outside of the image
  virtual float handleOutsideVoxel(float x,float y,float z);
  };

class BinaryImage;
class GrayImage : public ImageCube<unsigned short>
{
public:
  // Threshold the image
  void threshold(unsigned short intensity);

  // Override the parent's load method
  void loadFromFile(const char *file);

  // Load from a binary image
  void loadFromBinary(BinaryImage &image);

protected:
  // The shift value for converting to character data
  int charShift;

  // Compute the shift amount
  void computeDisplayShift();

  // Override the parent's method to compute the shift value
  void loadFromPaddedCube(const unsigned short *pixels, const int *size);

  unsigned char genSlicePixel(const unsigned short& pixel) 
    { return(unsigned char) (pixel >> charShift); }
};

class BinaryImage : public ImageCube<unsigned char>
{
public:

  // Load from a boundary which is a triangle mesh
  void loadFromBYUFile(const char *byuFile,int maxRes);

  // Crop the image around a segmented section
  //void autoCrop(float padding = 0.2,int resampling = 1,short background=0);

  // Load an image from a grayscale image
  void loadFromGray(GrayImage &image);

  // Threshold the image
  void threshold(unsigned char background);

protected:
  unsigned char genSlicePixel(const unsigned char& pixel) {
    return pixel;
  }
};

class DistanceTransform : public ImageCube<float>
{
public:

  // Generate the distance transform from a binary image file
  void loadFromBinaryImage(BinaryImage &bin,int dtMetric = 0);

  // Slice pixel computation
  unsigned char genSlicePixel(const float& pixel);

  // Mask inside or outside of the transform.  Sets pixels inside or outside to huge values, 
  // so forcing the object to fit the data from inside or outside (shrinkwrapping)
  void maskInside(float addVal = 1000.0f,float scaleVal = 0.0f);
  void maskOutside(float addVal = 1000.0f,float scaleVal = 0.0f);

protected:
  float handleOutsideVoxel(float x,float y,float z);

private:
  // This method does the transform
  // void applyToBinaryData(DataCube<unsigned char> &image);
  void applyToFilledBinaryData(DataCube<unsigned char> &image);

  float vlen(int x,int y,int z);

};

// Convert a distance map value to an RGB color
SMLVec3f convertDistanceToRGB(float pixel);

// Convert a t-score to RGB color 
void convertTScoreToRGB(float pixel, float &r, float &g, float &b);

#include "imaging.txx"

#endif 
