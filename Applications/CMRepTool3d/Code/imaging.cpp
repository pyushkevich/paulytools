#include "Danielsson.h"
#include <mmintrin.h>
#include "imaging.h"
#include "DrawTriangles.h"
#include "DrawFillInside.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>

#include "imaging.txx"

using namespace std;

template class DataCube<unsigned char>;
// template class DataCube<unsigned short>;
// template class DataCube<char>;
// template class DataCube<short>;
template class DataCube<float>;

template class ImageCube<unsigned char>;
// template class ImageCube<unsigned short>;
// template class ImageCube<char>;
// template class ImageCube<short>;
template class ImageCube<float>;

#ifndef M_PI
extern const double M_PI;
#endif

// If intel is not available, define the f32Vec4

void AbstractImage3D::makeDefaultTransform(SMLVec3f vd) 
{  
  // Save the transform
  this->vd = vd;

  // Scaled image dimensions, in voxel coordinates
  SMLVec3f dvs(vd[0] * dim[0], vd[1] * dim[1], vd[2] * dim[2]);

  // Direction of the greatest extent of the image
  int imax = dvs[0] < dvs[1] ? (dvs[1] < dvs[2] ? 2 : 1) : (dvs[0] < dvs[2] ? 2 : 0);
  float dmax = dvs[imax];

  // The 'left' border of the image, when placed in the unit box
  SI.set_identity();
  for (int i=0;i<3;i++)
    {
    mmData[i] = dmax / vd[i];
    mmData[i+4] = (dim[i] - mmData[i]) / 2.0f;
    mmData[i+8] = dim[i] - 1.0001f;
    mmData[i+12] = 1.0 / mmData[i];

    SI.put(i,i,mmData[i]);
    SI.put(i,3,mmData[i+4]);
    }   
  mmData[3]=mmData[7]=mmData[11]=mmData[15]=0.0f;

  // Invert the matrix
  IS = vnl_matrix_inverse<float>(SI);
/*      
    // Find the maximum image dimension
    int max = dim[0] > dim[1] ? (dim[0] > dim[2] ? dim[0] : dim[2]) : (dim[1] > dim[2] ? dim[1] : dim[2]);

    // Initialize matrices
    SI.Identity();

    // Scale by max dimension
    for(int i=0;i<3;i++) {
        SI.Set(i,i,max);
        SI.Set(i,3,(dim[i]-max)>>1);
    }

    // Save the fast values for intel image access
    mmData[0]  = mmData[1] = mmData[2] = max;
    mmData[3]  = 0.0f;
    mmData[4]  = SI.Get(0,3);
    mmData[5]  = SI.Get(1,3);
    mmData[6]  = SI.Get(2,3);
    mmData[7]  = 0.0f;
    mmData[8]  = dim[0]-1.0f;
    mmData[9]  = dim[1]-1.0f;
    mmData[10] = dim[2]-1.0f;
    mmData[11] = 0.0f;
    mmData[12]  = mmData[13] = mmData[14] = mmData[15] = 1.0f / max;

    // Invert the matrix
    IS.Invert(SI);
*/
}

/**
 Scale the image by an amount
 */
void AbstractImage3D::setVoxelSize(float sx,float sy,float sz) {
  makeDefaultTransform(SMLVec3f(sx,sy,sz));
}

/*
C:\Program Files\Intel\Compiler50\ia32\include\xmm_func.h
typedef struct {
    float m128_f32[4];
} __m128;
*/

void AbstractImage3D::doInterpolation(const __m128 &r0,const __m128 &in0,const __m128 &in1, float *rtn) 
{
  // Registers we are going to use
  __m128 r1,r2,r3,r4,r5,r6,r8,i0,i1;                      //  3       2       1       0

  // Perform the interpolation
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

  // Scale the intensities by multipliers
  i0 = _mm_mul_ps(in0,r8);                                 
  i1 = _mm_mul_ps(in1,r3);                                 

  // Now add up all elements
  i0 = _mm_add_ps(i0,i1);
  i1 = _mm_shuffle_ps(i0,i0,0xb1);
  i0 = _mm_add_ps(i0,i1);
  i1 = _mm_shuffle_ps(i0,i0,0x4e);
  i0 = _mm_add_ps(i0,i1);

  // The result is in r1
  _mm_store_ss(rtn, i0);
}

inline ostream &operator << (ostream &sout, __m128 r)
{
  ALIGN_PRE float f[4] ALIGN_POST;
  _mm_storeu_ps(f, r);

  sout << "[" 
    << f[0] << "," 
    << f[1] << "," 
    << f[2] << "," 
    << f[3] << "]";

  return sout;
} 

void GrayImage::computeDisplayShift() { 
  // OR all voxels in the image
  unsigned short d = 0;
  for (int c=0;c<cube.sizeOfCube();c++)
    {
    DataCube<unsigned short> &dc = cube(c);
    for (int i=0;i<dc.sizeOfCube();i++)
      {
      d |= dc(i);
      }
    }

  // Right shift until d = 0
  charShift = 0;
  d >>= 8;
  while (d > 0)
    {
    charShift++;
    d >>= 1;
    }
}

/*
void DistanceTransform::applyToBinaryData(DataCube<unsigned char> &image) {

    // Create black/white images
    for(int z=0;z<image.size(2);z++) {
        for(int y=0;y<image.size(1);y++) {
            for(int x=0;x<image.size(0);x++) {
                image(x,y,z) = image(x,y,z) ? 1 : 0;
            }
        }
    }

    for(int d=0;d<3;d++) {
        // Compute cube dimensions
        cd[d] = image.size(d) >> LDCUBE;
    }

    // Distance map values
    short *od[3];

    // Perform the distance transform
    edtSetVoxelSize(vd.x,vd.y,vd.z);
    edt3ddan((char *)image.root(),image.size(0),image.size(1),image.size(2),0,od+0,od+1,od+2);

    // Compute the distance, store into od
    for(int i=0;i<image.sizeOfCube();i++) {
        od[0][i] = vlen(od[0][i],od[1][i],od[2][i]);
    }

    // Free unused data
    free(od[1]);
    free(od[2]);

    // Create image cubes from the transform
    cube.resize(cd[0],cd[1],cd[2]);
    populateCubes(od[0]);

    // Delete the data
    free(od[0]);
}
*/
/*
template class<T>
void ImageCube<T>::blur(float sigma, int kernelSize)
{
    double * kernel;
    T * data;
    T * newData;
    T * dataPtr;
    T * newDataPtr;
    int minKernelIndex, maxKernelIndex;
    int xDim, yDim, zDim;
    int imageSize;
    int i, j, xIndex, yIndex, zIndex;
    int stepSize;
    double val;
    double factor;

    cout << "Starting Gauss blur" << endl;

    // Make sure we have an odd-sized kernel
    if(kernelSize % 2 == 0)
        kernelSize++;

    maxKernelIndex = kernelSize / 2;
    minKernelIndex = -maxKernelIndex;

    kernel = new double[kernelSize];
    factor = 1.0 / (sqrt(R_TWO_PI) * sigma);
    for(i = minKernelIndex; i <= maxKernelIndex; i++)
        kernel[i - minKernelIndex] = factor * exp((-0.5 / (sigma * sigma)) * (double)(i * i));

    data = image.getVoxels();

    xDim = image.getXDim();
    yDim = image.getYDim();
    zDim = image.getZDim();

    imageSize = xDim * yDim * zDim;

    newData = new GreyValue[imageSize];
    if(newData == NULL)
    {
        cerr << "Cannot allocate space for image" << endl;
        return;
    }

    // Convolve in X
    dataPtr = data;
    newDataPtr = newData;
    for(i = 0; i < imageSize; i++)
    {
        xIndex = i % xDim;

        val = 0;
        for(j = minKernelIndex; j <= maxKernelIndex; j++)
        {
            if(xIndex + j < 0)
                val += (double)(*(dataPtr - xIndex)) * kernel[j - minKernelIndex];
            else if(xIndex + j >= xDim)
                val += (double)(*(dataPtr + xDim - xIndex - 1)) * kernel[j - minKernelIndex];
            else
                val += (double)(*(dataPtr + j)) * kernel[j - minKernelIndex];
        }

        (*newDataPtr) = (GreyValue) val;
        dataPtr++;
        newDataPtr++;
    }

    // Convolve in Y
    dataPtr = newData;
    newDataPtr = data;
    for(i = 0; i < imageSize; i++)
    {
        xIndex = i % xDim;
        yIndex = ((i - xIndex) / xDim) % yDim;

        val = 0;
        for(j = minKernelIndex; j <= maxKernelIndex; j++)
        {
            if(yIndex + j < 0)
                val += (double)(*(dataPtr - yIndex * xDim)) * kernel[j - minKernelIndex];
            else if(yIndex + j >= yDim)
                val += (double)(*(dataPtr + (yDim - yIndex - 1) * xDim)) * kernel[j - minKernelIndex];
            else
                val += (double)(*(dataPtr + j * xDim)) * kernel[j - minKernelIndex];
        }

        (*newDataPtr) = (GreyValue) val;
        dataPtr++;
        newDataPtr++;
    }

    // Convolve in Z
    dataPtr = data;
    newDataPtr = newData;
    stepSize = xDim * yDim;
    for(i = 0; i < imageSize; i++)
    {
        xIndex = i % xDim;
        yIndex = ((i - xIndex) / xDim) % yDim;
        zIndex = (i - xIndex - xDim * yIndex) / stepSize;

        val = 0;
        for(j = minKernelIndex; j <= maxKernelIndex; j++)
        {
            if(zIndex + j < 0)
                val += (double)(*(dataPtr - zIndex * stepSize)) * kernel[j - minKernelIndex];
            else if(zIndex + j >= zDim)
                val += (double)(*(dataPtr + (zDim - zIndex - 1) * stepSize)) * kernel[j - minKernelIndex];
            else
                val += (double)(*(dataPtr + j * stepSize)) * kernel[j - minKernelIndex];
        }

        (*newDataPtr) = (GreyValue) val;
        dataPtr++;
        newDataPtr++;
    }

    delete [] kernel;

    image.setVoxels(newData, xDim, yDim, zDim);
    cout << "done" << endl;
}
*/


void DistanceTransform::applyToFilledBinaryData(DataCube<unsigned char> &image) {

  // Create black/white images
  DataCube<char> negative(image.size(0),image.size(1),image.size(2),127);
  for (int z=0;z<image.size(2);z++)
    {
    for (int y=0;y<image.size(1);y++)
      {
      for (int x=0;x<image.size(0);x++)
        {
        image(x,y,z) = image(x,y,z) ? 127 : 0;
        negative(x,y,z) = 127 - image(x,y,z);
        }
      }
    }

  for (int d=0;d<3;d++)
    {
    // Compute cube dimensions
    cd[d] = image.size(d) >> LDCUBE;
    }

  // Distance map values
  short *od[3],*id[3];

  // Perform the distance transform
  edtSetVoxelSize(vd[0],vd[1],vd[2]);
  edt3ddan((char *)image.root(),image.size(0),image.size(1),image.size(2),0,od+0,od+1,od+2);
  edt3ddan(negative.slice(0),image.size(0),image.size(1),image.size(2),0,id+0,id+1,id+2);

  // Deallocate the negative image
  negative.resize(2,2,2);

  // Allocate a float image
  float *d2 = new float[image.sizeOfCube()];

  // Compute the distance, store into od, id arrays   
  for (int i=0;i<image.sizeOfCube();i++)
    {
    float dist1 = vlen(od[0][i],od[1][i],od[2][i]);
    float dist2 = vlen(id[0][i],id[1][i],id[2][i]);

    // We are computing signed distance: the distance inside is negative.  It's important to
    // make sure we fabs() this distance for important measurements
    d2[i] = dist1 - dist2;
    // od[0][i] = dist1 > dist2 ? dist1 : dist2;
    }

  // Free the unneeded arrays 
  free(od[0]);
  free(od[1]);
  free(od[2]);
  free(id[0]);
  free(id[1]);
  free(id[2]);

  // Create image cubes from the transform
  cube.resize(cd[0],cd[1],cd[2]);
  populateCubes(d2);

  // Remove the data
  delete d2;
}

void DistanceTransform::loadFromBinaryImage(BinaryImage &bin,int metric) {

  // Resize to match the binary image
  dim[0] = bin.size(0);
  dim[1] = bin.size(1);
  dim[2] = bin.size(2);

  // Create an unpacked image
  DataCube<unsigned char> image(dim[0],dim[1],dim[2],0);
  bin.unpackCubes(image.root());

  // Generate a transform, maybe later load it from the image
  makeDefaultTransform(bin.getVoxelSize());

  // Apply the distance transform
  applyToFilledBinaryData(image);

  // Create a slice from the image
  allocateSlice();
}

float DistanceTransform::handleOutsideVoxel(float x,float y,float z) {
  return(float)0xffff;
}

unsigned char DistanceTransform::genSlicePixel(const float& pixel) {
  // This returns the pixel color index
  const float scale1 = 256.0f / 3.1416f;
  const float scale2 = 1.0f / 12.0f;

  // return (pixel > 255) ? 255 : (unsigned char)pixel;
  unsigned char c = (unsigned char)(128 + (atan(pixel * scale2)*scale1));
  return c;
}

void DistanceTransform::maskInside(float addVal,float scaleVal) {
  for (int z=0;z<cd[2];z++)
    {
    for (int y=0;y<cd[1];y++)
      {
      for (int x=0;x<cd[0];x++)
        {
        DataCube<float> &C = this->cube(x,y,z);
        for (int k=0;k<=res;k++)
          {
          for (int j=0;j<=res;j++)
            {
            for (int i=0;i<=res;i++)
              {
              float &f = C(i,j,k);
              if (f < 0)
                f = addVal + scaleVal * f;
              }
            }
          }
        }
      }
    }
}

void DistanceTransform::maskOutside(float addVal,float scaleVal) {
  for (int z=0;z<cd[2];z++)
    {
    for (int y=0;y<cd[1];y++)
      {
      for (int x=0;x<cd[0];x++)
        {
        DataCube<float> &C = this->cube(x,y,z);
        for (int k=0;k<=res;k++)
          {
          for (int j=0;j<=res;j++)
            {
            for (int i=0;i<=res;i++)
              {
              float &f = C(i,j,k);
              if (f > 0)
                f = addVal + scaleVal * f;
              }
            }
          }
        }
      }
    }
}

void convertTScoreToRGB(float pixel, float &r, float &g, float &b)
{
  // Define the color map 
  const int nColors = 7;
  const float cmap[nColors][3] = {
    {0.0f,1.0f,0.0f},
    {0.0f,1.0f,1.0f},
    {0.0f,0.0f,1.0f},
    {0.0f,0.0f,0.0f},
    {1.0f,0.0f,0.0f},
    {1.0f,1.0f,0.0f},
    {1.0f,1.0f,1.0f}};

  // Scale the value
  float val = 3.0f + 0.5f * pixel;

  // Find the range for scaling
  int i0 = (int) val; 
  int i1 = i0 + 1; 
  
  if(i0 < 0) 
    i1 = i0 = 0;
  if(i1 >= nColors) 
    i0 = i1 = nColors-1;

  // Interpolate
  float rem = val - i0;
  r = rem * cmap[i1][0] + (1 - rem) * cmap[i0][0];
  g = rem * cmap[i1][1] + (1 - rem) * cmap[i0][1];
  b = rem * cmap[i1][2] + (1 - rem) * cmap[i0][2];
}


SMLVec3f convertDistanceToRGB(float pixel) 
{
  // This returns the pixel color index
  const float scale1 = 1.0f / 3.1416f;
  const float scale2 = 1.0f / 12.0f;
  float i = atan(pixel * scale2)*scale1;  
  float i2 = i+i;

  if (i < 0)
    {
    return SMLVec3f(1.0f+i2,0.5f+i,0.5f+i);
    }
  else if (i > 0)
    {
    return SMLVec3f(0.5f-i,0.5f-i,1.0f-i2);
    }
  else
    {
    return SMLVec3f(1.0f,1.0f,1.0f);
    }
}

inline float DistanceTransform::vlen(int x,int y,int z)  {
  float x1 = vd[0] * x;
  float y1 = vd[1] * y;
  float z1 = vd[2] * z;
  return x1*x1 + y1*y1 + z1*z1;
}

void BinaryImage::loadFromBYUFile(const char *byuFile,int vmaxRes) {
  int i,index;
  int numVertices,numTriangles,numEdges;

  FILE *fp = fopen(byuFile, "rb");
  if (fp == NULL)
    {
    throw "Can not open file for reading";
    }

  // Read the number of triangles
  fscanf(fp, "1 %d %d %d\n", &numVertices, &numTriangles, &numEdges);
  fscanf(fp, "1 %d", &numEdges);

  // Allocate the vertex list and triangle list
  double *vertexList = new double[numVertices * 3];
  double **vertexTable = new double *[numTriangles * 3];

  // XYZ bounds
  SMLVec3f vmax,vmin;

  // Read in the vertices
  index = 0;
  for (i = 0; i < numVertices; i++)
    {
    float x,y,z;
    fscanf(fp, "%f %f %f\n", &x,&y,&z);

    if (i==0)
      {
      vmin[0] = vmax[0] = x;
      vmin[1] = vmax[1] = y;
      vmin[2] = vmax[2] = z;
      }
    else
      {
      vmax[0] = x > vmax[0] ? x : vmax[0];
      vmax[1] = y > vmax[1] ? y : vmax[1];
      vmax[2] = z > vmax[2] ? z : vmax[2];
      vmin[0] = x < vmin[0] ? x : vmin[0];
      vmin[1] = y < vmin[1] ? y : vmin[1];
      vmin[2] = z < vmin[2] ? z : vmin[2];
      }

    vertexList[index++] = x;
    vertexList[index++] = y;
    vertexList[index++] = z;
    }

  // Read the triangles
  index = 0;
  for (i = 0; i < numTriangles; i++)
    {
    int a,b,c;
    fscanf(fp, "%d %d -%d\n", &a, &b, &c);
    vertexTable[index++] = &(vertexList[3*(a-1)]);
    vertexTable[index++] = &(vertexList[3*(b-1)]);
    vertexTable[index++] = &(vertexList[3*(c-1)]);
    }

  fclose(fp);

  // Compute the triangle to image transform
  SMLVec3f C = (vmax + vmin) * 0.5f;
  SMLVec3f D = (vmax - vmin) * 1.5f;

  // The largest dimension
  int mdix = D[0] > D[1] ? (D[0] > D[2] ? 0 : 2) : (D[1] > D[2] ? 1 : 2);   
  float md = D[mdix];

  // Compute the dimensions of the new image
  for (int d=0;d<3;d++)
    {
    dim[d] = (mdix==d) ? vmaxRes 
             : (1 << (int)ceil(log(vmaxRes * D[d] / md) / log(2.0)));
    cd[d] = dim[d] >> LDCUBE;
    }

  // Update the coordinates
  for (i = 0; i < numVertices; i++)
    {
    vertexList[i*3] = (vertexList[i*3] - C[0]) * vmaxRes / md + dim[0] * 0.5;
    vertexList[i*3+1] = (vertexList[i*3+1] - C[1]) * vmaxRes / md + dim[1] * 0.5;
    vertexList[i*3+2] = (vertexList[i*3+2] - C[2]) * vmaxRes / md + dim[2] * 0.5;
    }

  // Create a binary image
  DataCube<unsigned char> image(dim[0],dim[1],dim[2],0);

  drawBinaryTrianglesSheetFilled(image.slice(0), dim, vertexTable, numTriangles);
  draw_fill_inside_image(image.slice(0), dim);

  delete vertexList;
  delete vertexTable;

  // Populate the data cubes
  cube.resize(cd[0],cd[1],cd[2]);
  populateCubes(image.root());    

  // Generate a transform, maybe later load it from the image
  makeDefaultTransform(SMLVec3f(1.0f,1.0f,1.0f));

  // Create a slice from the image
  allocateSlice();
}

void GrayImage::loadFromFile(const char *file) {
  // Call parent
  ImageCube<unsigned short>::loadFromFile(file);
  computeDisplayShift();
}

void GrayImage::loadFromBinary(BinaryImage &image)
{
  DataCube<unsigned char> cube(image.size(0),image.size(1),image.size(2));
  image.unpackCubes(cube.root());

  DataCube<unsigned short> c2(image.size(0),image.size(1),image.size(2));
  for (int i=0;i<cube.sizeOfCube();i++)
    {
    c2(i) = (cube(i) == 0) ? 0 : 0x7fff;
    }

  loadFromPaddedCube(c2.root(), c2.size());
  makeDefaultTransform(image.getVoxelSize());
}

void GrayImage::loadFromPaddedCube(const unsigned short *pixels, const int *size) {
  // Call parent
  ImageCube<unsigned short>::loadFromPaddedCube(pixels,size);
  computeDisplayShift();
}

void BinaryImage::loadFromGray(GrayImage &image) {
  DataCube<unsigned short> cube(image.size(0),image.size(1),image.size(2));
  image.unpackCubes(cube.root());

  DataCube<unsigned char> c2(image.size(0),image.size(1),image.size(2));
  for (int i=0;i<cube.sizeOfCube();i++)
    {
    c2(i) = (cube(i) == 0) ? 0 : 255;
    }

  loadFromPaddedCube(c2.root(), c2.size());
  makeDefaultTransform(image.getVoxelSize());
}

void BinaryImage::threshold(unsigned char background) {
  for (int z=0;z<cd[2];z++)
    {
    for (int y=0;y<cd[1];y++)
      {
      for (int x=0;x<cd[0];x++)
        {
        DataCube<unsigned char> &C = this->cube(x,y,z);
        for (int k=0;k<res;k++)
          {
          for (int j=0;j<res;j++)
            {
            for (int i=0;i<res;i++)
              {
              unsigned char &c = C(i,j,k);
              c = (c==background) ? 0 : 255;
              }
            }
          }
        }
      }
    }
}

void GrayImage::threshold(unsigned short intensity) {
  for (int z=0;z<cd[2];z++)
    {
    for (int y=0;y<cd[1];y++)
      {
      for (int x=0;x<cd[0];x++)
        {
        DataCube<unsigned short> &C = this->cube(x,y,z);
        for (int k=0;k<res;k++)
          {
          for (int j=0;j<res;j++)
            {
            for (int i=0;i<res;i++)
              {
              unsigned short &c = C(i,j,k);
              c = (c < intensity) ? 0 : 32767;
              }
            }
          }
        }
      }
    }
}


/*
template<class T> 
void ImageCube<T>::gaussianBlur(float *sigma) 
{
    // Convert packed data to straight data
    DataCube<T> unpack(dim[0],dim[1],dim[2]);
    unpackCubes(unpack.root());

    // Allocate a fftw_complex array
    DataCube<fftw_complex> A(image.size(0),image.size(1),image.size(2));
    for(int i=0;i<unpack.sizeOfCube();i++) {
        A(i).re = unpack(i);
        A(i).im = 0;
    }

    // Create a plan for FFT
    fftwnd_plan plan = fftw3d_create_plan(A.size(0),A.size(1),A.size(2),FFTW_FORWARD,FFTW_IN_PLACE);

    // Perform fft on the complex data
    fftwnd_one(plan, A.root(), NULL);

    // Destroy the plan
    fftwnd_destroy_plan(plan);

    // Multiply by a Gaussian

    // Create a plan for IFFT
    //plan = fftw3d_create_plan(A.size(0),A.size(1),A.size(2),FFTW_INVERSE,FFTW_IN_PLACE);

    // Perform fft on the complex data
    //fftwnd_one(plan, A.root(), NULL);

    // Destroy the plan
    //fftwnd_destroy_plan(plan);

    // Convert the image back
    for(int i=0;i<unpack.sizeOfCube();i++) {
        unpack(i) = A(i).re;
    }

    populateCubes(unpack.root());
}
*/

