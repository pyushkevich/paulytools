// #include <RAWImageFile.h>
#include "Danielsson.h"
#include <mmintrin.h>
#include <ivec.h>
#include <fvec.h>
#include "imaging.h"
#include "DrawTriangles.h"
#include "DrawFillInside.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>

#include <itkImage.h>
#include <itkCastImageFilter.h>
#include <itkImageFileReader.h>
#include <itkDiscreteGaussianImageFilter.h>

using namespace std;
using namespace itk;

ImageCube<float> test1;
ImageCube<unsigned short> test3;
ImageCube<unsigned char> test2;

extern const double M_PI;

#pragma auto_inline(off)

/*
template <class T>
void DataCube<T>::alloc(int x,int y,int z) {
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
    for(int i=1;i<z;i++) {
        offSlice[i] = offSlice[i-1] + szSlice;
    }

    offRow[0] = 0;
    for(int j=1;j<y;j++) {
        offRow[j] = offRow[j-1] + szRow;
    }
}*/
/*
template <class T>
void DataCube<T>::dealloc(void) {
    delete[] offSlice;
    delete[] offRow;
    // _mm_free(voxels);
    delete[] voxels;
}*/

/**********************************************************
 IMAGE CUBE
 **********************************************************/

template <class T>
void ImageCube<T>::autoCrop(float padding,int resampling,unsigned short background) {
  // Extent vectors
  int mx[3] = {0,0,0},mn[3] = {0,0,0};
  bool firstTime = true;

  // Find the extent of the non-zero voxels
  for (int z=0;z<dim[2];z++)
    {
    for (int y=0;y<dim[1];y++)
      {
      for (int x=0;x<dim[0];x++)
        {
        if ((unsigned short) getVoxelNBC(x,y,z) > background)
          {
          if (firstTime)
            {
            mx[0] = x+1;
            mx[1] = y+1;
            mx[2] = z+1;
            mn[0] = x;
            mn[1] = y;
            mn[2] = z;
            firstTime = false;
            }
          else
            {
            mx[0] = vnl_math_max(x+1,mx[0]);
            mx[1] = vnl_math_max(y+1,mx[1]);
            mx[2] = vnl_math_max(z+1,mx[2]);
            mn[0] = vnl_math_min(x,mn[0]);
            mn[1] = vnl_math_min(y,mn[1]);
            mn[2] = vnl_math_min(z,mn[2]);                       
            }
          }
        }
      }
    }

  // Scale the image
  int span[3], C[3], CN[3], ndim[3], j[3], k[3];
  for (int d=0;d<3;d++)
    {
    span[d] = mx[d] - mn[d];
    C[d] = (mx[d] + mn[d]) >> 1;

    ndim[d] = roundup((int)((span[d]*resampling) * (1.0 + 2.0 * padding)),res);
    CN[d] = ndim[d] >> 1;

    std::cout << "     dim " << d << endl;
    std::cout << "     min " << mn[d] << endl;
    std::cout << "     max " << mx[d] << endl;
    std::cout << "new size " << ndim[d] << endl;
    }

  if (span[0] > 0 && span[1] > 0 && span[2] > 0)
    {
    // Create a temporary image to hold the cropped image
    DataCube<T> paddedCube(ndim[0],ndim[1],ndim[2]);

    for (int t=0;t<paddedCube.sizeOfCube();t++)
      {
      paddedCube.root()[t] = 0;
      }

    // Transform the voxel into the new image
    for (j[2]=mn[2];j[2]<mx[2];j[2]++)
      {
      for (j[1]=mn[1];j[1]<mx[1];j[1]++)
        {
        for (j[0]=mn[0];j[0]<mx[0];j[0]++)
          {
          T vox = (T) getVoxelNBC(j[0],j[1],j[2]);

          int x0[3];
          for (int d=0;d<3;d++)
            {
            x0[d] = CN[d] + resampling * (j[d]-C[d]);
            }

          for (k[2]=0;k[2]<resampling;k[2]++)
            {
            for (k[1]=0;k[1]<resampling;k[1]++)
              {
              for (k[0]=0;k[0]<resampling;k[0]++)
                {
                paddedCube(k[0]+x0[0],k[1]+x0[1],k[2]+x0[2]) = vox;
                }
              }
            }
          }
        }
      }

    // Load the new image
    loadFromPaddedCube(paddedCube.root(), paddedCube.size());

    // Update the transform
    makeDefaultTransform(vd);
    }
}

template <class T>
void ImageCube<T>::gaussianBlur(float *sigma) 
{
  typedef itk::Image<T,3> ImageType;
  typedef itk::Image<float,3> FloatImageType;
  typedef ImageType::RegionType RegionType;
  typedef itk::DiscreteGaussianImageFilter<FloatImageType,FloatImageType> FilterType;
  typedef itk::CastImageFilter<ImageType,FloatImageType> FloatCasterType;
  typedef itk::CastImageFilter<FloatImageType,ImageType> ReverseCasterType;

  // Initialize regions
  ImageType::SizeType size = {{dim[0],dim[1],dim[2]}};
  ImageType::RegionType region(size);
  
  // Create a flat image
  ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->Allocate();

  // Dump the image cube into a flat image
  unpackCubes(image->GetBufferPointer());

  // Create a variance vector
  SMLVec3f variance(
    sigma[0] * sigma[0], sigma[1] * sigma[1], sigma[2] * sigma[2]);

  // Construct and run a pipeline
  FloatCasterType::Pointer floatCaster = FloatCasterType::New();
  floatCaster->SetInput(image);

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput(floatCaster->GetOutput());
  filter->SetVariance(variance.data_block());
  
  ReverseCasterType::Pointer revCaster = ReverseCasterType::New();
  revCaster->SetInput(filter->GetOutput());
  revCaster->Update();

  // Stick the result back 
  loadFromPaddedCube(revCaster->GetOutput()->GetBufferPointer(),dim);
}

template <class T>
unsigned char *ImageCube<T>::getSlice(int d0,int sliceNo,int &outRes,float &trimS,float &trimT) {
  int c[3],x[3];

  // The fixed dimension
  c[d0] = sliceNo >> LDCUBE;
  x[d0] = sliceNo & (res-1);

  // Other dimensions
  int d1 = (d0+1) % 3;
  int d2 = (d0+2) % 3;

  // Slice data pointer
  unsigned char *sld = slice.getData();

  // Offsets at the beginning of each loop
  int off[4];

  off[0] = 0;
  for (c[d2]=0;c[d2]<cd[d2];c[d2]++)
    {
    off[1] = off[0];
    for (c[d1]=0;c[d1]<cd[d1];c[d1]++)
      {
      off[2] = off[1];

      // Cache the selected cube
      DataCube<T> &curCube = cube(c[0],c[1],c[2]);

      // Render the pixels in the cube
      for (x[d2]=0;x[d2]<res;x[d2]++)
        {
        off[3] = off[2];
        for (x[d1]=0;x[d1]<res;x[d1]++)
          {
          sld[off[3]] = genSlicePixel(curCube(x[0],x[1],x[2]));                   
          off[3]++;
          }
        off[2] += slice.width();
        }
      off[1] += res;
      }
    off[0] += res * slice.width();
    }

  // Set outgoing parameters
  outRes = slice.width();

  // Trim values are ratios of sizes to total size
  float rcp = 1.0f / outRes;
  trimS = dim[d1] * rcp;
  trimT = dim[d2] * rcp;

  return sld;
}


template <class T>
void ImageCube<T>::populateCubes(const T* dataZYX) {  
  int off[] = {0,0,0,0,0,0};

  for (int cz=0;cz<cd[2];cz++)
    {
    off[1] = off[0];
    for (int cy=0;cy<cd[1];cy++)
      {
      off[2] = off[1];
      for (int cx=0;cx<cd[0];cx++)
        {
        // The cube is padded by 1 voxel on each side.  For border cubes, the padded area
        // is not initialized
        cube(cx,cy,cz).resize(res+1,res+1,res+1);

        // These values indicate whether we want to fill res or res+1 voxels in each direction
        int wx = cx==cd[0]-1 ? res : res+1;
        int wy = cy==cd[1]-1 ? res : res+1;
        int wz = cz==cd[2]-1 ? res : res+1;

        off[3] = off[2];
        for (int z=0;z<wz;z++)
          {
          off[4] = off[3];
          for (int y=0;y<wy;y++)
            {
            off[5] = off[4];
            for (int x=0;x<wx;x++)
              {
              cube(cx,cy,cz)(x,y,z) = dataZYX[off[5]];
              off[5]++;
              }
            off[4] += dim[0];
            }
          off[3] += dim[0]*dim[1];
          }

        off[2] += res;
        }
      off[1] += dim[0] * res;
      }
    off[0] += dim[0] * dim[1] * res;
    }
}


/*
template <class T, class S>
void unpackCubes(ImageCube<T> &ic,S* dataZYX) { 
    int off[] = {0,0,0,0,0,0};
    for(int cz=0;cz<ic.cd[2];cz++) {
        off[1] = off[0];
        for(int cy=0;cy<ic.cd[1];cy++) {
            off[2] = off[1];
            for(int cx=0;cx<ic.cd[0];cx++) {
                off[3] = off[2];
                for(int z=0;z<ic.res;z++) {
                    off[4] = off[3];
                    for(int y=0;y<ic.res;y++) {
                        off[5] = off[4];
                        for(int x=0;x<ic.res;x++) {
                            dataZYX[off[5]] = ic.cube(cx,cy,cz)(x,y,z);
                            off[5]++;
                        }
                        off[4] += ic.dim[0];
                    }
                    off[3] += ic.dim[0]*ic.dim[1];
                }                           
                off[2] += ic.res;
            }
            off[1] += ic.dim[0] * ic.res;
        }
        off[0] += ic.dim[0] * ic.dim[1] * ic.res;
    }
}
*/

template <class T>
void ImageCube<T>::unpackCubes(T* dataZYX) {    
  int off[] = {0,0,0,0,0,0};
  for (int cz=0;cz<cd[2];cz++)
    {
    off[1] = off[0];
    for (int cy=0;cy<cd[1];cy++)
      {
      off[2] = off[1];
      for (int cx=0;cx<cd[0];cx++)
        {
        off[3] = off[2];
        for (int z=0;z<res;z++)
          {
          off[4] = off[3];
          for (int y=0;y<res;y++)
            {
            off[5] = off[4];
            for (int x=0;x<res;x++)
              {
              dataZYX[off[5]] = cube(cx,cy,cz)(x,y,z);
              off[5]++;
              }
            off[4] += dim[0];
            }
          off[3] += dim[0]*dim[1];
          }                           
        off[2] += res;
        }
      off[1] += dim[0] * res;
      }
    off[0] += dim[0] * dim[1] * res;
    }
}

template<class T>
void ImageCube<T>::allocateSlice() {
  // Compute the slice dimension
  int sd = res;
  while (sd < dim[0] || sd < dim[1] || sd < dim[2])
    sd <<= 1;

  // Allocate the slice
  slice.resize(sd,sd);
}

/*
template<class T>
inline int ImageCube<T>::roundup(int x)  {
    if(x & (res-1))
        x += res - (x & (res-1));
    return x;
}
*/

inline int roundup(int x,int res)  {
  if (x & (res-1))
    x += res - (x & (res-1));
  return x;
}

template<class T>
float ImageCube<T>::handleOutsideVoxel(float x,float y,float z)  {
  return 0.0f;
}

/** 
 * File Format:
 
 First, a 10 byte preheader, of format:
    0..3    file format indicator, 'PYIC'
    4..5    major version, ASCII
    6..7    minor version, ASCII
    8       header length, msb
    9       header length, lsb 
 
 Second a header of length 256*msb+lsb of format 
    key=value;

    required keys are:
        cubeResolution,
        cubeCountX,
        cubeCountY,
        cubeCountZ,
        voxelSize

 Finally, data, written using fwrite, in cubes, in Z/Y/X order
 */
template<class T>
void ImageCube<T>::saveToFile(const char *file) {
  // Create a header string
  ostringstream oss;
  oss << "cubeResolution=" << res << ";";
  oss << "cubeCountX=" << cd[0] << ";";
  oss << "cubeCountY=" << cd[1] << ";";
  oss << "cubeCountZ=" << cd[2] << ";";
  oss << "voxelDimX=" << vd[0] << ";";
  oss << "voxelDimY=" << vd[1] << ";";
  oss << "voxelDimZ=" << vd[2] << ";";
  oss << "voxelSize=" << sizeof(T) << ";";
  oss << "cubePadding=1;";
  string header = oss.str();

  // I use traditional files
  FILE *f = fopen(file,"wb");

  // Generate a fixed length pre-header
  string preheader = "PYIC0103";
  preheader += (unsigned char)(header.length() >> 8);
  preheader += (unsigned char)(header.length() & 0xff);

  // Write header into file
  fwrite(preheader.c_str(),1,10,f);
  fwrite(header.c_str(),1,header.length(),f);

  // Write cube data into the file
  for (int i=0;i<cube.sizeOfCube();i++)
    {
    fwrite(cube.root()[i].root(),sizeof(T),cube.root()[i].sizeOfCube(),f);
    }

  // Close the file
  fclose(f);
}

template<class T>
void ImageCube<T>::loadFromFile(const char *file) {
  // Whether or not the cubes are padded 
  bool paddedCubes = false;

  // I use traditional files
  FILE *f = fopen(file,"rb");
  if (!f)
    throw "File does not exist or is not accessible";

  // Read fixed length header
  char preheader[10];
  fread(preheader,1,10,f);

  // Compare code
  if (strncmp(preheader,"PYIC",4))
    throw "File is not a PY Image Cube";

  // Check header length
  unsigned char l1,l2;
  l1 = preheader[8];
  l2 = preheader[9];
  int hl = (l1 << 8) + l2;

  // Set the vox dimensinos
  vd[0] = vd[1] = vd[2] = 1.0f;

  // Read the header
  char *header = new char[hl+1];
  fread(header,1,hl,f);
  header[hl]=0;

  char *token = strtok(header,"=");
  while (token)
    {
    if (!strcmp(token,"cubeResolution"))
      {
      token = strtok(NULL,";");
      sscanf(token,"%d",&res);
      }
    else if (!strcmp(token,"cubeCountX"))
      {
      token = strtok(NULL,";");
      sscanf(token,"%d",cd+0);
      }
    else if (!strcmp(token,"cubeCountY"))
      {
      token = strtok(NULL,";");
      sscanf(token,"%d",cd+1);
      }
    else if (!strcmp(token,"cubeCountZ"))
      {
      token = strtok(NULL,";");
      sscanf(token,"%d",cd+2);
      }
    else if (!strcmp(token,"voxelDimX"))
      {
      token = strtok(NULL,";");
      sscanf(token,"%g",&vd[0]);
      }
    else if (!strcmp(token,"voxelDimY"))
      {
      token = strtok(NULL,";");
      sscanf(token,"%g",&vd[1]);
      }
    else if (!strcmp(token,"voxelDimZ"))
      {
      token = strtok(NULL,";");
      sscanf(token,"%g",&vd[2]);
      }
    else if (!strcmp(token,"voxelSize"))
      {
      token = strtok(NULL,";");
      int vs;
      sscanf(token,"%d",&vs);

      if (vs != sizeof(T))
        {
        throw "Incorrect voxel size, can not load image!";
        }
      }
    else if (!strcmp(token,"cubePadding"))
      {
      int pc;
      token = strtok(NULL,";");
      sscanf(token,"%d",&pc);
      paddedCubes = (pc==1);
      }
    token = strtok(NULL,"=");
    }

  delete header;

  // Initialize the dimensions array
  for (int i=0;i<3;i++)
    dim[i] = cd[i] * res;

  // Resize the cube
  cube.resize(cd[0],cd[1],cd[2]);

  // Read the cubes
  for (int j=0;j<cube.sizeOfCube();j++)
    {
    int cw = (paddedCubes) ? res+1 : res;
    cube.root()[j].resize(cw,cw,cw);
    fread(cube.root()[j].root(),sizeof(T),cube.root()[j].sizeOfCube(),f);       
    }

  // Close the file
  fclose(f);

  // If the cubes that we read were not padded, we need to unpack and repack
  if (!paddedCubes)
    {
    T *tmp = new T[dim[0]*dim[1]*dim[2]];
    unpackCubes(tmp);
    populateCubes(tmp);
    delete tmp;
    }

  // Generate a transform, maybe later load it from the image
  makeDefaultTransform(vd);

  // Create a slice from the image
  allocateSlice();
}

void AbstractImage3D::makeDefaultTransform(SMLVec3f vd) {
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

inline float AbstractImage3D::doInterpolation(const F32vec4 &r0,const F32vec4 &in0,const F32vec4 &in1) {
  // Registers we are going to use
  F32vec4 r1,r2,r3,r4,r5,r6,r8,i0,i1;                         //      3       2       1       0

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
  return i0[0];
}

/**
 * This method retrieves eight voxels around the position xyz in object space
 * r0:  xyz inside the voxel (between 0,1)
 * r1:  first 4 voxels
 * r2:  last 4 voxels
 */
template <class T>
bool ImageCube<T>::getEightVoxels(float x,float y,float z,F32vec4 &r0,F32vec4 &r1,F32vec4 &r2) {
  F32vec4 r3;                                         //      3       2       1       0
  Is16vec4 i0,i1,i2;
  Is32vec2 j0,j1;

  // Store the voxel  
  r2 = _mm_load_ps(mmData);
  r3 = _mm_load_ps(mmData+4); 
  r1 = _mm_load_ps(mmData+8);

  r0 = _mm_loadu_ps(&x);                              //      ?       z       y       x
  r0 *= r2;                                           //      ?       z       y       x   
  r0 += r3;                                           //      ?       z       y       x

  // Check whether the voxel falls inside the image cube  
  r1 -= r0;

  // These two integers reflect the sign.  They should both be 0 or we are outside of the image
  int maskMin = move_mask(r0);
  int maskMax = move_mask(r1);

  // If the voxel is outside, use special method
  if ((maskMin | maskMax) & 0x07)
    {
    return false;
    }

  // Round down the values and cast back to floating point
  r1 = _mm_shuffle_ps(r0,r0,0x0E);                    //      x       x       ?       z   
  j0 = F32vec4ToIs32vec2(r0);                         //      tr(y)   tr(x)
  j1 = F32vec4ToIs32vec2(r1);                         //      ?       tr(z)
  r1 = Is32vec2ToF32vec4(r1,j1);                      //      ?       ?       ?       tr(z)
  r1 = _mm_shuffle_ps(r1,r1,0x00);                    //      tr(z)   tr(z)   tr(z)   tr(z)   
  r1 = Is32vec2ToF32vec4(r1,j0);                      //      tr(z)   tr(z)   tr(y)   tr(x)   

  // Take integers to shorts
  i0 = (Is16vec4)pack_sat(j0,j1);

  // We still need the voxel index and cube index 
  // Shift right to find the cube index of the starting voxel
  i1 = i0 >> ImageCube::LDCUBE;

  // Prefetch the cube (is this really worth while?)
  int xc = _MM_4W(0,i1);int yc = _MM_4W(1,i1);int zc = _MM_4W(2,i1);
  DataCube<T> &C = cube(xc,yc,zc);
  _mm_prefetch((const char *)C.root(),_MM_HINT_T0);

  // Shift back and subtract to get the in-cube index of the voxel
  i2 = i1 << ImageCube::LDCUBE;
  i0 -= i2;

  // Increment i0 by 1
  i2 = (Is16vec4)_mm_set1_pi16(1);
  i2 += i0;

  // At this point, voxel cube i1 with offset i0 points to V000 and cube i2 with offset i4 points to V111
  // This is the slow point...    
  int x0 = _MM_4W(0,i0);int y0 = _MM_4W(1,i0);int z0 = _MM_4W(2,i0);
  int x1 = _MM_4W(0,i2);int y1 = _MM_4W(1,i2);int z1 = _MM_4W(2,i2);

  // Clear the MMX registers
  empty();

  // This places the fractional part of the voxel computation into r0
  r0 -= r1;

  // Clear the MMX registers
  _mm_empty();

  // Get the data from the cube
  float c000 = (float) C(x0,y0,z0);
  float c001 = (float) C(x0,y0,z1);
  float c010 = (float) C(x0,y1,z0);
  float c011 = (float) C(x0,y1,z1);
  float c100 = (float) C(x1,y0,z0);
  float c101 = (float) C(x1,y0,z1);
  float c110 = (float) C(x1,y1,z0);
  float c111 = (float) C(x1,y1,z1);

  // Create the return values
  r1 = F32vec4(c001,c010,c100,c000);
  r2 = F32vec4(c111,c011,c101,c110);

  // Clear the MMX registers
  _mm_empty();

  return true;
}

template <class T>
float ImageCube<T>::interpolateVoxel(float x,float y,float z) {
  F32vec4 r0,r1,r2;

  bool in = getEightVoxels(x,y,z,r0,r1,r2);
  return(in) ? doInterpolation(r0,r1,r2) : handleOutsideVoxel(x,y,z);
}

template <class T>
void ImageCube<T>::interpolateVoxelGradient(float xs,float ys,float zs, float *G) {
  F32vec4 r0,r1,r2,r3;                                    //      3       2       1       0
  F32vec4 r4,r5,r6,r7,r8;                                 //      3       2       1       0
  Is16vec4 i0,i1,i2,i3,i4,i5,i6;
  Is32vec2 j0,j1;

  if (!getEightVoxels(xs,ys,zs,r0,r4,r5))
    {
    G[0] = G[1] = G[2] = 0.0f;
    return;
    }

  // Now compute the interpolation
  r1 = _mm_set1_ps(1.0f);                             //  1           1           1           1
  r1 = _mm_sub_ps(r1,r0);                             //  ?           1-z         1-y         1-x
  r2 = _mm_shuffle_ps(r0,r0,0xD2);                    //  ?           y           x           z
  r0 = _mm_shuffle_ps(r0,r0,0xC9);                    //  ?           x           z           y
  r3 = _mm_shuffle_ps(r1,r1,0xD2);                    //  ?           Y           X           Z
  r1 = _mm_shuffle_ps(r1,r1,0xC9);                    //  ?           X           Z           Y

  r6 = _mm_mul_ps(r0,r2);                             //  ?           xy          zx          yz
  r7 = _mm_mul_ps(r1,r3);                             //  ?           XY          ZX          YZ
  r1 = _mm_mul_ps(r1,r2);                             //  ?           Xy          Zx          Yz
  r0 = _mm_mul_ps(r0,r3);                             //  ?           xY          zX          yZ

  r2 = _mm_shuffle_ps(r4,r4,0x00);                    //  ?           V000        V000        V000
  r3 = _mm_shuffle_ps(r4,r4,0x1B);                    //  ?           V001        V010        V100
  r3 = _mm_sub_ps(r3,r2);                             //  ?           V001-V000   V010-V000   V100-V000
  r7 = _mm_mul_ps(r7,r3);                             //  ?           *           *           *

  r2 = _mm_shuffle_ps(r4,r4,0x2D);                    //  ?           V010        V100        V001
  r3 = _mm_shuffle_ps(r5,r5,0x09);                    //  ?           V011        V110        V101
  r3 = _mm_sub_ps(r3,r2);                             //  ?           V011-V010   V110-V100   V101-V001
  r1 = _mm_mul_ps(r1,r3);                             //  ?           *           *           *

  r2 = _mm_shuffle_ps(r4,r4,0x36);                    //  ?           V100        V001        V010
  r3 = _mm_shuffle_ps(r5,r5,0x12);                    //  ?           V101        V011        V110
  r3 = _mm_sub_ps(r3,r2);                             //  ?           V101-V100   V011-V001   V110-V010
  r0 = _mm_mul_ps(r0,r3);                             //  ?           *           *           *

  // r2 = _mm_shuffle_ps(r5,r5,0x36);                 //  ?           V110        V101        V011
  r3 = _mm_shuffle_ps(r5,r5,0xFF);                    //  ?           V111        V111        V111
  r5 = _mm_sub_ps(r3,r5);                             //  ?           V111-V110   V111-V101   V111-V011
  r6 = _mm_mul_ps(r6,r5);                             //  ?           *           *           *

  // Add the results
  r1 = _mm_add_ps(r1,r7);
  r0 = _mm_add_ps(r0,r6);
  r0 = _mm_add_ps(r0,r1);

  // Scale by the max value
  // r0 = _mm_mul_ps(r0,r8);

  // SMLVec3f T;
  // float G0 = ((V100-V000)*Y*Z + (V101-V001)*Y*z + (V110-V010)*y*Z + (V111-V011)*y*z) * mmData[0];
  // float G1 = ((V010-V000)*Z*X + (V110-V100)*Z*x + (V011-V001)*z*X + (V111-V101)*z*x) * mmData[0];
  // float G2 = ((V001-V000)*X*Y + (V011-V010)*X*y + (V101-V100)*x*Y + (V111-V110)*x*y) * mmData[0];

  // G[0] = ((V100-V000)*Y*Z + (V101-V001)*Y*z + (V110-V010)*y*Z + (V111-V011)*y*z) * mmData[0];
  // G[1] = ((V010-V000)*X*Z + (V011-V001)*X*z + (V110-V100)*x*Z + (V111-V101)*x*z) * mmData[0];
  // G[2] = ((V001-V000)*X*Y + (V011-V010)*X*y + (V101-V100)*x*Y + (V111-V110)*x*y) * mmData[0];
  G[0] = r0[0];
  G[1] = r0[1];
  G[2] = r0[2];
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

SMLVec3f convertDistanceToRGB(float pixel) {
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

// This method loads the image from a zyx ordered array of shorts
template<class T, class S> 
void makePaddedDataCube(int dx,int dy,int dz,int res,S *array,DataCube<T> &DC) {

  // Image must be padded to 16 voxels in all dimensions
  int d0 = roundup(dx,res);
  int d1 = roundup(dy,res);
  int d2 = roundup(dz,res);

  // Rescale the image
  DC.resize(d0,d1,d2,0);

  // Check if there was any roundup performed and if so, create a larger image and copy to it
  for (int z=0;z<dz;z++)
    {
    for (int y=0;y<dy;y++)
      {
      for (int x=0;x<dx;x++)
        {
        DC(x,y,z) = *array;
        array++;
        }
      }
    }
}

// This method loads the image from a zyx ordered array of shorts
template<class T> 
void ImageCube<T>::loadFromPaddedCube(const T *pixels, const int *size) 
{
// Image must be padded to 16 voxels in all dimensions
  dim[0] = size[0];
  dim[1] = size[1];
  dim[2] = size[2];

  // Compute the CD values as well
  cd[0] = dim[0] >> LDCUBE;
  cd[1] = dim[1] >> LDCUBE;
  cd[2] = dim[2] >> LDCUBE;

  // Resize the array
  cube.resize(cd[0],cd[1],cd[2]);
  populateCubes(pixels);

  // Create a slice from the image
  allocateSlice();
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


template<class T>
void ImageCube<T>::loadFromITKReadableFile(const char *giplFile, std::ostream &out) 
{
  typedef itk::Image<T,3> ImageType;
  typedef ImageType::RegionType RegionType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  // Load an image
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(giplFile);
  reader->Update();

  // Get the dimensions
  ImageType::SizeType size = reader->GetOutput()->GetBufferedRegion().GetSize();
  
  // Create a padded cube
  DataCube<T> imageCube;
  makePaddedDataCube(size[0],size[1],size[2],res,
                     reader->GetOutput()->GetBufferPointer(),imageCube);
  loadFromPaddedCube(imageCube.root(),imageCube.size());
    
  // Compute the voxel transform
  makeDefaultTransform(SMLVec3f(1.0f,1.0f,1.0f));
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


DataCube<unsigned short> test4(16,16,16);
DataCube<float> test5(16,16,16);
DataCube<unsigned char> test6(16,16,16);



#pragma auto_inline(on)