#ifndef __imaging_txx_
#define __imaging_txx_

#include <itkImage.h>
#include <itkCastImageFilter.h>
#include <itkImageFileReader.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <iostream>
#include <sstream>

using namespace std;

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
  typedef typename ImageType::RegionType RegionType;
  typedef itk::DiscreteGaussianImageFilter<FloatImageType,FloatImageType> FilterType;
  typedef itk::CastImageFilter<ImageType,FloatImageType> FloatCasterType;
  typedef itk::CastImageFilter<FloatImageType,ImageType> ReverseCasterType;

  // Initialize regions
  typename ImageType::SizeType size = {{dim[0],dim[1],dim[2]}};
  typename ImageType::RegionType region(size);
  
  // Create a flat image
  typename ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->Allocate();

  // Dump the image cube into a flat image
  unpackCubes(image->GetBufferPointer());

  // Create a variance vector
  SMLVec3f variance(
    sigma[0] * sigma[0], sigma[1] * sigma[1], sigma[2] * sigma[2]);

  // Construct and run a pipeline
  typename FloatCasterType::Pointer floatCaster = FloatCasterType::New();
  floatCaster->SetInput(image);

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput(floatCaster->GetOutput());
  filter->SetVariance(variance.data_block());
  
  typename ReverseCasterType::Pointer revCaster = ReverseCasterType::New();
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

template<class T>
void ImageCube<T>::loadFromITKReadableFile(const char *giplFile, std::ostream &out) 
{
  typedef itk::Image<T,3> ImageType;
  typedef typename ImageType::RegionType RegionType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  // Load an image
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(giplFile);
  reader->Update();

  // Get the dimensions
  typename ImageType::SizeType size = reader->GetOutput()->GetBufferedRegion().GetSize();
  
  // Create a padded cube
  DataCube<T> imageCube;
  makePaddedDataCube(size[0],size[1],size[2],res,
                     reader->GetOutput()->GetBufferPointer(),imageCube);
  loadFromPaddedCube(imageCube.root(),imageCube.size());

  // Set the voxel size
  setVoxelSize(
    reader->GetOutput()->GetSpacing()[0],
    reader->GetOutput()->GetSpacing()[1],
    reader->GetOutput()->GetSpacing()[2]);
  out << "Setting voxel size to " << getVoxelSize() << endl;
    
  // Compute the voxel transform
  makeDefaultTransform(SMLVec3f(1.0f,1.0f,1.0f));
} 

/**
 * This method retrieves eight voxels around the position xyz in object space
 * r0:  xyz inside the voxel (between 0,1)
 * r1:  first 4 voxels
 * r2:  last 4 voxels
 */
template <class T>
inline bool 
ImageCube<T>
::getEightVoxels(float x,float y,float z,__m128 &out0,__m128 &out1,__m128 &out2) 
{
  register __m128 r0,r1,r2,r3;                        //      3       2       1       0
  __m64 i0,i1,i2;
  __m64 j0,j1;

  // Store the voxel  
  r2 = _mm_load_ps(mmData);
  r3 = _mm_load_ps(mmData+4); 
  r1 = _mm_load_ps(mmData+8);

  r0 = _mm_set_ps(0, z, y, x);                        //      ?       z       y       x
  r0 = _mm_mul_ps(r0,r2);                             //      ?       z       y       x   
  r0 = _mm_add_ps(r0,r3);                             //      ?       z       y       x

  // Check whether the voxel falls inside the image cube  
  r1 = _mm_sub_ps(r1,r0);

  // These two integers reflect the sign.  They should both be 0 or we are outside of the image
  int maskMin = _mm_movemask_ps(r0);
  int maskMax = _mm_movemask_ps(r1);

  // If the voxel is outside, use special method
  if ((maskMin | maskMax) & 0x07) return false;

  // Set the rounding state to truncation
  int flgRound = _mm_getcsr();
  _mm_setcsr(flgRound | 0x00006000);
  
  // Round down the values and cast back to floating point
  r1 = _mm_shuffle_ps(r0,r0,0x0E);                    //      x       x       ?       z   
  j0 = _mm_cvtps_pi32(r0);                            //      tr(y)   tr(x)
  j1 = _mm_cvtps_pi32(r1);                            //      ?       tr(z)
  r1 = _mm_cvtpi32_ps(r1,j1);                         //      ?       ?       ?       tr(z)
  r1 = _mm_shuffle_ps(r1,r1,0x00);                    //      tr(z)   tr(z)   tr(z)   tr(z)
  r1 = _mm_cvtpi32_ps(r1,j0);                         //      tr(z)   tr(z)   tr(y)   tr(x)

  // Take integers to shorts
  i0 = _mm_packs_pi32(j0, j1);

  // We still need the voxel index and cube index 
  // Shift right to find the cube index of the starting voxel
  i1 = _mm_srli_pi16(i0, ImageCube::LDCUBE);

  // This places the fractional part of the voxel computation into r0
  out0 = _mm_sub_ps(r0,r1);                           //      fr(z)   fr(z)   fr(y)   fr(x)

  // Clear the MMX registers
  _mm_empty();

  // Prefetch the cube (is this really worth while?)
  DataCube<T> *C = cube.voxel(*((short*)&i1 + 0),*((short*)&i1 + 1),*((short*)&i1 + 2));
  _mm_prefetch((const char *)C->root(),_MM_HINT_T0);

  // Shift back and subtract to get the in-cube index of the voxel
  i2 = _mm_slli_pi16(i1,ImageCube::LDCUBE);
  i0 = _mm_sub_pi16(i0,i2);

  // Restore rounding flags
  _mm_setcsr(flgRound | 0x00006000);

  // Clear the MMX registers
  _mm_empty();

  // Get the data from the cube
  ALIGN_PRE float voxels[8] ALIGN_POST;
  C->getEightVoxelCube(*((short*)&i0 + 0),*((short*)&i0 + 1),*((short*)&i0 + 2),voxels);
  
  // Store the data 
  out1 = _mm_load_ps(voxels);
  out2 = _mm_load_ps(voxels + 4);

  return true;
}

template <class T>
float ImageCube<T>::interpolateVoxel(float x,float y,float z) {
  __m128 z0,z1,z2;
  float val;

  bool in = getEightVoxels(x,y,z,z0,z1,z2);
  if(in)
    {
    // doInterpolation(r0,r1,r2,&val);

    // Registers we are going to use
    __m128 r1,r2,r3,r4,r5,r6,r8,i0,i1;                      //  3       2       1       0

    // Perform the interpolation
    r1 = _mm_shuffle_ps(z0,z0,0x19);                        //  x       y       z       y
    r2 = _mm_shuffle_ps(z0,z0,0x62);                        //  y       z       x       z
    r3 = _mm_mul_ps(r1,r2);                                 //  xy      yz      zx      yz
    r4 = _mm_shuffle_ps(z0,z0,0x84);                        //  z       x       y       x           
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
    i0 = _mm_mul_ps(z1,r8);                                 
    i1 = _mm_mul_ps(z2,r3);                                 

    // Now add up all elements
    i0 = _mm_add_ps(i0,i1);
    i1 = _mm_shuffle_ps(i0,i0,0xb1);
    i0 = _mm_add_ps(i0,i1);
    i1 = _mm_shuffle_ps(i0,i0,0x4e);
    i0 = _mm_add_ps(i0,i1);

    // The result is in r1
    _mm_store_ss(&val, i0);

    return val;
    }
  else
    { return handleOutsideVoxel(x,y,z); }
}

  template <class T>
void ImageCube<T>::interpolateVoxelGradient(float xs,float ys,float zs, float *G) 
{
  __m128 r0,r1,r2,r3;                                    //      3       2       1       0
  __m128 r4,r5,r6,r7;                                 //      3       2       1       0

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

  // Store the gradient
  ALIGN_PRE float gtemp[4] ALIGN_POST;
  _mm_store_ps(gtemp,r0);

  G[0] = gtemp[0];
  G[1] = gtemp[1];
  G[2] = gtemp[2];
}

#endif // __imaging_txx_
