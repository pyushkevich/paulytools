#ifndef __ImageCubeITK_txx_
#define __ImageCubeITK_txx_

#include "itkImageRegionConstIteratorWithIndex.h"

template<typename TPixel>
void
ImageCubeITK<TPixel>
::SetImage(ImageType *imgSource, const TPixel &xBackground)
{
  // Get the dimensions of the image
  typename ImageType::IndexType idx = imgSource->GetBufferedRegion().GetIndex();
  typename ImageType::SizeType size = imgSource->GetBufferedRegion().GetSize();

  // Allocate the appropriate number of cubes
  for(unsigned int d = 0; d < 3; d++)
    {
    unsigned int sz = size[d];
    if(sz % (1 << LDCUBE) > 0)
      cd[d] = 1 + (sz >> LDCUBE);
    else
      cd[d] = sz >> LDCUBE;
    dim[d] = cd[d] << LDCUBE;
    }

  // Resize the array and each of the internal cubes
  cube.resize(cd[0], cd[1], cd[2]);
  for(unsigned int iCube = 0; iCube < cube.sizeOfCube(); iCube++)
    { cube(iCube).resize(res+1, res+1, res+1, xBackground); }

  // Iterate through the image volume, filling everything with image data
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> Iterator;
  Iterator it(imgSource, imgSource->GetBufferedRegion());
  while(!it.IsAtEnd())
    {
    unsigned int ix = (unsigned int) (it.GetIndex()[0] - idx[0]);
    unsigned int iy = (unsigned int) (it.GetIndex()[1] - idx[1]);
    unsigned int iz = (unsigned int) (it.GetIndex()[2] - idx[2]);

    cube(ix >> LDCUBE, iy >> LDCUBE, iz >> LDCUBE)
      (ix % LDCUBE, iy % LDCUBE, iz % LDCUBE) = it.Value();
    
    ++it;
    }

  // Fill the border cells
  for(unsigned int cx = 0; cx < cd[0] - 1; cx++)
    for(unsigned int cy = 0; cy < cd[1] - 1; cy++)
      for(unsigned int cz = 0; cz < cd[2] - 1; cz++)
        {
        DataCube<TPixel> &dc = cube(cx, cy, cz);
        for(unsigned int u = 0; u < res; u++)
          {
          for(unsigned int v = 0; v < res; v++)
            {
            dc(u, v, res) = cube(cx, cy, cz+1)(u, v, 0);
            dc(u, res, v) = cube(cx, cy+1, cz)(u, 0, v);
            dc(res, u, v) = cube(cx+1, cy, cz)(0, u, v);            
            }
          dc(u, res, res) = cube(cx, cy+1, cz+1)(u, 0, 0);
          dc(res, u, res) = cube(cx+1, cy, cz+1)(0, u, 0);
          dc(res, res, u) = cube(cx+1, cy+1, cz)(0, 0, u);
          }
        dc(res, res, res) = cube(cx+1, cy+1, cz+1)(0, 0, 0);
        }

  // Set the voxel size
  SMLVec3f xVoxel(imgSource->GetSpacing()[0],
    imgSource->GetSpacing()[1],
    imgSource->GetSpacing()[2]);

  // Compute the default transform
  makeDefaultTransform(xVoxel);
}

#endif // __ImageCubeITK_txx_
