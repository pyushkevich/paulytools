/**
 * ImageToMETIS : Takes a binary image and constructs a METIS graph
 * that can be partitioned according to some rules
 */
#include "ReadWriteImage.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include <map>
#include <list>
#include <cstdio>

using namespace std;

unsigned long GetIndexOffsetInImage(
  const itk::ImageRegion<3> &region, 
  const itk::Index<3> &index)
{
  unsigned long o = 0;
  o += (index[2] - region.GetIndex(2)) * region.GetSize(2) * region.GetSize(1);
  o += (index[1] - region.GetIndex(1)) * region.GetSize(1);
  o += (index[0] - region.GetIndex(0));
  return o;
}

int main(int argc, char *argv[])
{
  if(argc != 4)
    {
    cerr << "usage binary2metis image.hdr outmesh.metis tempfile.txt" << endl;
    return -1;
    }
  
  // Read the image
  typedef itk::Image<short,3> ImageType;
  ImageType::Pointer img = ImageType::New();
  ReadImage(img,argv[1]);
  
  typedef list<unsigned int> VertexInfo;
  typedef list<VertexInfo> VertexList;

  VertexList lVertex;

  // Keep track of numbers of vertices and edges
  unsigned int nVertices = 0, nEdges = 0;

  // Create a map from absolute index to relative index
  map<unsigned int,unsigned int> xIndexMap;
  map<unsigned int,itk::Index<3> > xInverseMap;
  
  // Traverse all pixels in the image
  itk::ImageRegionConstIteratorWithIndex<ImageType> it(img,img->GetBufferedRegion());
  while(!it.IsAtEnd())
    {
    ImageType::IndexType idx = it.GetIndex();
    if(img->GetPixel(idx) != 0)
      {
      // Get the absolute offset
      unsigned int offset = 
        (unsigned int)GetIndexOffsetInImage(img->GetBufferedRegion(),idx);

      // Start a new vertex info
      VertexInfo vi;
      
      // Look at the neighbor pixels
      for(unsigned int d=0;d<3;d++)
        {
        for(int k = -1; k <= 1; k+=2)
          {
          ImageType::IndexType nidx = idx;
          nidx[d] += k;
          if(img->GetBufferedRegion().IsInside(nidx))
            {
            if(img->GetPixel(nidx) != 0)
              {
              // Pixel is a legitimate neighbor
              unsigned int iNbr = 
                (unsigned int) GetIndexOffsetInImage(img->GetBufferedRegion(),nidx);
              vi.push_back(iNbr);
              nEdges++;
              }
            }
          }
        }

      // Only consider vertex if it has neighbors (must be connected)
      if(vi.size())
        {
        // Add index to index map
        xIndexMap[offset] = nVertices;
        xInverseMap[nVertices] = idx;
        lVertex.push_back(vi);
        nVertices++;
        }
      }
    ++it;
    }

  // Now we have a list of lists representing vertices. Write it to a metis file
  FILE *fout = fopen(argv[2],"wt");
  FILE *fmap = fopen(argv[3],"wt");

  // The number of 'phantom' voxels  
  unsigned int nPhantom = 0;

  fprintf(fout,"%d %d 11\n",nVertices + nPhantom,nEdges/2 + nPhantom);
  for(VertexList::iterator it = lVertex.begin();it!=lVertex.end();++it)
    {
    // Write the next line out
    fprintf(fout,"1 ");
    for(VertexInfo::iterator iv = it->begin();iv!=it->end();++iv)
      {
      fprintf(fout,"%d 10 ",xIndexMap[*iv]+1);
      }
    fprintf(fout,"\n");
    }

  // Create 'phantom vertices'
  for(unsigned int i=0;i<nPhantom;i++)
    {
    fprintf(fout,"1 %d 1 %d 1\n",((i-1) % nPhantom) + nVertices, ((i+1) % nPhantom) + nVertices);
    }

  for(unsigned int j=0;j<nVertices;j++)
    {
    itk::Index<3> idx = xInverseMap[j];
    fprintf(fmap,"%d %d %d\n",(unsigned int)idx[0], (unsigned int)idx[1], (unsigned int)idx[2] );
    }
  
  fclose(fout);
  fclose(fmap);

  return 0;
}
