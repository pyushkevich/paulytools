#include "mathlink.h"
#include "BasicImaging/ShortestPath.h"
#include <map>

// A data chuck stored for each SP problem
typedef DijkstraShortestPath<double> SPath;
struct SPProblem {
  SPath *sp;
  unsigned int *xAdjIndex, *xAdj;
  double *xWgt;
};

// A container for shortest path algorithm instances
typedef std::map<int,SPProblem> SPathMap;

SPathMap gSPathMap;

// Shortest path functions
int BuildDijkstra(
  int *xAdjIndex, long xAdjIndexLen, 
  int *xAdj, long xAdjLen, 
  double *wAdj, long wAdjLen)
  {
  // Create an SP problem
  SPProblem P;

  // Copy the data
  P.xAdjIndex = new unsigned int[xAdjIndexLen];
  P.xAdj = new unsigned int[xAdjLen];
  P.xWgt = new double[xAdjLen];

  for(long i=0;i<xAdjIndexLen;i++)
    P.xAdjIndex[i] = xAdjIndex[i];

  for(long j=0;j<xAdjLen;j++)
    {
    P.xAdj[j] = xAdj[j];
    P.xWgt[j] = wAdj[j];
    }

  // Create a new shortest path algo
  P.sp = new SPath(xAdjIndexLen-1,P.xAdjIndex,P.xAdj,P.xWgt);

  // Create a handle
  int handle = rand();
  while(gSPathMap.find(handle) != gSPathMap.end())
    handle = rand();

  // Add the pointer to the map
  gSPathMap[handle] = P;
  
  // Return the pointer as a long
  return handle;
  }

int ReleaseDijkstra(int handle)
  {
  if(gSPathMap.find(handle) != gSPathMap.end())
    {
    delete gSPathMap[handle].sp;
    delete gSPathMap[handle].xAdj;
    delete gSPathMap[handle].xAdjIndex;
    delete gSPathMap[handle].xWgt;
    gSPathMap.erase(handle);
    }
  return 0;
  }

double DistanceDijkstra(int handle, int xStart, int xEnd)
  {
  if(gSPathMap.find(handle) == gSPathMap.end())
    return -1.0;

  SPath *sp = gSPathMap[handle].sp;
  sp->ComputePathsFromSource(xStart,xEnd);
  return sp->GetDistanceArray()[xEnd];
  }

// Functions defined by the METIS library
typedef int idxtype;
extern "C" {
  void METIS_WPartGraphKway(
  	int *, idxtype *, idxtype *, idxtype *, idxtype *, 
  	int *, int *, int *, float *, int *, int *, idxtype *);

  void METIS_WPartGraphRecursive(
  	int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, 
    int *, int *, float *, int *, int *, idxtype *);
}

void METISPartition( 
  int *xAdjIndex, long xAdjIndexLen, 
  int *xAdj, long xAdjLen, 
  int *wVertex, long wVertexLen, 
  int *wAdj, long wAdjLen, 
  double *wCluster, long wClusterLen, 
  int option, int method)
{
  // Set up the parameters for the call. We have to convert double arrays to float.
  int nVertices = (int) xAdjIndexLen - 1;
  int nClusters = (int) wClusterLen;
  float *fwCluster = new float[wClusterLen];
  int outEdgecut, *outPart = new int[nVertices];
  int wFlag = 3, fmtFlag = 1;

  // Copy the arrays
  for(long i = 0; i < wClusterLen; i++) fwCluster[i] = (float) wCluster[i];

  // Run the METIS method
  METIS_WPartGraphRecursive(
  	&nVertices, xAdjIndex, xAdj, wVertex, wAdj,
    &wFlag, &fmtFlag, &nClusters, fwCluster,
    &option, &outEdgecut, outPart);

  // Return to mathematica
	MLPutFunction(stdlink,"List",2);
	MLPutInteger(stdlink,outEdgecut);
	MLPutIntegerList(stdlink,outPart,nVertices);

  // Delete stuff
  delete outPart;
  delete fwCluster;
}

int main(int argc, char *argv[])
{
	return MLMain(argc, argv);
}