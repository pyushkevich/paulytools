#include "mathlink.h"

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