#include "mathlink.h"
#include <cluto.h>

void CLUTOClusterDirect(
  int *ai, long n_ai,
  int *a, long n_a,
  double *w, long n_w,
  const char *crit, int ntials, int niter, int seed, int nclust)  
{
  int nVertices = (int) n_ai - 1;
  int i;

  float *fw = new float[n_a];
  for(i = 0; i < n_a; i++) fw[i] = (float) w[i];

  // Decrement all the indices
  for(i = 0; i < n_a; i++) a[i]--;
  for(i = 0; i < n_ai; i++) ai[i]--;

  int *outPart = new int[nVertices];
  
  int icrit = CLUTO_CLFUN_E1;
  if(0==strcmp(crit,"g1"))
    icrit = CLUTO_CLFUN_G1;
  else if(0==strcmp(crit,"g1p"))
    icrit = CLUTO_CLFUN_G1P;
  else if(0==strcmp(crit,"h1"))
    icrit = CLUTO_CLFUN_H1;
  else if(0==strcmp(crit,"h2"))
    icrit = CLUTO_CLFUN_H2;
  else if(0==strcmp(crit,"i1"))
    icrit = CLUTO_CLFUN_I1;
  else if(0==strcmp(crit,"i2"))
    icrit = CLUTO_CLFUN_I2;

  CLUTO_SP_ClusterDirect(nVertices,ai,a,fw,icrit,ntials,niter,seed,0,nclust,outPart);

  // Return to mathematica
	MLPutFunction(stdlink,"List",1);
	MLPutIntegerList(stdlink,outPart,nVertices);

  // Delete stuff
  delete outPart;
  delete fw;
}

void CLUTOGraphClusterRB(
  int *ai, long n_ai,
  int *a, long n_a,
  double *w, long n_w,
  int nnbr, double eprun, double vprun, 
  int mincmp, int ntials, int seed, 
  const char *strat, int nclust)
{
  int nVertices = (int) n_ai - 1;
  int i;

  float *fw = new float[n_a];
  for(i = 0; i < n_a; i++) fw[i] = (float) w[i];

  // Decrement all the indices
  for(i = 0; i < n_a; i++) a[i]--;
  for(i = 0; i < n_ai; i++) ai[i]--;

  int *outPart = new int[nVertices];
  float outEdgecut;
  
  int istrat = CLUTO_CSTYPE_BESTFIRST;
  if(0==strcmp(strat,"SelectBest"))
    istrat = CLUTO_CSTYPE_BESTFIRST;
  else if(0==strcmp(strat,"SelectLargest"))
    istrat = CLUTO_CSTYPE_LARGEFIRST;
  else if(0==strcmp(strat,"SelectLargestSubspace"))
    istrat = CLUTO_CSTYPE_LARGESUBSPACEFIRST;

  CLUTO_SP_GraphClusterRB(nVertices, ai, a, fw,
              nnbr, (float) eprun, (float) vprun, mincmp, ntials, seed, istrat, 0, 
              nclust, outPart, &outEdgecut);

  // Return to mathematica
	MLPutFunction(stdlink,"List",2);
	MLPutReal(stdlink,outEdgecut);
	MLPutIntegerList(stdlink,outPart,nVertices);

  // Delete stuff
  delete outPart;
  delete fw;
}



int main(int argc, char *argv[])
{
	return MLMain(argc, argv);
}