#include "mspline.h"
#include <iostream>
#include <ctime>
#include <matrix.h>

using namespace std;

inline SMLVec3f& cast3f(const SMLVec4f &in) {
  return *((SMLVec3f*)&in);
}

// Perform simultaneous dot product of four vectors with another vector
inline __m128 mul_4x4_4x1(const __m128 &R0,const __m128 &R1,const __m128 &R2,const __m128 &R3,const __m128 &RB) {
  __m128 TM,A0,A1,A2,A3,B0,B1,B2,B3;

  // Perform the vector multiplication
  A0 = _mm_mul_ps(R0,RB);
  A1 = _mm_mul_ps(R1,RB);
  A2 = _mm_mul_ps(R2,RB);
  A3 = _mm_mul_ps(R3,RB);

  B1 = _mm_shuffle_ps(A0,A1,0x44);
  TM = _mm_shuffle_ps(A2,A3,0x44);
  B0 = _mm_shuffle_ps(B1,TM,0x88);
  B1 = _mm_shuffle_ps(B1,TM,0xDD);

  B3 = _mm_shuffle_ps(A0,A1,0xee);
  TM = _mm_shuffle_ps(A2,A3,0xee);
  B2 = _mm_shuffle_ps(B3,TM,0x88);
  B3 = _mm_shuffle_ps(B3,TM,0xDD);

  // Perform the additions
  B0 = _mm_add_ps(B0,B1);
  B2 = _mm_add_ps(B2,B3);

  // Perform the final addition
  B0 = _mm_add_ps(B0,B2);

  // Return the result in B0
  return B0;
}

/**
 * Medial point accesors
 */
SMLVec3f& MedialPoint::X() {
  return *((SMLVec3f*)(&F));
}

SMLVec3f& MedialPoint::Xu() {
  return *((SMLVec3f*)(&Fu));
}

SMLVec3f& MedialPoint::Xv() {
  return *((SMLVec3f*)(&Fv));
}

SMLVec3f& MedialPoint::Xuu() {
  return *((SMLVec3f*)(&Fuu));
}

SMLVec3f& MedialPoint::Xuv() {
  return *((SMLVec3f*)(&Fuv));
}

SMLVec3f& MedialPoint::Xvv() {
  return *((SMLVec3f*)(&Fvv));
}

SMLVec3f& MedialPoint::N() {
  return *((SMLVec3f*)(&F3));
}

SMLVec3f& MedialPoint::NRaw() {
  return *((SMLVec3f*)(&F3raw));
}

float& MedialPoint::R() {
  return F[3];
}

float& MedialPoint::Ru() {
  return Fu[3];
}

float& MedialPoint::Rv() {
  return Fv[3];
}

float& MedialPoint::Ruu() {
  return Fuu[3];
}

float& MedialPoint::Ruv() {
  return Fuv[3];
}

float& MedialPoint::Rvv() {
  return Fvv[3];
}

const SMLVec3f& MedialPoint::X() const {
  return *((SMLVec3f*)(&F));
}

const SMLVec3f& MedialPoint::Xu() const {
  return *((SMLVec3f*)(&Fu));
}

const SMLVec3f& MedialPoint::Xv() const {
  return *((SMLVec3f*)(&Fv));
}

const SMLVec3f& MedialPoint::Xuu() const {
  return *((SMLVec3f*)(&Fuu));
}

const SMLVec3f& MedialPoint::Xuv() const {
  return *((SMLVec3f*)(&Fuv));
}

const SMLVec3f& MedialPoint::Xvv() const {
  return *((SMLVec3f*)(&Fvv));
}

const SMLVec3f& MedialPoint::N() const {
  return *((SMLVec3f*)(&F3));
}

const SMLVec3f& MedialPoint::NRaw() const {
  return *((SMLVec3f*)(&F3raw));
}

float MedialPoint::R() const {
  return F[3];
}

float MedialPoint::Ru() const {
  return Fu[3];
}

float MedialPoint::Rv() const {
  return Fv[3];
}

float MedialPoint::Ruu() const {
  return Fuu[3];
}

float MedialPoint::Ruv() const {
  return Fuv[3];
}

float MedialPoint::Rvv() const {
  return Fvv[3];
}

// Forward Basis Computation Algorithm, Piegl 1995, p.70
// For now we are only using uniform B-splines.  So this is easier done with a 
// matrix multiplication than anything else
void DynamicBSpline2D::basisJet(int dim,int i,float u,SMLVec4f *N) {
  if (regKnots)
    {
    __m128 r0,r1,r2,r3,r4,r5,r6,r7;
    __m128 M0,M1,M2,M3;
    
    
    SMLVec4f *jm = jetMul[dim];
    float *ub = UB->data_block();

    // Load u                                                   //  3       2       1       0
    r0 = _mm_load_ss(&u);                                       //                          u
    r1 = _mm_load_ss(knots[dim]+i);                             //                          knot
    r0 = _mm_sub_ss(r0,r1);                                     //                          u-knot

    // Prefetch
    _mm_prefetch((const char *)jm,1);
    _mm_prefetch((const char *)ub,1);

    // Compute the vector u^3 u^2 u 1
    r1 = _mm_set_ss(1.0f);
    r1 = _mm_shuffle_ps(r0,r1,0x00);                            //  1       1       u       u
    r2 = _mm_shuffle_ps(r1,r1,0x80);                            //  1       u       u       u
    r3 = _mm_shuffle_ps(r1,r1,0xa8);                            //  1       1       1       u
    r1 = _mm_mul_ps(r1,r3);                                     //  1       1       u       uu
    r2 = _mm_mul_ps(r2,r1);                                     //  1       u       uu      uuu

    // Compute all the products
    r4 = _mm_load_ps(jm[0].data_block());
    r5 = _mm_load_ps(jm[1].data_block());
    r6 = _mm_load_ps(jm[2].data_block());
    r7 = _mm_load_ps(jm[3].data_block());

    // Compute the vector 3*u^2 2u 1 0
    r4 = _mm_mul_ps(r4,r2);
    r5 = _mm_mul_ps(r5,r1);
    r6 = _mm_mul_ps(r6,r3);

    // Load the BSpline matrix
    M0 = _mm_load_ps(ub);
    M1 = _mm_load_ps(ub+4);
    M2 = _mm_load_ps(ub+8);
    M3 = _mm_load_ps(ub+12);

    // That should be it because we will encorporate the BSpline matrix into the control point's matrix
    r0 = mul_4x4_4x1(M0,M1,M2,M3,r4);
    r1 = mul_4x4_4x1(M0,M1,M2,M3,r5);
    r2 = mul_4x4_4x1(M0,M1,M2,M3,r6);

    // Store the results
    _mm_store_ps(N[0].data_block(),r0);
    _mm_store_ps(N[1].data_block(),r1);
    _mm_store_ps(N[2].data_block(),r2);
    }
  else
    {
    // We hardcode p because it's more efficient
    const int p = 3;
    const int n = 2;

    float *U = this->knots[dim];

    // Left and right arrays are not dynallocated for effifiency.  Don't pass p >= 10
    float ndu[p+1][p+1],a[2][p+1],left[p+1],right[p+1];

    float saved,temp;
    int j,k,j1,j2,r;

    ndu[0][0]=1.0f;
    for (j=1; j<=p; j++ )
      {
      left[j]=u-U[i+1-j];
      right[j]=U[i+j]-u;
      saved=0.0f;

      for ( r=0; r<j; r++ )
        {
        ndu[j][r]=right[r+1]+left[j-r];
        temp=ndu[r][j-1]/ndu[j][r];

        ndu[r][j]=saved+right[r+1]*temp;
        saved=left[j-r]*temp;
        }
      ndu[j][j]=saved;
      }

    for ( j=0; j<=p; j++ )
      N[0][j]=ndu[j][p];  // basis function

    // compute derivatives
    for ( r=0; r<=p; r++ )
      {
      int s1=0, s2=1;  // alternate rows in a
      a[0][0]=1.;

      // compute k'th derivative
      for ( k=1; k<=n; k++ )
        {
        float d=0.0;
        int rk=r-k, pk=p-k;
        if ( r>=k )
          {
          a[s2][0]=a[s1][0]/ndu[pk+1][rk];
          d=a[s2][0]*ndu[rk][pk];
          }

        j1 = rk >= -1 ? 1 : -rk;
        j2 = (r-1<=pk) ? k-1 : p-r;

        for ( j=j1; j<=j2; j++ )
          {
          a[s2][j]=(a[s1][j]-a[s1][j-1])/ndu[pk+1][rk+j];
          d+=a[s2][j]*ndu[rk+j][pk];
          }
        if ( r<=pk )
          {
          a[s2][k]= -a[s1][k-1]/ndu[pk+1][r];
          d+=a[s2][k]*ndu[r][pk];
          }
        N[k][r]=d;
        j=s1; s1=s2; s2=j;  // switch rows
        }
      }

    r=p;
    for ( k=1; k<=n; k++ )
      {
      for ( j=0; j<=p; j++ )
        N[k][j]*=r;
      r*=p-k;
      }    
    }
}


// This method sets one of the control points
void DynamicBSpline2D::setControl(int i,int j,int d, float value) {
  // Change the value only if necessary
  if (P[d](i,j).x == value)
    return;

  // Set the value
  P[d](i,j).x = value;

  // Set the value in all adjacent matrices
  for (int k=0;k<4;k++)
    {
    for (int l=0;l<4;l++)
      {
      int iPatch = i-k;
      int jPatch = j-l;
      if (iPatch >= 0 && jPatch >= 0)
        {
        // Update the neighbour matrix
        SMLMatrix4f &nbr = P[d](iPatch,jPatch).nbr;
        nbr.put(l,k,value);

        P[d](iPatch,jPatch).mmNbrArray[(k << 2) + l] = value;
        P[d](iPatch,jPatch).mmNbrArrayT[(l << 2) + k] = value;

        // Update the MGM matrix
        // SMLMatrix4f T;
        // T.Multiply(nbr,MBspT);
        // P[d](iPatch,jPatch).MGM.Multiply(MBsp,T);
        // P[d](iPatch,jPatch).MGM.Transpose();

        // Update the timestamp on the patch
        if (iPatch < m[0]-2 && jPatch < m[1]-2)
          {
          tsPatch(iPatch,jPatch) = ++tsControl;
          }
        }
      }
    }
}

// Interpolate the spline or its derivatives over a grid with variable spacing.
// The out array must match the (dimensions of the grid)*outStride!!!
void DynamicBSpline2D::interpolateGridPoint(const SplineGridDefinition &grid,
                                            int u, int v,int uJet,int vJet,
                                            int iControlFirst,int iControlLast,
                                            float *out)

{
  SMLVec4f T;

  for (int k=iControlFirst;k<=iControlLast;k++)
    {
    // Find the control point
    BSplineControl &C = P[k](grid.knot(0,u)-3,grid.knot(1,v)-3);

    // This is the actual computation: a (1x4) * (4x4) * (4x1) matrix multiplication
    T = C.nbr * grid.W[0][u][uJet];

    // T.Transform(C.MGM,grid.W[0][u][uJet]);
    *out = dot_product(T,grid.W[1][v][vJet]); 

    // We just increment the pointer...
    out++;
    }
}

void DynamicBSpline2D::interpolatePoint(int iPatch,int jPatch,SMLVec4f *Wu,SMLVec4f *Wv,int uJet,int vJet,int iControlFirst,int iControlLast,float *out)
{
  SMLVec4f T;

  for (int k=iControlFirst;k<=iControlLast;k++)
    {
    // Find the control point
    BSplineControl &C = P[k](iPatch,jPatch);

    // This is the actual computation: a (1x4) * (4x4) * (4x1) matrix multiplication
    T = C.nbr * Wu[uJet];
    *out = dot_product(T,Wv[vJet]); 

    // We just increment the pointer...
    out++;
    }   
}




// Interpolate the spline or its derivatives over a grid with variable spacing.
// The out array must match the (dimensions of the grid)*outStride!!!
void DynamicBSpline2D::interpolateGrid(const SplineGridDefinition &grid,
                                       int uvFirst[2],int uvLast[2],int uvStep[2],int uJet,int vJet,
                                       int iControlFirst,int iControlLast,
                                       float *out)

{
  SMLVec4f T;

  int stride = 1 + iControlLast - iControlFirst;
  int iSkip = (stride)*(uvStep[0]-1);
  int jSkip = (stride)*(grid.size[0]*(uvStep[1]) - (uvStep[0]+uvLast[0]-uvFirst[0]));

  for (int j=uvFirst[1];j<=uvLast[1];j+=uvStep[1])
    {
    for (int i=uvFirst[0];i<=uvLast[0];i+=uvStep[0])
      {
      for (int k=iControlFirst;k<=iControlLast;k++)
        {
        // Find the control point
        BSplineControl &C = P[k](grid.knot(0,i)-3,grid.knot(1,j)-3);

        // This is the actual computation: a (1x4) * (4x4) * (4x1) matrix multiplication
        T = C.nbr * grid.W[0][i][uJet];
        *out = dot_product(T,grid.W[1][j][vJet]); 

        // We just increment the pointer...
        out++;
        }
      out += iSkip;
      }
    out += jSkip;
    }   
}

// This method creates a default m by n spline
DynamicBSpline2D::DynamicBSpline2D(int M,int N,int k,int zeroKnotMargin) {
  int d;

  // Store variables
  this->m[0] = M;
  this->m[1] = N;

  // Standard Spline Matrix 
  const float sixth = 1.0f/6.0f;
  float coeff[] = {-sixth, 0.5f, -0.5f, sixth, 0.5f, -1.0f, 0.5f, 0.0f, -0.5f, 0.0f, 0.5f, 0.0f, sixth, 4.0f*sixth, sixth, 0.0f};
  UB = (SMLMatrix4f *)_aligned_malloc(sizeof(SMLMatrix4f),16);
  *UB = SMLMatrix4f(coeff).transpose();

  // Initialize the jet multipliers
  for (d=0;d<2;d++)
    {
    float M = (float) m[d];
    jetMul[d] = (SMLVec4f*)_aligned_malloc(sizeof(SMLVec4f)*4,16);
    
    jetMul[d][0][0] = M*M*M;
    jetMul[d][0][1] = M*M;
    jetMul[d][0][2] = M;
    jetMul[d][0][3] = 1.0f;

    jetMul[d][1][0] = 3 * M*M*M;
    jetMul[d][1][1] = 2 * M*M;
    jetMul[d][1][2] = M;
    jetMul[d][1][3] = 0.0f;

    jetMul[d][2][0] = 6 * M*M*M;
    jetMul[d][2][1] = 2 * M*M;
    jetMul[d][2][2] = 0.0f;
    jetMul[d][2][3] = 0.0f;

    jetMul[d][3][0] = 6 * M*M*M;
    jetMul[d][3][1] = 0.0f;
    jetMul[d][3][2] = 0.0f;
    jetMul[d][3][3] = 0.0f;
    }

  // Allocate the control point array
  P = new BSplineControlArray[k];

  // Allocate the time stamp array
  tsControl = 0;
  tsPatch.resize(m[0]-2,m[1]-2);
  tsPatch.setAll(0);

  // Initialize the control point array
  for (int i=0;i<k;i++)
    {
    // Resise the array
    P[i].resize(m[0]+1,m[1]+1);

    // Set the values of the control points
    for (int x=0;x<=m[0];x++)
      {
      for (int y=0;y<=m[1];y++)
        {
        float value = 0;
        switch (i)
          {
          case 0 : value = ((float)x) / m[0];break;
          case 1 : value = ((float)y) / m[1];break;
          case 2 : value = 0.45f + (float)(0.1 * rand()) / RAND_MAX;break;
          case 3 : 
            if (x==0||y==0||x==m[0]||y==m[1])
              value = -0.7f;
            // value = 0.025 + (float)(0.05 * rand()) / RAND_MAX;
            else
              value = 0.025f + (float)(0.05 * rand()) / RAND_MAX;
            0.075;
            break;
          }

        setControl(x,y,i,value);
        }
      }       
    }

  // Loop over the two dimensions
  for (d=0;d<2;d++)
    {
    // Create knot arrays
    knots[d] = new float[dim(d)+5];
    float scale = 1.0f / dim(d);

    // Set the knot sequences
    for (int i=0;i<=dim(d)+4;i++)
      {
      if (i<zeroKnotMargin)
        knots[d][i] = 0.0f;
      else if (i > dim(d)+4-zeroKnotMargin)
        knots[d][i] = 1.0f;
      else
        knots[d][i] = (i-zeroKnotMargin) * scale;
      }

    regKnots = (zeroKnotMargin==2);
    }
}

DynamicBSpline2D::~DynamicBSpline2D() {
  for (int d=0;d<2;d++)
    {
    delete[] knots[d];
    }
  delete[] P;

  _aligned_free(UB);
  _aligned_free(jetMul[0]);
  _aligned_free(jetMul[1]);
}


void DynamicBSpline2D::push(int uKnot,int vKnot,float u,float v,int iControlFirst,int iControlLast,float *vec) {
  // Get the nearest knots' u and v values
  float uk = knots[0][uKnot];
  float vk = knots[1][vKnot];
  float uk1 = knots[0][uKnot+1];
  float vk1 = knots[1][vKnot+1];

  int i;

  // Compute the basis function at the position
  ALIGN_PRE SMLVec4f ALIGN_POST NJu[4],NJv[4];
  basisJet(0,uKnot,u,NJu);
  basisJet(1,vKnot,v,NJv);

  SMLVec4f &Nu = NJu[0],&Nv = NJv[0];
  //SMLVec4f Nu,Nv;
  //MBspT.Transform(NJu[0],Nu);
  //MBspT.Transform(NJv[0],Nv);

  // Get the length of the push vector
  float len = 0;
  for (i=0;i<=iControlLast-iControlFirst;i++)
    len += vec[i]*vec[i];
  len = sqrt(len);

  // Normalize the vector
  if (len > 0)
    {
    for (i=0;i<=iControlLast-iControlFirst;i++)
      vec[i] /= len;

    double g = (u - uk) / (uk1-uk);
    double d = (v - vk) / (vk1-vk);
    double a = len / (((1-g)*(1-d)* Nu[1]*Nv[1]) + ((g)*(1-d)* Nu[2]*Nv[1]) 
                      + ((1-g)*(d)* Nu[1]*Nv[2]) + ((d)*(g)* Nu[2]*Nv[2]));

    // Update the control points
    for (i=0;i<2;i++)
      {
      for (int j=0;j<2;j++)
        {
        for (int k=iControlFirst;k<=iControlLast;k++)
          {
          float p = getControl(uKnot-2+i,vKnot-2+j,k);
          p += (float)
            ((i?g:1-g) * (j?d:1-d) * a * vec[k-iControlFirst] * Nu.data_block()[i+1] * Nv.data_block()[j+1]);
          setControl(uKnot-2+i,vKnot-2+j,k,p);
          }
        }
      }
    }
}


void DynamicBSpline2D::getControlArray(Array2D<SMLVec3f> &out,int startIdx) {
  out.resize(P[0].width(),P[0].height());
  for (int i=0;i<P[0].width();i++)
    {
    for (int j=0;j<P[0].height();j++)
      {
      out(i,j)[0] = P[startIdx](i,j).x;
      out(i,j)[1] = P[startIdx+1](i,j).x;
      out(i,j)[2] = P[startIdx+2](i,j).x;
      }
    }
}

void DynamicBSpline2D::getControlArray(Array2D<SMLVec4f> &out,int startIdx) {
  out.resize(P[0].width(),P[0].height());
  for (int i=0;i<P[0].width();i++)
    {
    for (int j=0;j<P[0].height();j++)
      {
      out(i,j)[0] = P[startIdx](i,j).x;
      out(i,j)[1] = P[startIdx+1](i,j).x;
      out(i,j)[2] = P[startIdx+2](i,j).x;
      out(i,j)[3] = P[startIdx+3](i,j).x;
      }
    }
}

/*****************************************************************
 *  MSpline
 *****************************************************************/

void testMultiply() {
  float xx[] = {1,2,3,4};
  float mm[] = {1,0,3,0,0,1,2,0,4,0,1,2,0,4,0,1};
  SMLMatrix4f M(mm);
  SMLVec4f X(xx),R,R1;

  __m128 A0= _mm_loadu_ps(M.data_block());
  __m128 A1= _mm_loadu_ps(M.data_block()+4);
  __m128 A2= _mm_loadu_ps(M.data_block()+8);
  __m128 A3= _mm_loadu_ps(M.data_block()+12);

  __m128 B = _mm_loadu_ps(X.data_block());

  __m128 TM,B0,B1,B2,B3;

  // Perform the vector multiplication
  A0 = _mm_mul_ps(A0,B);
  A1 = _mm_mul_ps(A1,B);
  A2 = _mm_mul_ps(A2,B);
  A3 = _mm_mul_ps(A3,B);

  B1 = _mm_shuffle_ps(A0,A1,0x44);
  TM = _mm_shuffle_ps(A2,A3,0x44);
  B0 = _mm_shuffle_ps(B1,TM,0x88);
  B1 = _mm_shuffle_ps(B1,TM,0xDD);

  B3 = _mm_shuffle_ps(A0,A1,0xee);
  TM = _mm_shuffle_ps(A2,A3,0xee);
  B2 = _mm_shuffle_ps(B3,TM,0x88);
  B3 = _mm_shuffle_ps(B3,TM,0xDD);

  // Perform the additions
  B0 = _mm_add_ps(B0,B1);
  B2 = _mm_add_ps(B2,B3);

  // Perform the final addition
  B0 = _mm_add_ps(B0,B2);

  _mm_store_ps(R.data_block(),B0);

  R1 = M * X;

  //_mm_storeu_ps(dat,B0);
  cout << R[0] << ", "  << R[1] << ", "  << R[2] << ", "  << R[3] << endl;
  cout << R1[0] << ", "  << R1[1] << ", "  << R1[2] << ", "  << R1[3] << endl;
}

// Perform the cross product 
inline __m128 crossProduct(__m128 A,__m128 B) {
  __m128 A1,B1;

  // Shuffle the elements
  A1 = _mm_shuffle_ps(A,A,0xd2);
  B1 = _mm_shuffle_ps(B,B,0xc9);
  A  = _mm_shuffle_ps(A,A,0xc9);
  B  = _mm_shuffle_ps(B,B,0xd2);

  // Multiply pairwise
  A = _mm_mul_ps(A,B);
  A1 = _mm_mul_ps(A1,B1);

  // Add the elements
  return _mm_sub_ps(A,A1);
}

// Perform the cross product 
inline void triangleAreaVec(float *A,float *B,float *C, float *outR) {
  __m128 r0,r1,r2;
  ALIGN_PRE float N[4] ALIGN_POST;

  // Load A, B and C
  r0 = _mm_loadu_ps(A);
  r1 = _mm_loadu_ps(C);
  r2 = _mm_loadu_ps(B);

  // Compute differences
  r0 = _mm_sub_ps(r0,r2);
  r1 = _mm_sub_ps(r1,r2);

  // Compute cross product
  r0 = crossProduct(r0,r1);

  // Pull out the area vector
  _mm_store_ps(N, r0);

  // Copy into the actual N
  outR[0] = N[0];
  outR[1] = N[1];
  outR[2] = N[2];
}

inline void triangleArea(float *A,float *B,float *C,float *outN, float *outArea) 
{
  __m128 r0,r1,r2;
  ALIGN_PRE float N[4] ALIGN_POST;

  // Load A, B and C
  r0 = _mm_loadu_ps(A);
  r1 = _mm_loadu_ps(C);
  r2 = _mm_loadu_ps(B);

  // Compute differences
  r0 = _mm_sub_ps(r0,r2);
  r1 = _mm_sub_ps(r1,r2);

  // Compute cross product
  r0 = crossProduct(r0,r1);

  // Pull out the area vector
  _mm_store_ps(N, r0);

  // Copy into the actual N
  outN[0] = N[0];
  outN[1] = N[1];
  outN[2] = N[2];
  
  // outN[0] = r0[0];
  // outN[1] = r0[1];
  // outN[2] = r0[2];

  // Square the elements
  r0 = _mm_mul_ps(r0,r0);

  // Add the elements
  r1 = _mm_shuffle_ps(r0,r0,0x01);
  r2 = _mm_shuffle_ps(r0,r0,0x02);
  r0 = _mm_add_ps(r0,r1);
  r0 = _mm_add_ps(r0,r2);

  // Take the fast square root
  r0 = _mm_sqrt_ss(r0);

  // We have the area vector for the triangle
  _mm_store_ss(outArea,r0);
}

// Normalize a vector
inline __m128 normalize(__m128 A) {
  __m128 B,C;

  // Square the vector
  B = _mm_mul_ps(A,A);

  // Add the elements
  C = _mm_shuffle_ps(B,B,0xb1);
  B = _mm_add_ps(B,C);
  C = _mm_shuffle_ps(B,B,0x0a);
  B = _mm_add_ps(B,C);

  // Take square root recipient of B
  B = _mm_rsqrt_ps(B);

  // Scale the vector A
  A = _mm_mul_ps(A,B);

  return A;
}

/*
void MSpline::interpolateMedialSurfacePatch(const SplineGridDefinition &grid,int iPatch,int jPatch,int step,
                                                          Aux2D<SMLVec4f> &X,Aux2D<SMLVec4f> &Xu,Aux2D<SMLVec4f> &Xv,
                                                          Aux2D<SMLVec4f> &NRaw,Aux2D<SMLVec4f> &N)
{   
    // Some vectors to use
    SMLVec4f **Wu = grid.patchWeight[0][iPatch];
    SMLVec4f **Wv = grid.patchWeight[1][jPatch];
    
    // Load the matrices in advance
    SMLMatrix4f Cx = P[0](iPatch,jPatch).nbr;
    SMLMatrix4f Cy = P[1](iPatch,jPatch).nbr;
    SMLMatrix4f Cz = P[2](iPatch,jPatch).nbr;
    SMLMatrix4f Cr = P[3](iPatch,jPatch).nbr;

    Cx.Transpose();
    Cy.Transpose();
    Cz.Transpose();
    Cr.Transpose();

    // Rows from these matrices
    __m128 CX0 = _mm_loadu_ps(Cx.GetData());
    __m128 CX1 = _mm_loadu_ps(Cx.GetData()+4);
    __m128 CX2 = _mm_loadu_ps(Cx.GetData()+8);
    __m128 CX3 = _mm_loadu_ps(Cx.GetData()+12);

    __m128 CY0 = _mm_loadu_ps(Cy.GetData());
    __m128 CY1 = _mm_loadu_ps(Cy.GetData()+4);
    __m128 CY2 = _mm_loadu_ps(Cy.GetData()+8);
    __m128 CY3 = _mm_loadu_ps(Cy.GetData()+12);
    
    __m128 CZ0 = _mm_loadu_ps(Cz.GetData());
    __m128 CZ1 = _mm_loadu_ps(Cz.GetData()+4);
    __m128 CZ2 = _mm_loadu_ps(Cz.GetData()+8);
    __m128 CZ3 = _mm_loadu_ps(Cz.GetData()+12);
    
    __m128 CR0 = _mm_loadu_ps(Cr.GetData());
    __m128 CR1 = _mm_loadu_ps(Cr.GetData()+4);
    __m128 CR2 = _mm_loadu_ps(Cr.GetData()+8);
    __m128 CR3 = _mm_loadu_ps(Cr.GetData()+12);
    
    // Starting and ending positions
    int jPatchStart = grid.patchStart(1,jPatch);
    int iPatchStart = grid.patchStart(0,iPatch);
    int jPatchEnd = grid.patchStart(1,jPatch+1);
    int iPatchEnd = grid.patchStart(0,iPatch+1);
    int jLast = jPatchEnd - jPatchStart;
    int iLast = iPatchEnd - iPatchStart;
    int i,j;

    for(j=0;j<=jLast;j+=step) {
        
        // Get the weight factors
        __m128 WV0 = _mm_loadu_ps(Wv[j][0].data());
        __m128 WV1 = _mm_loadu_ps(Wv[j][1].data());

        // Multiply this by the matrices that we need
        __m128 LX0 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WV0);
        __m128 LX1 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WV1);
        __m128 LY0 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WV0);
        __m128 LY1 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WV1);
        __m128 LZ0 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WV0);
        __m128 LZ1 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WV1);
        __m128 LR0 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WV0);
        __m128 LR1 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WV1);

        // Get pointers to data rows for extra speed (maybe)
        SMLVec4f *rowX  = X.row(j);
        SMLVec4f *rowXu = Xu.row(j);
        SMLVec4f *rowXv = Xv.row(j);
        SMLVec4f *rowN = N.row(j);
        SMLVec4f *rowNraw = NRaw.row(j);

        for(i=0;i<=iLast;i+=step) {
            // Get the U, V data
            __m128 WU0 = _mm_loadu_ps(Wu[i][0].data());
            __m128 WU1 = _mm_loadu_ps(Wu[i][1].data());

            __m128 mX =  mul_4x4_4x1(LX0,LY0,LZ0,LR0,WU0);
            __m128 mXu = mul_4x4_4x1(LX0,LY0,LZ0,LR0,WU1);
            __m128 mXv = mul_4x4_4x1(LX1,LY1,LZ1,LR1,WU0);

            // Now the unit normal.
            __m128 mNRaw = crossProduct(mXu,mXv);
            __m128 mN = normalize(mNRaw);

            _mm_store_ps(rowX[i].data(),mX);
            _mm_store_ps(rowXu[i].data(),mXu);
            _mm_store_ps(rowXv[i].data(),mXv);
            _mm_store_ps(rowNraw[i].data(),mNRaw);
            _mm_store_ps(rowN[i].data(),mN);            
        }
    }
}*/

/*
void MSpline::interpolateMedialSurfacePatch(const SplineGridDefinition &grid,int iPatch,int jPatch,int step,
                                                          Aux2D<SMLVec4f> &X,Aux2D<SMLVec4f> &Xu,Aux2D<SMLVec4f> &Xv,
                                                          Aux2D<SMLVec4f> &Xuu,Aux2D<SMLVec4f> &Xuv,Aux2D<SMLVec4f> &Xvv,
                                                          Aux2D<SMLVec4f> &NRaw,Aux2D<SMLVec4f> &N)
{   
    // Some vectors to use
    SMLVec4f **Wu = grid.patchWeight[0][iPatch];
    SMLVec4f **Wv = grid.patchWeight[1][jPatch];
    
    // Load the matrices in advance
    SMLMatrix4f Cx = P[0](iPatch,jPatch).nbr;
    SMLMatrix4f Cy = P[1](iPatch,jPatch).nbr;
    SMLMatrix4f Cz = P[2](iPatch,jPatch).nbr;
    SMLMatrix4f Cr = P[3](iPatch,jPatch).nbr;

    Cx.Transpose();
    Cy.Transpose();
    Cz.Transpose();
    Cr.Transpose();

    // Rows from these matrices
    __m128 CX0 = _mm_loadu_ps(Cx.GetData());
    __m128 CX1 = _mm_loadu_ps(Cx.GetData()+4);
    __m128 CX2 = _mm_loadu_ps(Cx.GetData()+8);
    __m128 CX3 = _mm_loadu_ps(Cx.GetData()+12);

    __m128 CY0 = _mm_loadu_ps(Cy.GetData());
    __m128 CY1 = _mm_loadu_ps(Cy.GetData()+4);
    __m128 CY2 = _mm_loadu_ps(Cy.GetData()+8);
    __m128 CY3 = _mm_loadu_ps(Cy.GetData()+12);
    
    __m128 CZ0 = _mm_loadu_ps(Cz.GetData());
    __m128 CZ1 = _mm_loadu_ps(Cz.GetData()+4);
    __m128 CZ2 = _mm_loadu_ps(Cz.GetData()+8);
    __m128 CZ3 = _mm_loadu_ps(Cz.GetData()+12);
    
    __m128 CR0 = _mm_loadu_ps(Cr.GetData());
    __m128 CR1 = _mm_loadu_ps(Cr.GetData()+4);
    __m128 CR2 = _mm_loadu_ps(Cr.GetData()+8);
    __m128 CR3 = _mm_loadu_ps(Cr.GetData()+12);
    
    // Starting and ending positions
    int jPatchStart = grid.patchStart(1,jPatch);
    int iPatchStart = grid.patchStart(0,iPatch);
    int jPatchEnd = grid.patchStart(1,jPatch+1);
    int iPatchEnd = grid.patchStart(0,iPatch+1);
    int jLast = jPatchEnd - jPatchStart;
    int iLast = iPatchEnd - iPatchStart;
    int i,j;

    for(j=0;j<=jLast;j+=step) {
        
        // Get the weight factors
        __m128 WV0 = _mm_loadu_ps(Wv[j][0].data());
        __m128 WV1 = _mm_loadu_ps(Wv[j][1].data());
        __m128 WV2 = _mm_loadu_ps(Wv[j][2].data());

        // Multiply this by the matrices that we need
        __m128 LX0 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WV0);
        __m128 LY0 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WV0);
        __m128 LZ0 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WV0);
        __m128 LR0 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WV0);

        __m128 LX1 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WV1);
        __m128 LY1 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WV1);
        __m128 LZ1 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WV1);
        __m128 LR1 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WV1);

        __m128 LX2 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WV2);
        __m128 LY2 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WV2);
        __m128 LZ2 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WV2);
        __m128 LR2 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WV2);

        // Get pointers to data rows for extra speed (maybe)
        SMLVec4f *rowX  = X.row(j);
        SMLVec4f *rowXu = Xu.row(j);
        SMLVec4f *rowXv = Xv.row(j);

        SMLVec4f *rowXuu = Xuu.row(j);
        SMLVec4f *rowXuv = Xuv.row(j);
        SMLVec4f *rowXvv = Xvv.row(j);

        SMLVec4f *rowN = N.row(j);
        SMLVec4f *rowNraw = NRaw.row(j);

        for(i=0;i<=iLast;i+=step) {
            // Get the U, V data
            __m128 WU0 = _mm_loadu_ps(Wu[i][0].data());
            __m128 WU1 = _mm_loadu_ps(Wu[i][1].data());
            __m128 WU2 = _mm_loadu_ps(Wu[i][2].data());

            __m128 mX =  mul_4x4_4x1(LX0,LY0,LZ0,LR0,WU0);
            __m128 mXu = mul_4x4_4x1(LX0,LY0,LZ0,LR0,WU1);
            __m128 mXv = mul_4x4_4x1(LX1,LY1,LZ1,LR1,WU0);
            
            __m128 mXuu = mul_4x4_4x1(LX0,LY0,LZ0,LR0,WU2);
            __m128 mXuv = mul_4x4_4x1(LX1,LY1,LZ1,LR1,WU1);
            __m128 mXvv = mul_4x4_4x1(LX2,LY2,LZ2,LR2,WU0);

            // Now the unit normal.
            __m128 mNRaw = crossProduct(mXu,mXv);
            __m128 mN = normalize(mNRaw);

            _mm_store_ps(rowX[i].data(),mX);
            _mm_store_ps(rowXu[i].data(),mXu);
            _mm_store_ps(rowXv[i].data(),mXv);

            _mm_store_ps(rowXuu[i].data(),mXuu);
            _mm_store_ps(rowXuv[i].data(),mXuv);
            _mm_store_ps(rowXvv[i].data(),mXvv);

            _mm_store_ps(rowNraw[i].data(),mNRaw);
            _mm_store_ps(rowN[i].data(),mN);            
        }
    }
}
*/

void MSpline::interpolateMedialPoint02(int iPatch,int jPatch,
                                       SMLVec4f *Wu,SMLVec4f *Wv,MedialPoint &mp) 
{
  // Load the matrices in advance
  float *Cx = P[0](iPatch,jPatch).mmNbrArray;
  float *Cy = P[1](iPatch,jPatch).mmNbrArray;
  float *Cz = P[2](iPatch,jPatch).mmNbrArray;
  float *Cr = P[3](iPatch,jPatch).mmNbrArray;

  // Rows from these matrices
  __m128 CX0 = _mm_load_ps(Cx);
  __m128 CX1 = _mm_load_ps(Cx+4);
  __m128 CX2 = _mm_load_ps(Cx+8);
  __m128 CX3 = _mm_load_ps(Cx+12);

  __m128 CY0 = _mm_load_ps(Cy);
  __m128 CY1 = _mm_load_ps(Cy+4);
  __m128 CY2 = _mm_load_ps(Cy+8);
  __m128 CY3 = _mm_load_ps(Cy+12);

  __m128 CZ0 = _mm_load_ps(Cz);
  __m128 CZ1 = _mm_load_ps(Cz+4);
  __m128 CZ2 = _mm_load_ps(Cz+8);
  __m128 CZ3 = _mm_load_ps(Cz+12);

  __m128 CR0 = _mm_load_ps(Cr);
  __m128 CR1 = _mm_load_ps(Cr+4);
  __m128 CR2 = _mm_load_ps(Cr+8);
  __m128 CR3 = _mm_load_ps(Cr+12);

  // Get the weight factors
  __m128 WV0 = _mm_load_ps(Wv[0].data_block());
  __m128 WV1 = _mm_load_ps(Wv[1].data_block());
  __m128 WV2 = _mm_load_ps(Wv[2].data_block());

  // Multiply this by the matrices that we need
  __m128 LX0 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WV0);
  __m128 LY0 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WV0);
  __m128 LZ0 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WV0);
  __m128 LR0 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WV0);

  __m128 LX1 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WV1);
  __m128 LY1 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WV1);
  __m128 LZ1 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WV1);
  __m128 LR1 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WV1);

  __m128 LX2 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WV2);
  __m128 LY2 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WV2);
  __m128 LZ2 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WV2);
  __m128 LR2 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WV2);

  // Get the U, V data
  __m128 WU0 = _mm_load_ps(Wu[0].data_block());
  __m128 WU1 = _mm_load_ps(Wu[1].data_block());
  __m128 WU2 = _mm_load_ps(Wu[2].data_block());

  __m128 mX =  mul_4x4_4x1(LX0,LY0,LZ0,LR0,WU0);
  __m128 mXu = mul_4x4_4x1(LX0,LY0,LZ0,LR0,WU1);
  __m128 mXv = mul_4x4_4x1(LX1,LY1,LZ1,LR1,WU0);

  __m128 mXuu = mul_4x4_4x1(LX0,LY0,LZ0,LR0,WU2);
  __m128 mXuv = mul_4x4_4x1(LX1,LY1,LZ1,LR1,WU1);
  __m128 mXvv = mul_4x4_4x1(LX2,LY2,LZ2,LR2,WU0);

  // Now the unit normal.
  __m128 mNRaw = crossProduct(mXu,mXv);
  __m128 mN = normalize(mNRaw);

  // Access the medial point
  _mm_store_ps(mp.F.data_block(),mX);
  _mm_store_ps(mp.Fu.data_block(),mXu);
  _mm_store_ps(mp.Fv.data_block(),mXv);

  _mm_store_ps(mp.Fuu.data_block(),mXuu);
  _mm_store_ps(mp.Fuv.data_block(),mXuv);
  _mm_store_ps(mp.Fvv.data_block(),mXvv);

  _mm_store_ps(mp.F3raw.data_block(),mNRaw);
  _mm_store_ps(mp.F3.data_block(),mN);         
}

void MSpline::interpolateMedialCrestD1(int iPatch,int jPatch,SMLVec4f *Wu,SMLVec4f *Wv,MedialPoint &mp) 
{
  // Load the matrices in advance
  float *Cx = P[0](iPatch,jPatch).mmNbrArray;
  float *Cy = P[1](iPatch,jPatch).mmNbrArray;
  float *Cz = P[2](iPatch,jPatch).mmNbrArray;
  float *Cr = P[3](iPatch,jPatch).mmNbrArray;

  // Rows from these matrices
  __m128 CX0 = _mm_load_ps(Cx);
  __m128 CX1 = _mm_load_ps(Cx+4);
  __m128 CX2 = _mm_load_ps(Cx+8);
  __m128 CX3 = _mm_load_ps(Cx+12);
  __m128 CY0 = _mm_load_ps(Cy);
  __m128 CY1 = _mm_load_ps(Cy+4);
  __m128 CY2 = _mm_load_ps(Cy+8);
  __m128 CY3 = _mm_load_ps(Cy+12);
  __m128 CZ0 = _mm_load_ps(Cz);
  __m128 CZ1 = _mm_load_ps(Cz+4);
  __m128 CZ2 = _mm_load_ps(Cz+8);
  __m128 CZ3 = _mm_load_ps(Cz+12);
  __m128 CR0 = _mm_load_ps(Cr);
  __m128 CR1 = _mm_load_ps(Cr+4);
  __m128 CR2 = _mm_load_ps(Cr+8);
  __m128 CR3 = _mm_load_ps(Cr+12);

  // Get the weight factors
  __m128 WU0 = _mm_load_ps(Wu[0].data_block());
  __m128 WU1 = _mm_load_ps(Wu[1].data_block());
  __m128 WU2 = _mm_load_ps(Wu[2].data_block());
  __m128 WV0 = _mm_load_ps(Wv[0].data_block());
  __m128 WV1 = _mm_load_ps(Wv[1].data_block());

  // Multiply this by the matrices that we need
  __m128 LX0 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WV0);
  __m128 LY0 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WV0);
  __m128 LZ0 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WV0);
  __m128 LR0 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WV0);

  __m128 LX1 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WV1);
  __m128 LY1 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WV1);
  __m128 LZ1 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WV1);
  __m128 LR1 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WV1);

  // Get the U, V data
  __m128 mXu = mul_4x4_4x1(LX0,LY0,LZ0,LR0,WU1);
  __m128 mXv = mul_4x4_4x1(LX1,LY1,LZ1,LR1,WU0);
  __m128 mXuu = mul_4x4_4x1(LX0,LY0,LZ0,LR0,WU2);
  __m128 mXuv = mul_4x4_4x1(LX1,LY1,LZ1,LR1,WU1);

  // Now the unit normal.
  __m128 mNRaw = crossProduct(mXu,mXv);

  // Access the medial point
  _mm_store_ps(mp.Fu.data_block(),mXu);
  _mm_store_ps(mp.Fv.data_block(),mXv);
  _mm_store_ps(mp.Fuu.data_block(),mXuu);
  _mm_store_ps(mp.Fuv.data_block(),mXuv);
  _mm_store_ps(mp.F3raw.data_block(),mNRaw);
}


void MSpline::interpolateMedialCrestD2(int iPatch,int jPatch,SMLVec4f *Wu,SMLVec4f *Wv,MedialPoint &mp) 
{
  // Load the matrices in advance
  float *Cx = P[0](iPatch,jPatch).mmNbrArrayT;
  float *Cy = P[1](iPatch,jPatch).mmNbrArrayT;
  float *Cz = P[2](iPatch,jPatch).mmNbrArrayT;
  float *Cr = P[3](iPatch,jPatch).mmNbrArrayT;

  // Rows from these matrices
  __m128 CX0 = _mm_load_ps(Cx);
  __m128 CX1 = _mm_load_ps(Cx+4);
  __m128 CX2 = _mm_load_ps(Cx+8);
  __m128 CX3 = _mm_load_ps(Cx+12);
  __m128 CY0 = _mm_load_ps(Cy);
  __m128 CY1 = _mm_load_ps(Cy+4);
  __m128 CY2 = _mm_load_ps(Cy+8);
  __m128 CY3 = _mm_load_ps(Cy+12);
  __m128 CZ0 = _mm_load_ps(Cz);
  __m128 CZ1 = _mm_load_ps(Cz+4);
  __m128 CZ2 = _mm_load_ps(Cz+8);
  __m128 CZ3 = _mm_load_ps(Cz+12);
  __m128 CR0 = _mm_load_ps(Cr);
  __m128 CR1 = _mm_load_ps(Cr+4);
  __m128 CR2 = _mm_load_ps(Cr+8);
  __m128 CR3 = _mm_load_ps(Cr+12);

  // Get the weight factors
  __m128 WU0 = _mm_load_ps(Wu[0].data_block());
  __m128 WU1 = _mm_load_ps(Wu[1].data_block());
  __m128 WV0 = _mm_load_ps(Wv[0].data_block());
  __m128 WV1 = _mm_load_ps(Wv[1].data_block());
  __m128 WV2 = _mm_load_ps(Wv[2].data_block());

  // Multiply this by the matrices that we need
  __m128 LX0 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WU0);
  __m128 LY0 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WU0);
  __m128 LZ0 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WU0);
  __m128 LR0 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WU0);

  __m128 LX1 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WU1);
  __m128 LY1 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WU1);
  __m128 LZ1 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WU1);
  __m128 LR1 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WU1);

  // Get the U, V data
  __m128 mXv = mul_4x4_4x1(LX0,LY0,LZ0,LR0,WV1);
  __m128 mXu = mul_4x4_4x1(LX1,LY1,LZ1,LR1,WV0);
  __m128 mXuv = mul_4x4_4x1(LX1,LY1,LZ1,LR1,WV1);
  __m128 mXvv = mul_4x4_4x1(LX0,LY0,LZ0,LR0,WV2);

  // Now the unit normal.
  __m128 mNRaw = crossProduct(mXu,mXv);

  // Access the medial point
  _mm_store_ps(mp.Fu.data_block(),mXu);
  _mm_store_ps(mp.Fv.data_block(),mXv);
  _mm_store_ps(mp.Fuv.data_block(),mXuv);
  _mm_store_ps(mp.Fvv.data_block(),mXvv);
  _mm_store_ps(mp.F3raw.data_block(),mNRaw);
}

void MSpline::interpolateMedialCrestD1Missing(int iPatch,int jPatch,
                                              SMLVec4f *Wu,SMLVec4f *Wv,MedialPoint &mp) 
{
  // Load the matrices in advance
  float *Cx = P[0](iPatch,jPatch).mmNbrArrayT;
  float *Cy = P[1](iPatch,jPatch).mmNbrArrayT;
  float *Cz = P[2](iPatch,jPatch).mmNbrArrayT;
  float *Cr = P[3](iPatch,jPatch).mmNbrArrayT;

  // Rows from these matrices
  __m128 CX0 = _mm_load_ps(Cx);
  __m128 CX1 = _mm_load_ps(Cx+4);
  __m128 CX2 = _mm_load_ps(Cx+8);
  __m128 CX3 = _mm_load_ps(Cx+12);

  __m128 CY0 = _mm_load_ps(Cy);
  __m128 CY1 = _mm_load_ps(Cy+4);
  __m128 CY2 = _mm_load_ps(Cy+8);
  __m128 CY3 = _mm_load_ps(Cy+12);

  __m128 CZ0 = _mm_load_ps(Cz);
  __m128 CZ1 = _mm_load_ps(Cz+4);
  __m128 CZ2 = _mm_load_ps(Cz+8);
  __m128 CZ3 = _mm_load_ps(Cz+12);

  __m128 CR0 = _mm_load_ps(Cr);
  __m128 CR1 = _mm_load_ps(Cr+4);
  __m128 CR2 = _mm_load_ps(Cr+8);
  __m128 CR3 = _mm_load_ps(Cr+12);

  // Get the weight factors
  __m128 WV0 = _mm_load_ps(Wv[0].data_block());
  __m128 WV2 = _mm_load_ps(Wv[2].data_block());
  __m128 WU0 = _mm_load_ps(Wu[0].data_block());

  // Multiply this by the matrices that we need
  __m128 LX0 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WU0);
  __m128 LY0 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WU0);
  __m128 LZ0 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WU0);
  __m128 LR0 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WU0);

  // Get the U, V data
  __m128 mX =  mul_4x4_4x1(LX0,LY0,LZ0,LR0,WV0);
  __m128 mXvv = mul_4x4_4x1(LX0,LY0,LZ0,LR0,WV2);

  // Now the unit normal.
  __m128 mNRaw = _mm_load_ps(mp.F3raw.data_block());
  __m128 mN = normalize(mNRaw);

  // Access the medial point
  _mm_store_ps(mp.F.data_block(),mX);
  _mm_store_ps(mp.Fvv.data_block(),mXvv);
  _mm_store_ps(mp.F3.data_block(),mN);         
}


void MSpline::interpolateMedialCrestD2Missing(int iPatch,int jPatch,
                                              SMLVec4f *Wu,SMLVec4f *Wv,MedialPoint &mp) 
{
  // Load the matrices in advance
  float *Cx = P[0](iPatch,jPatch).mmNbrArray;
  float *Cy = P[1](iPatch,jPatch).mmNbrArray;
  float *Cz = P[2](iPatch,jPatch).mmNbrArray;
  float *Cr = P[3](iPatch,jPatch).mmNbrArray;

  // Rows from these matrices
  __m128 CX0 = _mm_load_ps(Cx);
  __m128 CX1 = _mm_load_ps(Cx+4);
  __m128 CX2 = _mm_load_ps(Cx+8);
  __m128 CX3 = _mm_load_ps(Cx+12);

  __m128 CY0 = _mm_load_ps(Cy);
  __m128 CY1 = _mm_load_ps(Cy+4);
  __m128 CY2 = _mm_load_ps(Cy+8);
  __m128 CY3 = _mm_load_ps(Cy+12);

  __m128 CZ0 = _mm_load_ps(Cz);
  __m128 CZ1 = _mm_load_ps(Cz+4);
  __m128 CZ2 = _mm_load_ps(Cz+8);
  __m128 CZ3 = _mm_load_ps(Cz+12);

  __m128 CR0 = _mm_load_ps(Cr);
  __m128 CR1 = _mm_load_ps(Cr+4);
  __m128 CR2 = _mm_load_ps(Cr+8);
  __m128 CR3 = _mm_load_ps(Cr+12);

  // Get the weight factors
  __m128 WV0 = _mm_load_ps(Wv[0].data_block());
  __m128 WU0 = _mm_load_ps(Wu[0].data_block());
  __m128 WU2 = _mm_load_ps(Wu[2].data_block());

  // Multiply this by the matrices that we need
  __m128 LX0 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WV0);
  __m128 LY0 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WV0);
  __m128 LZ0 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WV0);
  __m128 LR0 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WV0);

  // Get the U, V data
  __m128 mX   = mul_4x4_4x1(LX0,LY0,LZ0,LR0,WU0);
  __m128 mXuu = mul_4x4_4x1(LX0,LY0,LZ0,LR0,WU2);

  // Now the unit normal.
  __m128 mNRaw = _mm_load_ps(mp.F3raw.data_block());
  __m128 mN = normalize(mNRaw);

  // Access the medial point
  _mm_store_ps(mp.F.data_block(),mX);
  _mm_store_ps(mp.Fuu.data_block(),mXuu);
  _mm_store_ps(mp.F3.data_block(),mN);         
}

void MSpline::interpolateMedialPoint02(const SplineGridDefinition &grid,int u,int v,MedialPoint &mp) {
  int iPatch = grid.knot(0,u)-3;
  int jPatch = grid.knot(1,v)-3;
  int i = u - grid.patchStart(0,iPatch);
  int j = v - grid.patchStart(1,jPatch);
  interpolateMedialPoint02(iPatch,jPatch,grid.patchWeight[0][iPatch][i],grid.patchWeight[1][jPatch][j],mp);
}

void MSpline::interpolateMedialSurfacePatch(const SplineGridDefinition &grid,
                                            int iPatch,int jPatch,
                                            int step, Aux2D<MedialPoint> *MP)
{
  // Some vectors to use
  SMLVec4f **Wu = grid.patchWeight[0][iPatch];
  SMLVec4f **Wv = grid.patchWeight[1][jPatch];

  // Load the matrices in advance
  float *Cx = P[0](iPatch,jPatch).mmNbrArray;
  float *Cy = P[1](iPatch,jPatch).mmNbrArray;
  float *Cz = P[2](iPatch,jPatch).mmNbrArray;
  float *Cr = P[3](iPatch,jPatch).mmNbrArray;

  // Rows from these matrices
  __m128 CX0 = _mm_load_ps(Cx);
  __m128 CX1 = _mm_load_ps(Cx+4);
  __m128 CX2 = _mm_load_ps(Cx+8);
  __m128 CX3 = _mm_load_ps(Cx+12);

  __m128 CY0 = _mm_load_ps(Cy);
  __m128 CY1 = _mm_load_ps(Cy+4);
  __m128 CY2 = _mm_load_ps(Cy+8);
  __m128 CY3 = _mm_load_ps(Cy+12);

  __m128 CZ0 = _mm_load_ps(Cz);
  __m128 CZ1 = _mm_load_ps(Cz+4);
  __m128 CZ2 = _mm_load_ps(Cz+8);
  __m128 CZ3 = _mm_load_ps(Cz+12);

  __m128 CR0 = _mm_load_ps(Cr);
  __m128 CR1 = _mm_load_ps(Cr+4);
  __m128 CR2 = _mm_load_ps(Cr+8);
  __m128 CR3 = _mm_load_ps(Cr+12);

  // Starting and ending positions
  int jPatchStart = grid.patchStart(1,jPatch);
  int iPatchStart = grid.patchStart(0,iPatch);
  int jPatchEnd = grid.patchStart(1,jPatch+1);
  int iPatchEnd = grid.patchStart(0,iPatch+1);
  int jLast = jPatchEnd - jPatchStart;
  int iLast = iPatchEnd - iPatchStart;
  int i,j;

  for (j=0;j<=jLast;j+=step)
    {

    // Get the weight factors
    __m128 WV0 = _mm_load_ps(Wv[j][0].data_block());
    __m128 WV1 = _mm_load_ps(Wv[j][1].data_block());

    // Multiply this by the matrices that we need
    __m128 LX0 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WV0);
    __m128 LY0 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WV0);
    __m128 LZ0 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WV0);
    __m128 LR0 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WV0);

    __m128 LX1 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WV1);
    __m128 LY1 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WV1);
    __m128 LZ1 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WV1);
    __m128 LR1 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WV1);

    // Get pointers to data rows for extra speed (maybe)
    // MedialPoint *row = MP->row(j);

    for (i=0;i<=iLast;i+=step)
      {
      // Get the U, V data
      __m128 WU0 = _mm_load_ps(Wu[i][0].data_block());
      __m128 WU1 = _mm_load_ps(Wu[i][1].data_block());

      __m128 mX =  mul_4x4_4x1(LX0,LY0,LZ0,LR0,WU0);
      __m128 mXu = mul_4x4_4x1(LX0,LY0,LZ0,LR0,WU1);
      __m128 mXv = mul_4x4_4x1(LX1,LY1,LZ1,LR1,WU0);

      // Now the unit normal.
      __m128 mNRaw = crossProduct(mXu,mXv);
      __m128 mN = normalize(mNRaw);

      // Access the medial point
      MedialPoint &mp = (*MP)(i,j);

      _mm_store_ps(mp.F.data_block(),mX);
      _mm_store_ps(mp.Fu.data_block(),mXu);
      _mm_store_ps(mp.Fv.data_block(),mXv);
      _mm_store_ps(mp.F3raw.data_block(),mNRaw);
      _mm_store_ps(mp.F3.data_block(),mN);         
      }
    }
}

void MSpline::interpolateMedialSurfacePatch02(const SplineGridDefinition &grid,
                                              int iPatch,int jPatch,
                                              int step, Aux2D<MedialPoint> *MP)
{
  // Some vectors to use
  SMLVec4f **Wu = grid.patchWeight[0][iPatch];
  SMLVec4f **Wv = grid.patchWeight[1][jPatch];

  // Load the matrices in advance
  float *Cx = P[0](iPatch,jPatch).mmNbrArray;
  float *Cy = P[1](iPatch,jPatch).mmNbrArray;
  float *Cz = P[2](iPatch,jPatch).mmNbrArray;
  float *Cr = P[3](iPatch,jPatch).mmNbrArray;

  // Rows from these matrices
  __m128 CX0 = _mm_load_ps(Cx);
  __m128 CX1 = _mm_load_ps(Cx+4);
  __m128 CX2 = _mm_load_ps(Cx+8);
  __m128 CX3 = _mm_load_ps(Cx+12);

  __m128 CY0 = _mm_load_ps(Cy);
  __m128 CY1 = _mm_load_ps(Cy+4);
  __m128 CY2 = _mm_load_ps(Cy+8);
  __m128 CY3 = _mm_load_ps(Cy+12);

  __m128 CZ0 = _mm_load_ps(Cz);
  __m128 CZ1 = _mm_load_ps(Cz+4);
  __m128 CZ2 = _mm_load_ps(Cz+8);
  __m128 CZ3 = _mm_load_ps(Cz+12);

  __m128 CR0 = _mm_load_ps(Cr);
  __m128 CR1 = _mm_load_ps(Cr+4);
  __m128 CR2 = _mm_load_ps(Cr+8);
  __m128 CR3 = _mm_load_ps(Cr+12);

  // Starting and ending positions
  int jPatchStart = grid.patchStart(1,jPatch);
  int iPatchStart = grid.patchStart(0,iPatch);
  int jPatchEnd = grid.patchStart(1,jPatch+1);
  int iPatchEnd = grid.patchStart(0,iPatch+1);
  int jLast = jPatchEnd - jPatchStart;
  int iLast = iPatchEnd - iPatchStart;
  int i,j;

  for (j=0;j<=jLast;j+=step)
    {

    // Get the weight factors
    __m128 WV0 = _mm_load_ps(Wv[j][0].data_block());
    __m128 WV1 = _mm_load_ps(Wv[j][1].data_block());
    __m128 WV2 = _mm_load_ps(Wv[j][2].data_block());

    // Multiply this by the matrices that we need
    __m128 LX0 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WV0);
    __m128 LY0 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WV0);
    __m128 LZ0 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WV0);
    __m128 LR0 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WV0);

    __m128 LX1 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WV1);
    __m128 LY1 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WV1);
    __m128 LZ1 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WV1);
    __m128 LR1 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WV1);

    __m128 LX2 = mul_4x4_4x1(CX0,CX1,CX2,CX3,WV2);
    __m128 LY2 = mul_4x4_4x1(CY0,CY1,CY2,CY3,WV2);
    __m128 LZ2 = mul_4x4_4x1(CZ0,CZ1,CZ2,CZ3,WV2);
    __m128 LR2 = mul_4x4_4x1(CR0,CR1,CR2,CR3,WV2);

    // Get pointers to data rows for extra speed (maybe)
    // MedialPoint *row = MP->row(j);

    for (i=0;i<=iLast;i+=step)
      {
      // Get the U, V data
      __m128 WU0 = _mm_load_ps(Wu[i][0].data_block());
      __m128 WU1 = _mm_load_ps(Wu[i][1].data_block());
      __m128 WU2 = _mm_load_ps(Wu[i][2].data_block());

      __m128 mXuu = mul_4x4_4x1(LX0,LY0,LZ0,LR0,WU2);
      __m128 mXuv = mul_4x4_4x1(LX1,LY1,LZ1,LR1,WU1);
      __m128 mXvv = mul_4x4_4x1(LX2,LY2,LZ2,LR2,WU0);

      // Access the medial point
      MedialPoint &mp = (*MP)(i,j);
      _mm_store_ps(mp.Fuu.data_block(),mXuu);
      _mm_store_ps(mp.Fuv.data_block(),mXuv);
      _mm_store_ps(mp.Fvv.data_block(),mXvv);
      }
    }
}

/*
// Get the radius information for a spline grid
void MSpline::interpolateRadius(const SplineGridDefinition &grid,int uvFirst[2],int uvLast[2],int uvStep[2],
                                Array2D<float> &R,Array2D<float> &Ru,Array2D<float> &Rv)
{
    interpolateGrid(grid,uvFirst,uvLast,uvStep,0,0,3,3,R.row(uvFirst[1])+uvFirst[0]);
    interpolateGrid(grid,uvFirst,uvLast,uvStep,1,0,3,3,Ru.row(uvFirst[1])+uvFirst[0]);
    interpolateGrid(grid,uvFirst,uvLast,uvStep,0,1,3,3,Rv.row(uvFirst[1])+uvFirst[0]);
}

// Get the radius information for a spline grid point
void MSpline::interpolateRadiusPoint(const SplineGridDefinition &grid,int u,int v,float &r,float &ru,float &rv) {
    interpolateGridPoint(grid,u,v,0,0,3,3,&r);
    interpolateGridPoint(grid,u,v,1,0,3,3,&ru);
    interpolateGridPoint(grid,u,v,0,1,3,3,&rv);
}
*/

MSpline::MSpline(int m,int n) : DynamicBSpline2D(m,n,4) {

}

void MSpline::interpolateBoundaryCrest(MedialPoint &MP) {

  // Get a handle on the variables we will need
  SMLVec3f &X   = MP.X();
  SMLVec3f &Xu  = MP.Xu();
  SMLVec3f &Xv  = MP.Xv();
  SMLVec3f &N   = MP.N();
  float &R   = MP.R();
  float &Ru  = MP.Ru();
  float &Rv  = MP.Rv();

  // Elements of the first fundamental form
  float &E = MP.IE = Xu.squared_magnitude();
  float &G = MP.IG = Xv.squared_magnitude();
  float &F = MP.IF = dot_product(Xu,Xv);

  // Determinant of the one form
  float dI = E*G - F*F;
  float _dI = 1.0f / dI;

  // Terms going into gradR
  float CXv = (E*Rv-F*Ru)*_dI;
  float CXu = (G*Ru-F*Rv)*_dI;

  // Compute Grad R
  MP.GradR = Xv * CXv + Xu * CXu;

  // This is a crest point so sinTheta here is considered to be 0
  MP.sinTheta2 = 0.0f;

  // Compute the position and normals of boundary points
  MP.bp[0].N = MP.bp[1].N = -MP.GradR;
  MP.bp[0].X = MP.bp[1].X = X + MP.bp[0].N * R;
}

bool MSpline::interpolateBoundary(MedialPoint &MP) {

  // Get a handle on the variables we will need
  SMLVec3f &X   = MP.X();
  SMLVec3f &Xu  = MP.Xu();
  SMLVec3f &Xv  = MP.Xv();
  SMLVec3f &N   = MP.N();
  float &R   = MP.R();
  float &Ru  = MP.Ru();
  float &Rv  = MP.Rv();

  // Elements of the first fundamental form
  // float &E = MP.IE = Xu.LengthSquared();
  // float &G = MP.IG = Xv.LengthSquared();
  // float &F = MP.IF = Xu.Dot(Xv);

  float &E = MP.IE = Xu[0] * Xu[0] + Xu[1] * Xu[1] + Xu[2] * Xu[2];
  float &G = MP.IG = Xv[0] * Xv[0] + Xv[1] * Xv[1] + Xv[2] * Xv[2];
  float &F = MP.IF = Xu[0] * Xv[0] + Xu[1] * Xv[1] + Xu[2] * Xv[2];

  // Determinant of the one form
  float dI = E*G - F*F;
  float _dI = 1.0f / dI;

  // Terms going into gradR
  float CXv = (E*Rv-F*Ru)*_dI;
  float CXu = (G*Ru-F*Rv)*_dI;

  // Compute Grad R
  // MP.GradR = Xv * CXv + Xu * CXu;
  MP.GradR[0] = Xv[0] * CXv + Xu[0] * CXu;
  MP.GradR[1] = Xv[1] * CXv + Xu[1] * CXu;
  MP.GradR[2] = Xv[2] * CXv + Xu[2] * CXu;

  // Compute squared length of GradR
  // MP.sinTheta2 = 1.0 - MP.GradR.Dot(MP.GradR);
  MP.sinTheta2 = 1.0f - (MP.GradR[0] * MP.GradR[0] + MP.GradR[1] * MP.GradR[1] + MP.GradR[2] * MP.GradR[2]);

  // If the sine exists
  if (MP.sinTheta2 >= 0)
    {

    // Compute the derivatives of sine theta
    float sT = sqrtf(MP.sinTheta2);

    // Compute the normal and tangent components of the boundaries
    SMLVec3f CN = N * sT;
    SMLVec3f CT = - MP.GradR;

    // Compute the position and normals of boundary points
    MP.bp[0].N = CT - CN;
    MP.bp[1].N = CT + CN;
    MP.bp[0].X = X + MP.bp[0].N * R;
    MP.bp[1].X = X + MP.bp[1].N * R;

    // Return success
    return true;
    }

  // The point is outside of the boundary
  return false;
}
/*
bool MSpline::interpolateBoundaryO2(MedialPoint &MP,bool isCrest) {
    // Get a handle on the variables we will need
    SMLVec3f &X   = MP.X();
    SMLVec3f &Xu  = MP.Xu();
    SMLVec3f &Xv  = MP.Xv();
    SMLVec3f &Xuu = MP.Xuu();
    SMLVec3f &Xuv = MP.Xuv();
    SMLVec3f &Xvv = MP.Xvv();   
    SMLVec3f &N   = MP.NRaw();
    float &R   = MP.R();
    float &Ru  = MP.Ru();
    float &Rv  = MP.Rv();
    float &Ruu = MP.Ruu();
    float &Ruv = MP.Ruv();
    float &Rvv = MP.Rvv();

    // Elements of the first fundamental form
    float &E = MP.IE = Xu.LengthSquared();
    float &G = MP.IG = Xv.LengthSquared();
    float &F = MP.IF = Xu.Dot(Xv);

    // While we are here, lets compute the second fundamental form on the medial surface
    float IIL = Xuu.Dot(MP.N());
    float IIM = Xuv.Dot(MP.N());
    float IIN = Xvv.Dot(MP.N());

    // The Gaussian and mean curvatures
    float K = (IIL*IIN-IIM*IIM) / (E*G-F*F);
    float H = 0.5 * (E*IIN-2*F*IIM+G*IIL) / (E*G-F*F);

    // The principal curvatures
    MP.kappa1 = H + sqrtf(H*H-K);
    MP.kappa2 = H - sqrtf(H*H-K);

    // Derivatives of these elements 
    float Eu = 2 * Xu.Dot(Xuu);
    float Ev = 2 * Xu.Dot(Xuv);
    float Gu = 2 * Xv.Dot(Xuv);
    float Gv = 2 * Xv.Dot(Xvv);
    float Fu = Xuu.Dot(Xv) + Xu.Dot(Xuv);
    float Fv = Xuv.Dot(Xv) + Xu.Dot(Xvv);

    // Determinant of the one form
    float dI = E*G - F*F;
    float _dI = 1.0 / dI;

    // Derivative of the determinant and the multipliers
    float dIu = Eu*G + E*Gu - 2 * F * Fu;
    float dIv = Ev*G + E*Gv - 2 * F * Fv;
    float _dI2 = - (_dI*_dI); 
    float _dIu = _dI2 * dIu;
    float _dIv = _dI2 * dIv;

    // The direction of the medial sail in the tangent plane, so called 'B' vector 
    SMLVec3f B1 = (Xu*Ru*G - (Xu*Rv + Xv*Ru) * F + Xv*Rv * E);
    MP.GradR = B1 * _dI;
    SMLVec3f &B = MP.GradR;

    // The partial derivatives of this vector
    SMLVec3f B1u = (Xuu*Ru*G + Xu*Ruu*G + Xu*Ru*Gu) - 
                   ((Xuu*Rv + Xu*Ruv + Xuv*Ru + Xv*Ruu) * F + (Xu*Rv + Xv*Ru) * Fu) + 
                   (Xuv*Rv*E + Xv*Ruv*E + Xv*Rv*Eu);

    SMLVec3f B1v = (Xuv*Ru*G + Xu*Ruv*G + Xu*Ru*Gv) - 
                   ((Xuv*Rv + Xu*Rvv + Xvv*Ru + Xv*Ruv) * F + (Xu*Rv + Xv*Ru) * Fv) + 
                   (Xvv*Rv*E + Xv*Rvv*E + Xv*Rv*Ev);
    
    // The derivatives of the B vector
    SMLVec3f Bu = B1 * _dIu + B1u * _dI;
    SMLVec3f Bv = B1 * _dIv + B1v * _dI;

    // Sine of theta squared
    float &sT2 = MP.sinTheta2 = (isCrest) ? 0 : dI - (G * Ru*Ru - 2 * F * Ru*Rv + E * Rv*Rv);
    
    // If the sine exists
    if(!isCrest && sT2 > 0) {
        
        // Compute the derivatives of sine theta
        float sT = sqrtf(sT2);
        float sT2u = dIu - (Gu*Ru*Ru + 2*G*Ru*Ruu - 2*(Fu*Ru*Rv + F*Ruu*Rv + F*Ru*Ruv) + Eu*Rv*Rv + 2*E*Rv*Ruv);
        float sT2v = dIv - (Gv*Ru*Ru + 2*G*Ru*Ruv - 2*(Fv*Ru*Rv + F*Ruv*Rv + F*Ru*Rvv) + Ev*Rv*Rv + 2*E*Rv*Rvv);

        // Normal multiplier
        float k = _dI * sT;
        float ku = _dIu * sT + 0.5 * sT2u * _dI / sT;
        float kv = _dIv * sT + 0.5 * sT2v * _dI / sT;

        // Derivatives of the normal
        SMLVec3f Nu = Xu.Cross(Xuv) + Xuu.Cross(Xv);
        SMLVec3f Nv = Xu.Cross(Xvv) + Xuv.Cross(Xv);

        // Compute information for each boundary
        for(int d=0;d<2;d++) {
            
            // The boundary point
            BoundaryPoint &bp = MP.bp[d];
        
            // The multiplier for the normal component
            float D = d ? 1.0f : -1.0f;

            // The U vector, which is the normal to the boundary surface
            bp.N = - B + N * k * D;
            SMLVec3f &U = bp.N;

            // The derivatives
            SMLVec3f Uu = - Bu + (N * ku + Nu * k) * D;
            SMLVec3f Uv = - Bv + (N * kv + Nv * k) * D;
            bp.Nu = Uu;bp.Nv = Uv;

            // Left and right positions
            bp.X = X + U * R;
            SMLVec3f &Y = bp.X;

            // The first derivatives of the boundary vectors
            SMLVec3f Yu = Xu + U * Ru + Uu * R;
            SMLVec3f Yv = Xv + U * Rv + Uv * R;
            bp.Xu = Yu;bp.Xv = Yv;

            // The first fundamental form at the boundary
            float Eb = Yu.Dot(Yu);
            float Fb = Yu.Dot(Yv);
            float Gb = Yv.Dot(Yv);

            // The second fundamental form at the boundary
            float Lb = - Uu.Dot(Yu);
            float Mb = - Uv.Dot(Yu);
            float Nb = - Uv.Dot(Yv);

            // The Gaussian and mean curvatures
            float K = (Lb*Nb-Mb*Mb) / (Eb*Gb-Fb*Fb);
            float H = 0.5 * (Eb*Nb-2*Fb*Mb+Gb*Lb) / (Eb*Gb-Fb*Fb);

            // The principal curvatures
            bp.kappa1 = H + sqrtf(H*H-K);
            bp.kappa2 = H - sqrtf(H*H-K);
        }

        return true;
    }
    else if(isCrest || sT2 == 0) {

        // The boundary point
        BoundaryPoint &bp0 = MP.bp[0];
        BoundaryPoint &bp1 = MP.bp[1];
        
        // The U vector, which is the normal to the boundary surface
        bp0.N = bp1.N = - B;
        SMLVec3f &U = bp0.N;

        // Left and right positions
        bp1.X = bp0.X = X + U * R;
        SMLVec3f &Y = bp1.X;

        // The principal curvatures
        bp1.kappa1 = bp0.kappa1 = -1.0 / R;
        bp1.kappa2 = bp0.kappa2 = -1.0 / R;

        return true;
    }

    return false;
}
*/

void MSpline::computeCurvatures(MedialPoint &MP,bool isCrest) {

  // Get a handle on the variables we will need
  SMLVec3f &X   = MP.X();
  SMLVec3f &Xu  = MP.Xu();
  SMLVec3f &Xv  = MP.Xv();
  SMLVec3f &Xuu = MP.Xuu();
  SMLVec3f &Xuv = MP.Xuv();
  SMLVec3f &Xvv = MP.Xvv();   
  SMLVec3f &N   = MP.NRaw();
  float &R   = MP.R();
  float &Ru  = MP.Ru();
  float &Rv  = MP.Rv();
  float &Ruu = MP.Ruu();
  float &Ruv = MP.Ruv();
  float &Rvv = MP.Rvv();

  // Elements of the first fundamental form
  float &E = MP.IE;
  float &G = MP.IG;
  float &F = MP.IF;

  // While we are here, lets compute the second fundamental form on the medial surface
  float IIL = dot_product(Xuu,MP.N());
  float IIM = dot_product(Xuv,MP.N());
  float IIN = dot_product(Xvv,MP.N());

  // The Gaussian and mean curvatures
  float K = (IIL*IIN-IIM*IIM) / (E*G-F*F);
  float H = 0.5f * (E*IIN-2*F*IIM+G*IIL) / (E*G-F*F);

  // The principal curvatures
  MP.kappa1 = H + sqrtf(H*H-K);
  MP.kappa2 = H - sqrtf(H*H-K);

  // Derivatives of these elements 
  float Eu = 2 * dot_product(Xu,Xuu);
  float Ev = 2 * dot_product(Xu,Xuv);
  float Gu = 2 * dot_product(Xv,Xuv);
  float Gv = 2 * dot_product(Xv,Xvv);
  float Fu = dot_product(Xuu,Xv) + dot_product(Xu,Xuv);
  float Fv = dot_product(Xuv,Xv) + dot_product(Xu,Xvv);

  // Determinant of the one form
  float dI = E*G - F*F;
  float _dI = 1.0f / dI;

  // Derivative of the determinant and the multipliers
  float dIu = Eu*G + E*Gu - 2 * F * Fu;
  float dIv = Ev*G + E*Gv - 2 * F * Fv;
  float _dI2 = - (_dI*_dI); 
  float _dIu = _dI2 * dIu;
  float _dIv = _dI2 * dIv;

  // The direction of the medial sail in the tangent plane, so called 'B' vector 
  SMLVec3f B1 = (Xu*Ru*G - (Xu*Rv + Xv*Ru) * F + Xv*Rv * E);
  MP.GradR = B1 * _dI;
  SMLVec3f &B = MP.GradR;

  // The partial derivatives of this vector
  SMLVec3f B1u = (Xuu*Ru*G + Xu*Ruu*G + Xu*Ru*Gu) - 
                 ((Xuu*Rv + Xu*Ruv + Xuv*Ru + Xv*Ruu) * F + (Xu*Rv + Xv*Ru) * Fu) + 
                 (Xuv*Rv*E + Xv*Ruv*E + Xv*Rv*Eu);

  SMLVec3f B1v = (Xuv*Ru*G + Xu*Ruv*G + Xu*Ru*Gv) - 
                 ((Xuv*Rv + Xu*Rvv + Xvv*Ru + Xv*Ruv) * F + (Xu*Rv + Xv*Ru) * Fv) + 
                 (Xvv*Rv*E + Xv*Rvv*E + Xv*Rv*Ev);

  // The derivatives of the B vector
  SMLVec3f Bu = B1 * _dIu + B1u * _dI;
  SMLVec3f Bv = B1 * _dIv + B1v * _dI;

  // Sine of theta squared
  float &sT2 = MP.sinTheta2 = (isCrest) ? 0 : dI - (G * Ru*Ru - 2 * F * Ru*Rv + E * Rv*Rv);

  // If the sine exists
  if (!isCrest && sT2 > 0)
    {

    // Compute the derivatives of sine theta
    float sT = sqrtf(sT2);
    float sT2u = dIu - (Gu*Ru*Ru + 2*G*Ru*Ruu - 2*(Fu*Ru*Rv + F*Ruu*Rv + F*Ru*Ruv) + Eu*Rv*Rv + 2*E*Rv*Ruv);
    float sT2v = dIv - (Gv*Ru*Ru + 2*G*Ru*Ruv - 2*(Fv*Ru*Rv + F*Ruv*Rv + F*Ru*Rvv) + Ev*Rv*Rv + 2*E*Rv*Rvv);

    // Normal multiplier
    float k = _dI * sT;
    float ku = _dIu * sT + 0.5f * sT2u * _dI / sT;
    float kv = _dIv * sT + 0.5f * sT2v * _dI / sT;

    // Derivatives of the normal
    SMLVec3f Nu = vnl_cross_3d(Xu,Xuv) + vnl_cross_3d(Xuu,Xv);
    SMLVec3f Nv = vnl_cross_3d(Xu,Xvv) + vnl_cross_3d(Xuv,Xv);

    // Compute information for each boundary
    for (int d=0;d<2;d++)
      {

      // The boundary point
      BoundaryPoint &bp = MP.bp[d];

      // The multiplier for the normal component
      float D = d ? 1.0f : -1.0f;

      // The U vector, which is the normal to the boundary surface
      bp.N = - B + N * k * D;
      SMLVec3f &U = bp.N;

      // The derivatives
      SMLVec3f Uu = - Bu + (N * ku + Nu * k) * D;
      SMLVec3f Uv = - Bv + (N * kv + Nv * k) * D;
      bp.Nu = Uu;bp.Nv = Uv;

      // At this point we have Uu and Uv, and we can compute Damon's shape operator
      // The skewed frame
      float FskewValues[] = { 
        U[0],-Xu[0],-Xv[0],
        U[1],-Xu[1],-Xv[1],
        U[2],-Xu[2],-Xv[2]
      };

      SMLMatrix3f Fskew(FskewValues);
      SMLMatrix3f FskewInv = vnl_inverse(Fskew);
      SMLVec3f Fdu = FskewInv * Uu;
      SMLVec3f Fdv = FskewInv * Uv;

      // Shape operator pieces
      float a1 = Fdu[0];
      float b11 = Fdu[1];
      float b21 = Fdu[2];
      float a2 = Fdv[0];
      float b12 = Fdv[1];
      float b22 = Fdv[2];

      // Now we have the shape operator
      float trm1 = 0.5f * (b11+b22), trm2 = sqrt(b12*b21+0.25f*(b11-b22)*(b11-b22));
      float ev1 = trm1 - trm2, ev2 = trm1 + trm2;         

      // Left and right positions
      bp.X = X + U * R;
      SMLVec3f &Y = bp.X;

      // The first derivatives of the boundary vectors
      SMLVec3f Yu = Xu + U * Ru + Uu * R;
      SMLVec3f Yv = Xv + U * Rv + Uv * R;
      bp.Xu = Yu;bp.Xv = Yv;

      // The first fundamental form at the boundary
      float Eb = dot_product(Yu,Yu);
      float Fb = dot_product(Yu,Yv);
      float Gb = dot_product(Yv,Yv);

      // The second fundamental form at the boundary
      float Lb = - dot_product(Uu,Yu);
      float Mb = - dot_product(Uv,Yu);
      float Nb = - dot_product(Uv,Yv);

      // The Gaussian and mean curvatures
      float K = (Lb*Nb-Mb*Mb) / (Eb*Gb-Fb*Fb);
      float H = 0.5f * (Eb*Nb-2.0f*Fb*Mb+Gb*Lb) / (Eb*Gb-Fb*Fb);

      // The principal curvatures
      bp.kappa1 = H + sqrtf(H*H-K);
      bp.kappa2 = H - sqrtf(H*H-K);

      // double gc = bp.kappa1 * bp.kappa2;
      bp.kappa1 = 1.0f / (1.0f / ev1 - R);
      bp.kappa2 = 1.0f / (1.0f / ev2 - R);
      }
    }
  else if (isCrest || sT2 == 0)
    {

    // The boundary point
    BoundaryPoint &bp0 = MP.bp[0];
    BoundaryPoint &bp1 = MP.bp[1];

    // The U vector, which is the normal to the boundary surface
    bp0.N = bp1.N = - B;
    SMLVec3f &U = bp0.N;

    // Left and right positions
    bp1.X = bp0.X = X + U * R;
    SMLVec3f &Y = bp1.X;

    // The principal curvatures
    bp1.kappa1 = bp0.kappa1 = -1.0f / R;
    bp1.kappa2 = bp0.kappa2 = -1.0f / R;

    // Figure out the direction of tangent vector v1
    const float faFrame[] = {
      Xu[0],Xv[0],
      Xu[1],Xv[1],
      Xu[2],Xv[2]
    };
    const float faIInv[] = {
      _dI * G,  -_dI * F,
      -_dI * F, dI * E
    };
    const float faBuvT[] = {
      Bu[0],Bu[1],Bu[2],
      Bv[0],Bv[1],Bv[2]
    };
       
    vnl_matrix_fixed<float,3,2> mFrame(faFrame);
    vnl_matrix_fixed<float,2,2> mIInv(faIInv);
    vnl_matrix_fixed<float,2,3>  mBuvT(faBuvT);
    
    vnl_matrix_fixed<float,3,2> mFII = mFrame * mIInv;               
    SMLMatrix3f mDB = mFII * mBuvT;
    SMLVec3f v1 = mDB * B;
    v1 = vnl_cross_3d(v1,N).normalize();

    SMLVec3f vDBDV1 = - mDB * v1;
    SMLVec3f vDBDVn = - mDB * U;

    float faFrame3[] = { 
      U[0],-N[0],-v1[0],
      U[1],-N[1],-v1[1],
      U[2],-N[2],-v1[2]
    };
    
    SMLMatrix3f mFAFrame(faFrame3);
    SMLMatrix3f mFAFrameInv = vnl_inverse(mFAFrame);
    
    SMLVec3f mDBDv1 = mFAFrameInv * vDBDV1;
    SMLVec3f mDBDvn = mFAFrameInv * vDBDVn;

    float a1 = mDBDv1[0];
    float c1 = mDBDv1[1];
    float b1 = mDBDv1[2];
    float a2 = mDBDvn[0];
    float c2 = mDBDvn[1];
    float b2 = mDBDvn[2];

    float kappa2 = 1.0f / (1.0f / b1 - R);
    bp1.kappa2 = bp0.kappa2 = kappa2;
    }
}

/*
void MSpline::interpolate(const SplineGridDefinition &grid,int u,int v,
                          SMLVec3f &X,SMLVec3f &Xu,SMLVec3f &Xv,SMLVec3f &N,
                          float &r,SMLVec3f &XB0,SMLVec3f &NB0,SMLVec3f &XB1,SMLVec3f &NB1,float& sinTheta2)
{
    SMLVec3f NRaw;
    float ru,rv;

    interpolateMedialSurfacePoint(grid,u,v,X,Xu,Xv,NRaw,N);
    interpolateRadiusPoint(grid,u,v,r,ru,rv);
    interpolateBoundaryPoint(X,Xu,Xv,NRaw,r,ru,rv,XB0,NB0,XB1,NB1,sinTheta2);
}
*/

/*****************************************************************
 *  SplineGridDefinition
 *****************************************************************/
SplineGridDefinition::SplineGridDefinition() {  
}

SplineGridDefinition::~SplineGridDefinition() {
  for (int d=0;d<2;d++)
    {
    // Patch weights array
    for (int iPatch=0;iPatch<dim[d];iPatch++)
      {
      int res = patchStart(d,iPatch+1) - patchStart(d,iPatch);
      for (int iRes=0;iRes<=res;iRes++)
        {
        _aligned_free( patchWeight[d][iPatch][iRes] );
        }
      delete[] patchWeight[d][iPatch];
      }
    delete [] patchWeight[d];

    delete[] gridValues[d];
    delete[] knotValues[d];
    delete[] patchIndex[d];
    for (int i=0;i<size[d];i++)
      _aligned_free( W[d][i] );

    delete[] W[d];
    }
}

SplineGridDefinition *SplineGridDefinition::newSGDForCache(DynamicBSpline2D *spline,int resLevel) {
  SplineGridDefinition *sgd = new SplineGridDefinition();

  // Compute the resolution
  int res = 1 << resLevel;

  // Compute the size of the grid
  for (int d=0;d<2;d++)
    {
    // Initialize the sizes
    sgd->dim[d] = spline->dim(d)-2;
    sgd->size[d] = res*sgd->dim[d]+1;

    // Initialize the grids
    sgd->gridValues[d] = new float[sgd->size[d]];
    sgd->knotValues[d] = new int[sgd->size[d]];
    sgd->patchIndex[d] = new int[sgd->dim[d]+1];

    int iGrid = 0;
    int iPatch = 0;
    for (iPatch=0;iPatch<sgd->dim[d];iPatch++)
      {
      float u0 = spline->getKnotUV(d,iPatch+3);
      float u1 = spline->getKnotUV(d,iPatch+4);

      sgd->patchIndex[d][iPatch] = iGrid;

      for (int iRes=0;iRes<res;iRes++)
        {
        sgd->knotValues[d][iGrid] = iPatch+3;
        sgd->gridValues[d][iGrid] = (u0 * (res-iRes) + u1 * iRes) / res;
        iGrid++;
        }           
      }

    // Last location is special
    sgd->knotValues[d][iGrid] = iPatch+2;
    sgd->gridValues[d][iGrid] = spline->getKnotUV(d,iPatch+3);
    sgd->patchIndex[d][iPatch] = iGrid;

    // Initialize the weights array
    sgd->W[d] = new SMLVec4f*[sgd->size[d]];
    for (iGrid=0;iGrid<sgd->size[d];iGrid++)
      {
      sgd->W[d][iGrid] = 
        (SMLVec4f *) _aligned_malloc(sizeof(SMLVec4f) * 4, 16);
      }

    // Initialize the patch weights
    sgd->patchWeight[d] = new SMLVec4f**[sgd->dim[d]];
    for (iPatch=0;iPatch<sgd->dim[d];iPatch++)
      {
      sgd->patchWeight[d][iPatch] = new SMLVec4f*[res+1];
      for (int iRes=0;iRes<=res;iRes++)
        {
        sgd->patchWeight[d][iPatch][iRes] = 
          (SMLVec4f *) _aligned_malloc(sizeof(SMLVec4f) * 4, 16);
        }
      }
    }

  // Initialize the weights time stamp
  sgd->tsWeights = -1;
  sgd->updateWeights(spline);

  return sgd;
}

void SplineGridDefinition::updateWeights(DynamicBSpline2D *spline) {
  for (int d=0;d<2;d++)
    {
    int iWeight=0;
    for (int iGrid=0;iGrid<size[d];iGrid++)
      {
      spline->basisJet(d,knot(d,iGrid),uv(d,iGrid),W[d][iGrid]);
      }

    for (int iPatch=0;iPatch<dim[d];iPatch++)
      {
      int res = patchStart(d,iPatch+1) - patchStart(d,iPatch);
      for (int iRes=0;iRes<=res;iRes++)
        {
        spline->basisJet(d,iPatch+3,uv(d,iPatch,iRes),patchWeight[d][iPatch][iRes]);
        }
      }
    }
}

int SplineGridDefinition::getIndexForUV(int d,float u) {
  int lb = 0;
  int ub = patchStart(d,dim[d])+1;
  while (ub - lb > 1)
    {
    int test = (ub + lb) >> 1;
    if (u < uv(d,test))
      ub = test;
    else
      lb = test;
    }
  return lb;
}

PatchDataCache::PatchDataCache(MSpline *spline,SplineGridDefinition *grid,int iPatch,int jPatch,int maxLevel) {
  int i,j;

  this->spline = spline;
  this->grid = grid;
  this->iPatch = iPatch;
  this->jPatch = jPatch;
  this->maxLevel = maxLevel;

  // Allocate the time stamps

  // Initialize the time stamp
  for (int type=0;type < NUM_TYPES;type++)
    {
    tsPatch[type] = new long[maxLevel+1];       
    for (i=0;i<=maxLevel;i++)
      {
      tsPatch[type][i] = -1;
      }
    }

  // Get the width and the height
  int w = (1 << maxLevel) + 1;
  int aux = 16*w;

  // Initialize all the X arrays
  in.resize(w,w,aux);
  MP.resize(w,w,aux);
  BM[0].resize(w,w,aux);
  BM[1].resize(w,w,aux);
  MM.resize(w,w,aux);

  // Set the u and v values of the MP points  
  for (j=0;j<MP.height();j++)
    {
    for (i=0;i<MP.width();i++)
      {
      MP(i,j).u = grid->uv(0,iPatch,i);
      MP(i,j).v = grid->uv(1,jPatch,j);
      }
    }

  // Initialize the triangle array.  This part is trivial and only performed once
  nTrianglesReg = (w-1)*(w-1)*2;
  nTrianglesTrim = 0;
  PT = new PatchTriangle[nTrianglesReg+6*aux];
  nTriangles = 0;
  for (j=0;j<w-1;j++)
    {
    for (i=0;i<w-1;i++)
      {
      int v00 = MP.offset(i,j);
      PT[nTriangles].v[0]   = v00;
      PT[nTriangles].v[1]   = v00+w;
      PT[nTriangles++].v[2] = v00+1;
      PT[nTriangles].v[0]   = v00+w;
      PT[nTriangles].v[1]   = v00+w+1;
      PT[nTriangles++].v[2] = v00+1;
      }
    }

  // Initialize other structures
  idxTrimStrip = new int[6 * aux];
  idxTrimStripIdx = 0;
}

PatchDataCache::~PatchDataCache() {
  for (int type=0;type < NUM_TYPES;type++)
    {
    delete tsPatch[type];
    }
  delete idxTrimStrip;
  delete PT;
}

SplineDataCache::SplineDataCache(MSpline *spline,int maxLevel) {

  // Store the max number of levels
  this->maxLevel = maxLevel;
  this->spline = spline;

  // Get the patch sizes
  p[0] = spline->dim(0)-2;
  p[1] = spline->dim(1)-2;

  // Create a new grid definition
  grid = SplineGridDefinition::newSGDForCache(spline,maxLevel);

  // Initialize the patch arrays
  patch.resize(p[0],p[1]);
  for (int j=0;j<p[1];j++)
    for (int i=0;i<p[0];i++)
      patch(i,j) = new PatchDataCache(spline,grid,i,j,maxLevel);

  // Initialize the time stamp arrays
  for (int type=0;type < NUM_TYPES;type++)
    {
    tsPatchUB[type] = new long[maxLevel+1];
    tsPatchLB[type] = new long[maxLevel+1];

    for (int i=0;i<=maxLevel;i++)
      {
      tsPatchUB[type][i] = -1;
      tsPatchLB[type][i] = -1;
      }
    }
}

SplineDataCache::~SplineDataCache() {
  for (int j=0;j<p[1];j++)
    for (int i=0;i<p[0];i++)
      delete patch(i,j);
  delete grid;
}

void SplineDataCache::refreshMedialAllPatches(int level) {
  // Check the time stamp of the data
  if (tsPatchLB[MEDIAL][level] < spline->getControlTimeStamp())
    {

    // Update all the individual patches
    for (int i=0;i<p[0];i++)
      {
      for (int j=0;j<p[1];j++)
        {
        refreshMedialPatch(i,j,level);
        patch(i,j)->tsPatch[MEDIAL][level] = spline->getControlTimeStamp();
        }
      }

    // Update the lower bound
    tsPatchLB[MEDIAL][level] = spline->getControlTimeStamp();
    }
}   

void SplineDataCache::refreshBoundaryAllPatches(int level) {
  // Check the time stamp of the data
  if (tsPatchLB[BOUNDARY][level] < tsPatchUB[MEDIAL][level])
    {

    // Update all the individual patches
    for (int i=0;i<p[0];i++)
      {
      for (int j=0;j<p[1];j++)
        {
        refreshBoundaryPatch(i,j,level);
        patch(i,j)->tsPatch[BOUNDARY][level] = tsPatchUB[MEDIAL][level];
        }
      }

    // Update the lower bound
    tsPatchLB[BOUNDARY][level] = tsPatchUB[MEDIAL][level];
    }
}   

void PatchDataCache::refreshMedial(int level) {
  // This is simple - just give all the data to the MSpline optimized (hopefully) method
  // spline->interpolateMedialSurfacePatch(*grid,iPatch,jPatch,1 << (maxLevel-level),F,Fu,Fv,Fuu,Fuv,Fvv,Nraw,N);
  spline->interpolateMedialSurfacePatch(*grid,iPatch,jPatch,1 << (maxLevel-level),&MP);
}

/*
float PatchDataCache::integrateMedialMeasure(MedialMeasure &measure,bool trimmed,float &area) {

    // Patch dimensions
    int h = MP.height();
    int w = MP.width();
    int sz = MP.auxOffset(0) + MP.auxSize();
    int u,v,t;

    // Compute the distance over the patch
    for(v=0;v<h;v++) {
        for(u=0;u<w;u++) {
            if((!trimmed) || (in(u,v) & 1)) {
                MM(u,v) = measure.computeMedialMeasure(MP(u,v));
            }
        }
    }

    // Compute over the trim curve if necessary
    if(trimmed) {
        for(t=0;t<getCrestSize();t++) {
            MM.aux(t) = measure.computeMedialMeasure(MP.aux(t));
        }
    }

    // Distance integral
    float integral32 = 0, integral18 = 0;
    float area8 = 0,area6 = 0;

    // Integrate the distance over the patch
    for(v=0;v<h-1;v++) {
        for(u=0;u<w-1;u++) {
            
            // Make sure that the quad is located inside the trim curve if trimmed
            if(trimmed && (!(in(u,v) & in(u+1,v) & in(u,v+1) & in(u+1,v+1) & 0x01)))
                continue;
                                        
            // Some quick references to the data
            MedialPoint &M00 = MP(u,v);
            MedialPoint &M10 = MP(u+1,v);
            MedialPoint &M01 = MP(u,v+1);
            MedialPoint &M11 = MP(u+1,v+1);

            // Compute the area of the quad (I'll cheat here and pretend that the quad is flat and
            // has a normal that's the average
            SMLVec3f N4 = M00.N() + M10.N() + M01.N() + M11.N();
            float A8 = fabs((M11.X() - M00.X()).Cross(M10.X()-M01.X()).Dot(N4));

            // Compute the average distance
            float d4 = MM(u,v) + MM(u+1,v) + MM(u,v+1) + MM(u+1,v+1);
            
            // Compute the integral
            integral32 += A8 * d4;
            area8 += A8;
        }
    }

    // Integrate the auxilary triangles
    if(trimmed) {
        for(t=0;t<idxTrimStripIdx;) {
            // Get the three vertices
            int v0 = idxTrimStrip[t++];
            int v1 = idxTrimStrip[t++];
            int v2 = idxTrimStrip[t++];

            // Quick-refs to the data
            MedialPoint &M0 = MP(v0);
            MedialPoint &M1 = MP(v1);
            MedialPoint &M2 = MP(v2);

            // Compute the area of the quad (I'll cheat here and pretend that the quad is flat and
            // has a normal that's the average
            SMLVec3f N3 = M0.N() + M1.N() + M2.N();
            float A6 = fabs((M2.X()-M0.X()).Cross(M0.X()-M1.X()).Dot(N3));

            // Compute the average distance
            float d3 = MM(v0) + MM(v1) + MM(v2);
            
            // Compute the integral
            integral18 += A6 * d3;
            area6 += A6;
        }
    }

    // Compute area
    area = area8 / 8.0f + area6 / 6.0f;

    // Adjust the integrals by 1/32 and 1/18
    return integral32 / 32.0f + integral18 / 18.0f;
}

// A triangle formed by three medial points in a patch
struct MedialTriangle {
    // Indices of the triangle vertices
    int v[3];

    // Boundary and medial area elements
    SMLVec3f AM,AB[2];

    // Boundary and medial areas
    float am,ab[2];
};
*/
/*
void PatchDataCache::refreshMedialAreaWeights(int level) {

    // Patch dimensions
    int h = MP.height();
    int w = MP.width();
    int sz = MP.auxOffset(0) + getCrestSize();
    int u,v,t;

    // Resize the area arrays
    awMed.resize(sz);
    awMedTrimmed.resize(sz);
    for(t=0;t<sz;t++) {
        awMed[t] = 0;
        awMedTrimmed[t] = 0;
    }

    awMedTotal = 0;
    awMedTrimmedTotal = 0;

    // Create a list of triangles and compute their areas
    for(v=0;v<h-1;v++) {
        for(u=0;u<w-1;u++) {            
            int v00 = MP.offset(u,v);
            int v10 = v00+1;
            int v01 = v00+w;
            int v11 = v10+w;

            float a1 = triangleArea(MP(v00).F.data(),MP(v01).F.data(),MP(v10).F.data());
            float a2 = triangleArea(MP(v01).F.data(),MP(v11).F.data(),MP(v10).F.data());
            float area = a1 + a2;
            
            // Add areas to the accumulators
            awMedTotal += area;
            awMed[v00] += a1;
            awMed[v01] += area;
            awMed[v10] += area;
            awMed[v11] += a2;
            
            // If the triangles are internal, add to the trimmed area
            if(in(v00) & in(v01) & in(v10) & in(v11) & 0x01) {
                awMedTrimmedTotal += area;
                awMedTrimmed[v00] += a1;
                awMedTrimmed[v01] += area;
                awMedTrimmed[v10] += area;
                awMedTrimmed[v11] += a2;
            }
        }
    }

    // Integrate the auxilary triangles
    for(t=0;t<idxTrimStripIdx;) {       
        int v0 = idxTrimStrip[t++];
        int v1 = idxTrimStrip[t++];
        int v2 = idxTrimStrip[t++];

        float area = triangleArea(MP(v0).F.data(),MP(v1).F.data(),MP(v2).F.data());
        awMedTrimmedTotal += area;
        awMedTrimmed[v0] += area;
        awMedTrimmed[v1] += area;
        awMedTrimmed[v2] += area;
    }
}
*/

void PatchDataCache::refreshMedialAreaWeights(int level) {
  int t;

  // Patch dimensions
  int sz = MP.auxOffset(0) + getCrestSize();

  // Resize the area arrays
  awMed.resize(sz);
  awMedTrimmed.resize(sz);
  for (t=0;t<sz;t++)
    {
    awMed[t] = 0;
    awMedTrimmed[t] = 0;
    }

  awMedTotal = 0;
  awMedTrimmedTotal = 0;

  // Create a list of triangles and compute their areas
  for (t=0;t<nTriangles;t++)
    {
    PatchTriangle &T = PT[t];   
    int v0 = T.v[0],v1 = T.v[1],v2 = T.v[2];

    // Compute the area and area vector
    triangleArea(
      MP(v0).F.data_block(),
      MP(v1).F.data_block(),
      MP(v2).F.data_block(),
      T.NMed.data_block(),
      &T.areaMed );

    // Add areas to the accumulators
    awMedTotal += T.areaMed;
    awMed[v0] += T.areaMed;
    awMed[v1] += T.areaMed;
    awMed[v2] += T.areaMed;      

    if (T.in)
      {
      awMedTrimmedTotal += T.areaMed;
      awMedTrimmed[v0] += T.areaMed;
      awMedTrimmed[v1] += T.areaMed;
      awMedTrimmed[v2] += T.areaMed;
      }
    }
}


/*
void PatchDataCache::refreshMedialAreaWeights(int level) {

    // Patch dimensions
    int h = MP.height();
    int w = MP.width();
    int sz = MP.auxOffset(0) + getCrestSize();
    int u,v,t;

    // Area accumulation arrays
    vector<SMLVec3f> AA;
    AA.resize(sz);
    SMLVec3f AB;

    // Resize the area arrays
    awMed.resize(sz);
    awMedTrimmed.resize(sz);

    // Create a list of triangles and compute their areas
    for(v=0;v<h-1;v++) {
        for(u=0;u<w-1;u++) {            
            int v00 = MP.offset(u,v);
            int v10 = v00+1;
            int v01 = v00+w;
            int v11 = v10+w;

            triangleAreaVec(MP(v00).F.data(),MP(v01).F.data(),MP(v10).F.data(),AB.data());
            AA[v00] += AB;
            AA[v01] += AB;
            AA[v10] += AB;
            
            triangleAreaVec(MP(v01).F.data(),MP(v11).F.data(),MP(v10).F.data(),AB.data());
            AA[v01] += AB;
            AA[v11] += AB;
            AA[v10] += AB;          
        }
    }

    // Integrate the auxilary triangles
    for(t=0;t<idxTrimStripIdx;) {       
        int v0 = idxTrimStrip[t++];
        int v1 = idxTrimStrip[t++];
        int v2 = idxTrimStrip[t++];

        triangleAreaVec(MP(v00).F.data(),MP(v01).F.data(),MP(v10).F.data(),AB.data());
        AA[v0] += AB;
        AA[v1] += AB;
        AA[v2] += AB;
    }

    // Compute area weight of each point for untrimmed computations
    awMedTotal = 0;
    for(t=0;t<AA.size();t++) {
        awMed[t] = MP(t).N().Dot(AA[t]);
        awMedTotal += awMed[t];
    }

    // Now add up only the weights of the trimmed points
    awMedTrimmedTotal = 0;  
    for(t=0;t<idxTrimmed.size();t++) {
        int idx = idxTrimmed[t];
        awMedTrimmed[idx] = MP(idx).N().Dot(AA[idx]);
        awMedTrimmedTotal += awMedTrimmed[idx];
    }
}
*/
/*
void PatchDataCache::refreshBoundaryAreaWeights(int level) {

    // Patch dimensions
    int h = MP.height();
    int w = MP.width();
    int sz = MP.auxOffset(0) + getCrestSize();
    int u,v,t;

    // Resize the area arrays
    awBnd[0].resize(sz);
    awBnd[1].resize(sz);
    awBndTotal = 0;

    // Clear the contents of the accumulation arrays
    for(t=0;t<sz;t++) {
        awBnd[0][t] = 0;
        awBnd[1][t] = 0;
    }

    // Create a list of triangles and compute their areas
    for(v=0;v<h-1;v++) {
        for(u=0;u<w-1;u++) {            
            int v00 = MP.offset(u,v);
            int v10 = v00+1;
            int v01 = v00+w;
            int v11 = v10+w;

            // Make sure that the quad is located inside
            if((in(v00) & in(v01) & in(v10) & in(v11) & 0x01)) {
                for(int d=0;d<2;d++) {  
                    float a1 = triangleArea(MP(v00).bp[d].X.data(),MP(v01).bp[d].X.data(),MP(v10).bp[d].X.data());
                    float a2 = triangleArea(MP(v01).bp[d].X.data(),MP(v11).bp[d].X.data(),MP(v10).bp[d].X.data());
                    float area = a1 + a2;
                
                    // Add areas to the accumulators
                    awBndTotal += area;
                    awBnd[d][v00] += a1;
                    awBnd[d][v01] += area;
                    awBnd[d][v10] += area;
                    awBnd[d][v11] += a2;
                }
            }
        }
    }

    // Integrate the auxilary triangles
    for(t=0;t<idxTrimStripIdx;) {       
        int v0 = idxTrimStrip[t++];
        int v1 = idxTrimStrip[t++];
        int v2 = idxTrimStrip[t++];

        for(int d=0;d<2;d++) {
            float area = triangleArea(MP(v0).bp[d].X.data(),MP(v1).bp[d].X.data(),MP(v2).bp[d].X.data());

            // Add areas to the accumulators
            awBndTotal += area;
            awBnd[d][v0] += area;
            awBnd[d][v1] += area;
            awBnd[d][v2] += area;
        }
    }
}
*/

void PatchDataCache::refreshBoundaryAreaWeights(int level) {
  int t;

  // Patch dimensions
  int sz = MP.auxOffset(0) + getCrestSize();

  // Resize the area arrays
  awBnd[0].resize(sz);
  awBnd[1].resize(sz);
  awBndTotal = 0;

  // Clear the contents of the accumulation arrays
  for (t=0;t<sz;t++)
    {
    awBnd[0][t] = 0;
    awBnd[1][t] = 0;
    }

  // Go through the list of triangles and compute triangle normals and area
  for (t=0;t<nTriangles;t++)
    {
    PatchTriangle &T = PT[t];
    if (T.in)
      {
      int v0 = T.v[0],v1 = T.v[1],v2 = T.v[2];
      for (int d=0;d<2;d++)
        {
        float area;
        triangleArea(
          MP(v0).bp[d].X.data_block(),
          MP(v1).bp[d].X.data_block(),
          MP(v2).bp[d].X.data_block(),
          T.NBnd[d].data_block(),
          &area);

        // Add areas to the accumulators
        awBndTotal += area;
        awBnd[d][v0] += area;
        awBnd[d][v1] += area;
        awBnd[d][v2] += area;
        T.areaBnd[d] = area;
        }
      }
    }
}


/*
// This is the version for normal weighted signed areas
void PatchDataCache::refreshBoundaryAreaWeights(int level) {

    // Patch dimensions
    int h = MP.height();
    int w = MP.width();
    int sz = MP.auxOffset(0) + getCrestSize();
    int u,v,t;

    SMLVec3f AB;

    // Area accumulation arrays
    vector<SMLVec3f> AA[2];
    AA[0].resize(sz);
    AA[1].resize(sz);

    // Resize the area arrays
    awBnd[0].resize(sz);
    awBnd[1].resize(sz);

    // Create a list of triangles and compute their areas
    for(v=0;v<h-1;v++) {
        for(u=0;u<w-1;u++) {            
            int v00 = MP.offset(u,v);
            int v10 = v00+1;
            int v01 = v00+w;
            int v11 = v10+w;

            // Make sure that the quad is located inside
            if((in(v00) & in(v01) & in(v10) & in(v11) & 0x01)) {
                for(int d=0;d<2;d++) {
                    triangleAreaVec(MP(v00).bp[d].X.data(),MP(v01).bp[d].X.data(),MP(v10).bp[d].X.data(),AB.data());
                    AA[d][v00] += AB;
                    AA[d][v01] += AB;
                    AA[d][v10] += AB;
                    
                    triangleAreaVec(MP(v01).bp[d].X.data(),MP(v11).bp[d].X.data(),MP(v10).bp[d].X.data(),AB.data());
                    AA[d][v01] += AB;
                    AA[d][v11] += AB;
                    AA[d][v10] += AB;
                }
            }
        }
    }

    // Integrate the auxilary triangles
    for(t=0;t<idxTrimStripIdx;) {       
        int v0 = idxTrimStrip[t++];
        int v1 = idxTrimStrip[t++];
        int v2 = idxTrimStrip[t++];

        for(int d=0;d<2;d++) {
            triangleAreaVec(MP(v0).bp[d].X.data(),MP(v1).bp[d].X.data(),MP(v2).bp[d].X.data(),AB.data());
            AA[d][v0] += AB;
            AA[d][v1] += AB;
            AA[d][v2] += AB;
        }
    }

    // Multiply each normal by its area and also compute total area
    awBndTotal = 0;
    for(t=0;t<idxTrimmed.size();t++) {
        int idx = idxTrimmed[t];

        float A0 = - MP(idx).bp[0].N.Dot(AA[0][idx]);
        float A1 = MP(idx).bp[1].N.Dot(AA[1][idx]);
        awBnd[0][idx] = A0;
        awBnd[1][idx] = A1;
        awBndTotal += A0;
        awBndTotal += A1;

        // This indicates a kink in the boundary
        // if(A0 < 0 || A1 < 0) {
        //  cout << "Negative area at patch " << iPatch << "," << jPatch << "  Index " << idx << endl;
        // }
    }
}
*/

float PatchDataCache::integrateMedialMeasure(MedialMeasure &measure,float &area) {
  int t;
  float total = 0;
  const float scale = 1.0f / 6.0f;

  // Make sure that the area weights have been computed
  if (tsPatch[AW_MEDIAL][maxLevel] < tsPatch[BOUNDARY][maxLevel])
    {
    refreshMedialAreaWeights(maxLevel);
    tsPatch[AW_MEDIAL][maxLevel] = tsPatch[BOUNDARY][maxLevel];
    }

  // Integrate the measure
  int sz = MP.auxOffset(0);
  for (t=0;t<sz;t++)
    {
    total += awMed[t] * measure.computeMedialMeasure(MP(t));
    }

  area = awMedTotal * 0.5f;
  return total * scale;
}

float PatchDataCache::integrateTrimmedMedialMeasure(MedialMeasure &measure,float &area) {
  int t;
  float total = 0;
  const float scale = 1.0f / 6.0f;

  // Make sure that the area weights have been computed
  if (tsPatch[AW_MEDIAL][maxLevel] < tsPatch[BOUNDARY][maxLevel])
    {
    refreshMedialAreaWeights(maxLevel);
    tsPatch[AW_MEDIAL][maxLevel] = tsPatch[BOUNDARY][maxLevel];
    }

  // Integrate the measure over all internal points
  for (t=0;t<idxTrimmedInteriorSize;t++)
    {
    int idx = idxTrimmed[t];
    total += awMed[idx] * measure.computeMedialMeasure(MP(idx));
    }

  area = awMedTrimmedTotal * 0.5f;
  return total * scale;
}

float PatchDataCache::integrateBoundaryMeasure(BoundaryMeasure &measure,float &area) {
  int t;
  float total = 0;
  const float scale = 1.0f / 6.0f;

  // Make sure that the area weights have been computed
  if (tsPatch[AW_BOUNDARY][maxLevel] < tsPatch[BOUNDARY][maxLevel])
    {
    refreshBoundaryAreaWeights(maxLevel);
    tsPatch[AW_BOUNDARY][maxLevel] = tsPatch[BOUNDARY][maxLevel];
    }

  // Integrate the measure over all internal points
  for (t=0;t<idxTrimmedInteriorSize;t++)
    {
    int idx = idxTrimmed[t];
    float f0 = measure.computeBoundaryMeasure(MP(idx),0);
    float f1 = measure.computeBoundaryMeasure(MP(idx),1);
    total += awBnd[0][idx] * f0 + awBnd[1][idx] * f1;
    }

  // Integrate the crest measure over all crest points
  for (t=idxTrimmedInteriorSize;t<idxTrimmed.size();t++)
    {
    int idx = idxTrimmed[t];
    float f0 = measure.computeCrestBoundaryMeasure(MP(idx));
    total += (awBnd[0][idx] + awBnd[1][idx]) * f0;
    }

  // Compute area
  area = awBndTotal * 0.5f;

  // Adjust the integrals by 1/32 and 1/18
  return total * scale;
}

/*
float PatchDataCache::integrateBoundaryMeasure(BoundaryMeasure &measure,float &area) {

    // Patch dimensions
    int h = MP.height();
    int w = MP.width();
    int sz = MP.auxOffset(0) + MP.auxSize();
    int u,v,t;

    // Compute the distance over the patch
    for(v=0;v<h;v++) {
        for(u=0;u<w;u++) {
            if(in(u,v) & 1) {
                BM[0](u,v) = measure.computeBoundaryMeasure(MP(u,v),0);
                BM[1](u,v) = measure.computeBoundaryMeasure(MP(u,v),1);
            }
        }
    }

    // Compute over the trim curve
    for(t=0;t<getCrestSize();t++) {
        BM[0].aux(t) = BM[1].aux(t) = measure.computeCrestBoundaryMeasure(MP.aux(t));
    }

    // Distance integral
    float integral32 = 0, integral18 = 0;
    float area8 = 0,area6 = 0;

    // Integrate the distance over the patch
    for(v=0;v<h-1;v++) {
        for(u=0;u<w-1;u++) {
            
            // Make sure that the quad is located inside
            if(!(in(u,v) & in(u+1,v) & in(u,v+1) & in(u+1,v+1) & 0x01)) {
                continue;
            }
                            
            // Perform computation for both sides of the object
            for(int d=0;d<2;d++) {
                // Some quick references to the data
                BoundaryPoint &B00 = MP(u,v).bp[d];
                BoundaryPoint &B10 = MP(u+1,v).bp[d];
                BoundaryPoint &B01 = MP(u,v+1).bp[d];
                BoundaryPoint &B11 = MP(u+1,v+1).bp[d];

                // Compute the area of the quad (I'll cheat here and pretend that the quad is flat and
                // has a normal that's the average
                SMLVec3f N4 = B00.N + B10.N + B01.N + B11.N;
                float A8 = fabs((B11.X - B00.X).Cross(B10.X-B01.X).Dot(N4));

                // Compute the average distance
                float d4 = BM[d](u,v) + BM[d](u+1,v) + BM[d](u,v+1) + BM[d](u+1,v+1);
                
                // Compute the integral
                integral32 += A8 * d4;
                area8 += A8;
            }
        }
    }

    // Integrate the auxilary triangles
    for(t=0;t<idxTrimStripIdx;) {
        // Get the three vertices
        int v0 = idxTrimStrip[t++];
        int v1 = idxTrimStrip[t++];
        int v2 = idxTrimStrip[t++];

        // Compute the two triangles
        for(int d=0;d<2;d++) {
            // Quick-refs to the data
            BoundaryPoint &B0 = MP(v0).bp[d];
            BoundaryPoint &B1 = MP(v1).bp[d];
            BoundaryPoint &B2 = MP(v2).bp[d];

            // Compute the area of the quad (I'll cheat here and pretend that the quad is flat and
            // has a normal that's the average
            SMLVec3f N3 = B0.N + B1.N + B2.N;
            float A6 = fabs((B2.X-B0.X).Cross(B0.X-B1.X).Dot(N3));

            // Compute the average distance
            float d3 = BM[d](v0) + BM[d](v1) + BM[d](v2);
            
            // Compute the integral
            integral18 += A6 * d3;
            area6 += A6;
        }
    }

    // Compute area
    area = area8 / 8.0f + area6 / 6.0f;

    // Adjust the integrals by 1/32 and 1/18
    return integral32 / 32.0f + integral18 / 18.0f;
}
*/

void PatchDataCache::refreshBoundary(int level) {

  // Prepare to load the data from the spline 
  int jRes = grid->patchStart(1,jPatch+1) - grid->patchStart(1,jPatch);
  int iRes = grid->patchStart(0,iPatch+1) - grid->patchStart(0,iPatch);
  int step = 1 << (maxLevel-level);

  // Clear the array of inside points
  idxTrimmed.clear();

  // Go through all the points int the patch and interpolate the boudnary
  for (int j=0;j<=jRes;j+=step)
    {
    for (int i=0;i<=iRes;i+=step)
      {
      int iVtx = MP.offset(i,j);
      if (spline->interpolateBoundary(MP(iVtx)))
        {
        in(iVtx) = 1;
        idxTrimmed.push_back(iVtx);
        }
      else
        {
        in(iVtx) = 0;
        }
      }
    }

  // Go through the triangles and mark the trimmed ones
  for (int t=0;t<nTrianglesReg;t++)
    {
    PatchTriangle &T = PT[t];
    T.in = in(T.v[0]) & in(T.v[1]) & in(T.v[2]) & 0x01;
    }

  // Record the number of interior trimmed points
  idxTrimmedInteriorSize = idxTrimmed.size();

  // Compute the trim curve
  computeTrimCurve(level);
}

// Compute the second order medial information over the patch. 
// This is time consuming and should not be done for optimization
void PatchDataCache::refreshCurvatures(int level) {

  if (tsPatch[CURVATURES][maxLevel] < tsPatch[BOUNDARY][maxLevel])
    {

    // Set the new timestamp
    tsPatch[CURVATURES][maxLevel] = tsPatch[BOUNDARY][maxLevel];

    // Compute the second surface derivatives where needed
    spline->interpolateMedialSurfacePatch02(*grid,iPatch,jPatch,1,&MP);

    // Compute curvatures everywhere
    for (int v=0;v<MP.height();v++)
      for (int u=0;u<MP.width();u++)
        spline->computeCurvatures(MP(u,v),false);
    for (int t=0;t<iCrest;t++)
      spline->computeCurvatures(MP.aux(t),true);
    }

}




bool SplineDataCache::refreshMedialPatch(int iPatch,int jPatch,int level) {

  // Check the time stamp of the patch
  if (patch(iPatch,jPatch)->tsPatch[MEDIAL][level] < spline->getPatchTS(iPatch,jPatch))
    {

    // Do the actual rebuild of the cache
    patch(iPatch,jPatch)->refreshMedial(level);

    // Update the timestamps
    for (int l=0;l<=level;l++)
      {
      patch(iPatch,jPatch)->tsPatch[MEDIAL][l] = spline->getPatchTS(iPatch,jPatch);
      tsPatchUB[MEDIAL][l] = tsPatchUB[MEDIAL][l] > patch(iPatch,jPatch)->tsPatch[MEDIAL][l] ? tsPatchUB[MEDIAL][l] : patch(iPatch,jPatch)->tsPatch[MEDIAL][l];
      }

    return true;
    }

  // If nothing changed, return false
  return false;
}


bool SplineDataCache::refreshBoundaryPatch(int iPatch,int jPatch,int level) {
  // Make sure the medial data is fresh
  refreshMedialPatch(iPatch,jPatch,level);

  // Check the time stamp of the boundary patch
  if (patch(iPatch,jPatch)->tsPatch[BOUNDARY][level] < patch(iPatch,jPatch)->tsPatch[MEDIAL][level])
    {

    patch(iPatch,jPatch)->refreshBoundary(level);

    // Update the timestamps
    for (int l=0;l<=level;l++)
      {
      patch(iPatch,jPatch)->tsPatch[BOUNDARY][l] = patch(iPatch,jPatch)->tsPatch[MEDIAL][level];
      tsPatchUB[BOUNDARY][l] = tsPatchUB[BOUNDARY][l] > patch(iPatch,jPatch)->tsPatch[BOUNDARY][l] 
                               ? tsPatchUB[BOUNDARY][l] : patch(iPatch,jPatch)->tsPatch[BOUNDARY][l];
      }

    // Patch has been modified
    return true;
    }

  // If nothing changed, return false
  return false;
}

void PatchDataCache::computeCrestFunction(CrestSearch &cs,float coordVal,float *f,float *df) {
  // Reference to the current medial point
  MedialPoint &M = MP.aux(iCrest);

  // Get the X and R partials

  // Compute the appropriate derivative
  if (cs.dir==0)
    {
    cs.u = coordVal;

    // Compute the weights
    spline->basisJet(0,iPatch+3,cs.u,cs.Wu);
    spline->basisJet(1,jPatch+3,cs.v,cs.Wv);

    spline->interpolateMedialCrestD1(iPatch,jPatch,cs.Wu,cs.Wv,M);
    //spline->interpolateMedialPoint02(iPatch,jPatch,cs.Wu,cs.Wv,M);

    SMLVec3f &Xu = M.Xu();
    SMLVec3f &Xv = M.Xv();
    SMLVec3f &Xuv = M.Xuv();
    SMLVec3f &Xuu = M.Xuu();
    SMLVec3f &N = M.NRaw();
    float Ru = M.Ru();
    float Rv = M.Rv();
    float Ruv = M.Ruv();
    float Ruu = M.Ruu();

    SMLVec3f B = Xu * Rv - Xv * Ru;
    *f = N.squared_magnitude() - B.squared_magnitude();
    *df = 2 * (dot_product(N,vnl_cross_3d(Xuu,Xv)+vnl_cross_3d(Xu,Xuv)) - 
               dot_product(B,Xuu*Rv + Xu*Ruv - (Xuv*Ru + Xv*Ruu)));
    }
  else
    {
    cs.v = coordVal;

    // Compute the weights
    spline->basisJet(0,iPatch+3,cs.u,cs.Wu);
    spline->basisJet(1,jPatch+3,cs.v,cs.Wv);

    spline->interpolateMedialCrestD2(iPatch,jPatch,cs.Wu,cs.Wv,M);
    // spline->interpolateMedialPoint02(iPatch,jPatch,cs.Wu,cs.Wv,M);

    SMLVec3f &Xu = M.Xu();
    SMLVec3f &Xv = M.Xv();
    SMLVec3f &Xuv = M.Xuv();
    SMLVec3f &Xvv = M.Xvv();
    SMLVec3f &N = M.NRaw();
    float Ru = M.Ru();
    float Rv = M.Rv();
    float Ruv = M.Ruv();
    float Rvv = M.Rvv();

    SMLVec3f B = Xu * Rv - Xv * Ru;
    *f = N.squared_magnitude() - B.squared_magnitude();
    *df = 2 * (dot_product(N,vnl_cross_3d(Xuv,Xv)+vnl_cross_3d(Xu,Xvv)) - 
               dot_product(B,Xuv*Rv + Xu*Rvv - (Xvv*Ru + Xv*Ruv)));
    } 
}

void PatchDataCache::computeCrestMedialBoundary(CrestSearch &cs) {

  // Interpolate the missing data
  MedialPoint &mp = MP.aux(iCrest);
  mp.u = cs.u;mp.v = cs.v;

  // spline->interpolateMedialPoint02(iPatch,jPatch,cs.Wu,cs.Wv,mp);
  if (cs.dir==0)
    {
    spline->interpolateMedialCrestD1Missing(iPatch,jPatch,cs.Wu,cs.Wv,mp);
    }
  else
    {
    spline->interpolateMedialCrestD2Missing(iPatch,jPatch,cs.Wu,cs.Wv,mp);
    }

  // Compute basic boudnary information
  spline->interpolateBoundaryCrest(mp);
  in.aux(iCrest) = 0x03;

  // Add the point to the in-array
  idxTrimmed.push_back(in.auxOffset(iCrest));
}

/**
 Using a combination of Newton-Raphson and bisection, find the root of a function bracketed
 between x1 and x2. The root, returned as the function value rtsafe, will be refined until
 its accuracy is known within +-xacc. funcd is a user-supplied routine that returns both the
 function value and the rst derivative of the function.
*/
void PatchDataCache::rtsafe(CrestSearch &crestSearch,float x1, float x2,float xacc)
{
  int j;
  float df,dx,dxold,f,fh,fl;
  float temp,xh,xl,rts;

  computeCrestFunction(crestSearch,x1,&fl,&df);
  if (fl == 0.0) return;

  computeCrestFunction(crestSearch,x2,&fh,&df);
  if (fh == 0.0) return;

  if (fl < 0.0)
    {
    xl=x1;
    xh=x2;
    }
  else
    {
    xh=x1;
    xl=x2;
    }
  rts=0.5f*(x1+x2); 
  dxold=fabs(x2-x1); 
  dx=dxold; 
  computeCrestFunction(crestSearch,rts,&f,&df);
  for (j=1;j<=100;j++)
    {
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) 
        || (fabs(2.0*f) > fabs(dxold*df)))
      {
      dxold=dx;
      dx=0.5f*(xh-xl);
      rts=xl+dx;
      if (xl == rts) return;
      }
    else
      {
      dxold=dx;
      dx=f/df;
      temp=rts;
      rts -= dx;
      if (temp == rts) return;
      }
    if (fabs(dx) < xacc) return;
    computeCrestFunction(crestSearch,rts,&f,&df);

    if (f < 0.0)
      xl=rts;
    else
      xh=rts;
    }
  return; 
}


// This is just a very basic root finder
void PatchDataCache::findCrest(int dir,int iIn,int jIn,int iOut,int jOut) {
  // We will be using a crest search object
  CrestSearch cs;

  // Configure the object
  cs.dir = dir;
  cs.u = grid->uv(0,grid->patchStart(0,iPatch) + iIn);
  cs.v = grid->uv(1,grid->patchStart(1,jPatch) + jIn);

  // Perform the search
  if (dir==0)
    {
    rtsafe(cs,cs.u,grid->uv(0,grid->patchStart(0,iPatch) + iOut),0.000001f);
    }
  else
    {
    rtsafe(cs,cs.v,grid->uv(1,grid->patchStart(1,jPatch) + jOut),0.000001f);
    }

  // Store the medial atom in the outgoing location.
  computeCrestMedialBoundary(cs);

  // Next crest point
  iCrest++;
}



void PatchDataCache::startMarch(int i,int j,int inSide) {
  // Indices of start and end vertices of each enumerated edge.
  int V[8][2] = {{i,j},{i+1,j},{i+1,j+1},{i,j+1},
    {i,j},{i+1,j},{i+1,j+1},{i,j+1}};

  // Find the crossing on the starting edge
  findCrest(inSide%2,V[inSide][0],V[inSide][1],V[inSide+1][0],V[inSide+1][1]);    

  // Call the march method
  doMarch(i,j,inSide);

  // Store the curve end information
  idxTrimCurves[idxTrimCurvesIdx++] = iCrest;
}

void PatchDataCache::doMarch(int i,int j,int inSide) {

  // Indices of start and end vertices of each enumerated edge.
  int V[8][2] = {{i,j},{i+1,j},{i+1,j+1},{i,j+1},
    {i,j},{i+1,j},{i+1,j+1},{i,j+1}};

  int outSide;

  // Mark this edge as visited
  // in(V[inSide][0],V[inSide][1]) |= 0x02;
  // in(V[inSide+1][0],V[inSide+1][1]) |= 0x02;

  // Search for another crossing edge
  for (outSide=inSide+1;outSide<inSide+4;outSide++)
    {
    // Is there another crossing?
    if (in(V[outSide][0],V[outSide][1]) % 2 > in(V[outSide+1][0],V[outSide+1][1]) % 2)
      {
      findCrest(outSide%2,V[outSide][0],V[outSide][1],V[outSide+1][0],V[outSide+1][1]);   
      break;
      }
    else if (in(V[outSide][0],V[outSide][1]) % 2 < in(V[outSide+1][0],V[outSide+1][1]) % 2)
      {
      findCrest(outSide%2,V[outSide+1][0],V[outSide+1][1],V[outSide][0],V[outSide][1]);   
      break;
      }
    }

  // Mark the outgoing edges
  in(V[outSide][0],V[outSide][1]) |= 0x02;
  in(V[outSide+1][0],V[outSide+1][1]) |= 0x02;

  // Figure out the type of crossing, and save the strips
  switch (outSide - inSide)
    {
    case 1 :
      if (in(V[outSide][0],V[outSide][1]) % 2 == 1)
        {
        // This is a simple case
        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-2);
        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-1);
        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[outSide][0],V[outSide][1]);
        }
      else
        {
        // This is a fan case
        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[outSide+2][0],V[outSide+2][1]);
        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-2);
        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[inSide][0],V[inSide][1]);

        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[outSide+2][0],V[outSide+2][1]);
        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-1);
        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-2);

        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[outSide+2][0],V[outSide+2][1]);
        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[outSide+1][0],V[outSide+1][1]);
        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-1);
        }
      break;
    case 2:
      if (in(V[inSide][0],V[inSide][1]) % 2 == 1)
        {
        // This is a quad case
        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[inSide][0],V[inSide][1]);
        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-1);
        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-2);

        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[inSide][0],V[inSide][1]);
        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[outSide+1][0],V[outSide+1][1]);
        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-1);
        }
      else
        {
        // This is a quad case
        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-2);
        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[outSide][0],V[outSide][1]);
        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[inSide+1][0],V[inSide+1][1]);

        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-2);
        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-1);
        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[outSide][0],V[outSide][1]);
        }
      break;
    case 3:
      if (in(V[inSide][0],V[inSide][1]) % 2 == 1)
        {
        // This is a simple case
        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-1);
        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-2);
        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[inSide][0],V[inSide][1]);
        }
      else
        {
        // This is a fan case
        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[inSide+2][0],V[inSide+2][1]);
        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-1);
        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[outSide][0],V[outSide][1]);

        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[inSide+2][0],V[inSide+2][1]);
        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-2);
        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-1);

        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[inSide+2][0],V[inSide+2][1]);
        idxTrimStrip[idxTrimStripIdx++] = in.offset(V[inSide+1][0],V[inSide+1][1]);
        idxTrimStrip[idxTrimStripIdx++] = in.auxOffset(iCrest-2);
        }
      break;
    }

  // Continue marching, unless at an edge
  outSide %= 4;
  if (outSide == 0 && j > 0)
    {
    doMarch(i,j-step,2);
    }
  else if (outSide == 2 && j+step < last)
    {
    doMarch(i,j+step,0);
    }
  else if (outSide == 1 && i+step < last)
    {
    doMarch(i+step,j,3);
    }
  else if (outSide == 3 && i > 0)
    {
    doMarch(i-step,j,1);
    }
}

void PatchDataCache::computeTrimCurve(int level) {
  // Step and size for parsing the arrays
  step = 1 << (maxLevel-level);
  last = 1 << maxLevel;

  // Initialize iCrest
  iCrest = 0;
  idxTrimStripIdx = 0;
  idxTrimCurvesIdx = 0;

  // Look for crossings along the patch border
  for (int k=0;k<last;k++)
    {

    for (int s=0;s<4;s++)
      {

      int i0,j0,i1,j1;
      switch (s)
        {
        case 0:
          i0 = k;i1 = k+1;
          j0 = 0,j1 = 0;
          break;
        case 1:
          i0 = last;i1 = last;
          j0 = k,j1 = k+1;
          break;
        case 2:
          i0 = k;i1 = k+1;
          j0 = last,j1 = last;
          break;
        case 3:
          i0 = 0;i1 = 0;
          j0 = k,j1 = k+1;
          break;
        }

      if (((in(i0,j0) & 0x01) != (in(i1,j1) & 0x01)) &&
          ((in(i0,j0) & 0x02) | (in(i1,j1) & 0x02)) == 0)
        {

        int i2 = i0 == last ? i0-1 : i0;
        int j2 = j0 == last ? j0-1 : j0;
        startMarch(i2,j2,s);
        }
      }
    }

  // Now that the trim curve has been computed, resize the triangle array to match
  nTrianglesTrim = idxTrimStripIdx / 3;
  nTriangles = nTrianglesReg + nTrianglesTrim;
  for (int i=0,j=nTrianglesReg;i<idxTrimStripIdx;j++)
    {
    PT[j].v[0] = idxTrimStrip[i++];
    PT[j].v[1] = idxTrimStrip[i++];
    PT[j].v[2] = idxTrimStrip[i++];
    PT[j].in = 1;
    }
}

