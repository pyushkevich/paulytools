#include "bspline.h"
#include <vnl/algo/vnl_qr.h>
#include <xmmintrin.h>

#ifndef _WIN32
void *_mm_malloc(size_t size,size_t align) {
  void *ptr = NULL;
  posix_memalign(&ptr,size,align);
  return ptr;
}
void _mm_free(void *ptr) {
  free(ptr);
}
#endif

// Matrix multiplication routine
// Perform simultaneous dot product of four vectors with another vector
__m128 mul_4x4_4x1(const __m128 &R0,const __m128 &R1,const __m128 &R2,const __m128 &R3,const __m128 &RB) {
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

KnotVector::KnotVector(int m,int z) {
  int i;

  // Save the m
  this->m = m;
  this->nk = m + 5;

  // Assign basis function based on z
  if (z==2)
    {
    // Standard Spline Matrix 
    const float sixth = 1.0f/6.0f;
    float coeff[] = {
      -sixth, 0.5f, -0.5f, sixth, 
      0.5f, -1.0f, 0.5f, 0.0f, 
      -0.5f, 0.0f, 0.5f, 0.0f, 
      sixth, 4.0f*sixth, sixth, 0.0f};

    // Unitialize the regular spline matrix
    UB = (MySMLVec4f *)_mm_malloc(sizeof(MySMLVec4f)*4,32);
    UB[0].Set(-sixth,0.5f,-0.5f,sixth);
    UB[1].Set(0.5f,-1.0f,0.0f,4.0*sixth);
    UB[2].Set(-0.5f,0.5f,0.5f,sixth);
    UB[3].Set(sixth,0.0f,0.0f,0.0f);

    // Initialize the jet multipliers
    jetMul = (MySMLVec4f*)_mm_malloc(sizeof(MySMLVec4f)*4,32);
    jetMul[0].Set(1.0f*m*m*m,1.0f*m*m,1.0f*m,1.0f);
    jetMul[1].Set(3.0f*m*m*m,2.0f*m*m,1.0f*m,0.0f);
    jetMul[2].Set(6.0f*m*m*m,2.0f*m*m,0.0f,  0.0f);
    jetMul[3].Set(6.0f*m*m*m,0.0f,    0.0f,  0.0f);

    // Set the basis function
    basisJetFN = &KnotVector::basisJetMtx;
    }
  else
    {
    basisJetFN = &KnotVector::basisJetRcr;
    jetMul = NULL;
    UB = NULL;
    }

  // Create the knot array (by uniform interpolation)
  knots = new float[nk];
  float step = 1.0f / (nk+1-2*z);

  // Set the knot sequences
  for (i=0;i<z;i++)
    knots[i] = 0.0f;
  for (i=z;i<nk-z;i++)
    knots[i] = (i+1-z) * step;
  for (i=nk-z;i<nk;i++)
    knots[i] = 1.0f;
}

KnotVector::~KnotVector() {
  delete knots;
  if (UB)
    {
    _mm_free(UB);
    _mm_free(jetMul);
    }
}

void KnotVector::basisJetMtx(int i,float u,MySMLVec4f *N) {
  __m128 r0,r1,r2,r3,r4,r5,r6,r7;
  __m128 M0,M1,M2,M3;
  float *jm = (float *)jetMul;
  float *ub = (float *)UB;

  // Load u                                                   //  3       2       1       0
  r0 = _mm_load_ss(&u);                                       //                          u
  r1 = _mm_load_ss(knots+i);                                  //                          knot
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
  r4 = _mm_load_ps(jm);
  r5 = _mm_load_ps(jm+4);
  r6 = _mm_load_ps(jm+8);
  r7 = _mm_load_ps(jm+12);

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
  _mm_storeu_ps(N[0].data(),mul_4x4_4x1(M0,M1,M2,M3,r4));
  _mm_storeu_ps(N[1].data(),mul_4x4_4x1(M0,M1,M2,M3,r5));
  _mm_storeu_ps(N[2].data(),mul_4x4_4x1(M0,M1,M2,M3,r6));
}

void KnotVector::basisJetRcr(int i,float u,MySMLVec4f *N) {
  // We hardcode p because it's more efficient
  const int p = 3;
  const int n = 2;
  float *U = knots;

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
    N[0].data()[j]=ndu[j][p];  // basis function

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
      N[k].data()[r]=d;
      j=s1; s1=s2; s2=j;  // switch rows
      }
    }

  r=p;
  for ( k=1; k<=n; k++ )
    {
    for ( j=0; j<=p; j++ )
      N[k].data()[j]*=r;
    r*=p-k;
    }    
}

void KnotVector::getKnotIndexSequence(vector<float> &in,vector<int> &out) {
  out.clear();
  int ki = 0;
  for (int i=0;i<in.size();i++)
    {
    while (getParmAtKnot(ki+1) <= in[i] && in[i] < 1)
      ki++;
    out.push_back(ki);
    }
}

int KnotVector::getKnotAtParm(float u) {
  int n = nk-3;
  if (u == knots[n+1]) return n;
  int lo = 3,hi = n+1,mid = (lo+hi)/2;
  while (u < knots[mid] || u>=knots[mid+1])
    {
    if (u < knots[mid]) hi = mid;
    else lo = mid;
    mid = (lo+hi)/2;
    }
  return mid;
}


BSpline1D::BSpline1D(int m,int k,int z)
: kv(m,z)
{
  this->k = k;
  this->m = m;
  P = new ControlPoint1D*[k];
  for (int i=0;i<k;i++)
    P[i] = new ControlPoint1D[m+1];
}

BSpline1D::~BSpline1D()  {
  for (int i=0;i<k;i++)
    delete P[i];
  delete[] P;
}

void BSpline1D::setControl(int i,int d,float val) {
  // Set the value
  P[d][i].x = val;

  // Set the value in all neighbour matrices
  for (int l=0;l<4;l++)
    {
    int iPatch = i-l;
    if (iPatch >= 0)
      {
      // Update the neighbour matrix
      MySMLVec4f &nbr = P[d][iPatch].nbr;
      nbr.data()[l] = val;
      }
    }
}

void BSpline1D::interpolatePoint(int i,MySMLVec4f *W,int o,int d1,int d2,float *out) {
  for (int d=d1;d<=d2;d++)
    {
    // Find the control point
    ControlPoint1D &C = P[d][i];

    // This is the actual computation: a (1x4)*(4x1) multiplication
    *out = W[o].Dot(P[d][i].nbr);

    // We just increment the pointer...
    out++;
    }   
}

void BSpline1D::fitToPoints(const MatrixType &Q) 
{
  int i;

  // M is the index of the last point Q
  int m = Q.rows()-1;

  // N is the index of the last control point
  int n = this->m;

  // D is the length of the curve Q
  double l = 0;
  for (i=1;i<=m;i++)
    l += (Q.get_row(i) - Q.get_row(i-1)).two_norm();

  // Construct an array of u-values for the Q points
  vector<float> uk;
  uk.push_back(0.0f);
  for (i=1;i<m;i++)
    {
    uk.push_back(uk.back() + (Q.get_row(i) - Q.get_row(i-1)).two_norm() / l);
    }
  uk.push_back(1.0f);

  // Find knot indices of each uk
  vector<int> uki;
  kv.getKnotIndexSequence(uk,uki);

  // Constuct a 'band' matrix of N_i,p(u_j)
  MatrixType BN(m+1,n+1); BN.fill(0);
  for (int r=0;r<=m;r++)
    {
    MySMLVec4f W[4];
    kv.basisJet(uki[r],uk[r],W);
    int c0 = uki[r]-3;
    int c1 = c0+3 < n ? c0+3 : n;
    for (int c=c0;c<=c1;c++)
      {
      BN(r,c) = W[0].data()[c-c0];
      }
    }

  // Construct matrices N,R
  MatrixType R(n-1,Q.columns()); R.fill(0);
  MatrixType N(m-1,n-1); N.fill(0);

  // R is computed as a sum
  typedef vnl_vector<double> VectorType;
  for (int c=1;c<=n-1;c++)
    {
    for (int r=1;r<=m-1;r++)
      {
      VectorType L = Q.get_row(r)-Q.get_row(0)*BN(r,0)-Q.get_row(m)*BN(r,n);
      VectorType Rk = R.get_row(c-1);
      Rk += L * BN(r,c);
      R.set_column(c-1,Rk);
      }
    }

  // N is a chunk of BN
  N = BN.extract(m-1,n-1,1,1);  
  MatrixType NTN = N.transpose() * N;

  // Solve for P
  vnl_qr<double> qr(NTN);
  MatrixType R1 = qr.solve(R);

  // NTN.ipSolveLinearSystem(R);

  // Set the control point matrix
  for (int p=0;p<this->k;p++)
    {
    setControl(0,p,Q(0,p));
    for (int i=1;i<=n-1;i++)
      {
      setControl(i,p,R1(i-1,p));
      }
    setControl(n,p,Q(m,p));
    }
}
