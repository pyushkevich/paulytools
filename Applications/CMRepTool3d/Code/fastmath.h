#ifndef __fastmath_h_
#define __fastmath_h_

#include "align.h"

/**
 * Perform the cross product of two 3-vectors, result
 * stored in A
 */
inline __m128 crossProduct(__m128 A,__m128 B) 
{
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

/**
 * Compute an 'area vector' of a triangle. The magnitude of the vector
 * is equal to the area of the triangle, and the direction is normal to
 * the plane of the triangle 
 */
inline void triangleAreaVec(float *A,float *B,float *C, float *outR) 
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
  outR[0] = N[0];
  outR[1] = N[1];
  outR[2] = N[2];
}

/**
 * Compute the area triangle. Also returns the normal vector whose magnitude
 * is equal to the area of the triangle
 */
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

/**
 * Compute the area of a triangle. 
 */
inline float triangleArea(float *A,float *B,float *C) 
{
  __m128 r0,r1,r2;
  float xArea;

  // Load A, B and C
  r0 = _mm_loadu_ps(A);
  r1 = _mm_loadu_ps(C);
  r2 = _mm_loadu_ps(B);

  // Compute differences
  r0 = _mm_sub_ps(r0,r2);
  r1 = _mm_sub_ps(r1,r2);

  // Compute cross product
  r0 = crossProduct(r0,r1);

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
  _mm_store_ss(&xArea,r0);

  // Return the area value
  return xArea;
}

/**
 * Normalize a 3-vector
 */
inline __m128 normalize(__m128 A) 
{
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

/** 
 * Perform simultaneous dot product of four vectors with another vector
 */
inline __m128 mul_4x4_4x1(
  const __m128 &R0,const __m128 &R1,const __m128 &R2,
  const __m128 &R3,const __m128 &RB) 
{
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

#endif // __fastmath_h_
