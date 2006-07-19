#ifndef __AffineTransformDescriptor_h_
#define __AffineTransformDescriptor_h_

class AffineTransformDescriptor
{
public:
  // Vector and matrix definition
  typedef vnl_vector<double> Vec;
  typedef vnl_matrix<double> Mat;

  /** Get center of rotation for affine transforms */
  virtual SMLVec3d GetCenterOfRotation(const Vec &C) const = 0;

  /** 
   * Apply the affine transform to the coefficients. That is, solve the
   * following problem for C':
   * 
   *   Model(C') = AffineTransform[ Model(C); A, b, ctr ];
   *
   * The solution depends on the definition of the model. Hence this method is
   * virtual and implemented differently for different surfaces.
   */
  virtual Vec ApplyAffineTransform(
    const Vec &C, const Mat &A, const Vec &b, const Vec &ctr) const = 0;

  /** 
   * Use the Jacobian of the transformation C'(A,b) to map a direction in the
   * affine transform space into the corresponding direction in the
   * coefficient space
   */
  virtual Vec ApplyJacobianInParameters(
    const Vec &C, const Mat &dA, const Vec &db) const = 0;

  /**
   * Use the Jacobian of the transformation C'(C) to map a direction in
   * the source coefficient space to target coefficient space
   */
  virtual Vec ApplyJacobianInCoefficients(
    const Vec &dC, const Mat &A, const Vec &b, const Vec &ctr) const = 0;
};

class FourierAffineTransformDescriptor : public AffineTransformDescriptor
{
public:
  typedef vnl_vector<double> Vec;
  typedef vnl_matrix<double> Mat;

  /** 
   * Get the center of rotation given a coefficient vector. This is the first
   * coefficient C[0,0]
   */
  SMLVec3d GetCenterOfRotation(const Vec &C) const
    { return SMLVec3d(C.extract(3)); }

  /** 
   * Apply the transform to the coefficients. The rule is:
   *   C'[ij] = C[ij] + b         if i = j = 0
   *   C'[ij] = A * C[ij]         o/w
   */
  Vec ApplyAffineTransform(const Vec &C, const Mat &A, const Vec &b, const Vec &ctr) const
    {
    Vec CPrime = C;

    // For the first coefficient, add the b vector
    CPrime[0] += b[0]; CPrime[1] += b[1]; CPrime[2] += b[2];

    // For the rest of the coefficients, multiply by A
    for(size_t i = 4; i < C.size(); i+=4)
      CPrime.update(A * C.extract(3, i), i);

    return CPrime;
    }

  /**
   * Compute direction in coefficient space corresponding to a direction in
   * the affine transform space. By differentiating the expression above we
   * have:
   *
   *   dC'[ijk]/db[a] = a         if i=j=0 and a=k
   *                  = 0         o/w
   *
   *   dC'[ijk]/dA[ab] = C[ijb]   if a=k
   *                   = 0        o/w
   */
  Vec ApplyJacobianInParameters(const Vec &C, const Mat &dA, const Vec &db) const
    {
    Vec dCPrime(C.size(), 0.0);

    // Apply the variation in the first coefficient
    dCPrime[0] = db[0]; dCPrime[1] = db[1]; dCPrime[2] = db[2];

    // Apply the variation in the rest of the coefficients
    for(size_t i = 4; i < C.size(); i+=4)
      dCPrime.update(dA * C.extract(3, i), i);

    return dCPrime;
    }

  /** 
   * Compute direction in output coefficient space corresponding to a
   * direction in input coefficient space (v_C' = J(C'(C)) * v_C)
   */
  Vec ApplyJacobianInCoefficients(
    const Vec &dC, const Mat &A, const Vec &b, const Vec &ctr) const
    {
    Vec dCPrime(dC.size(), 0.0);

    // Apply the variation in the first coefficient
    dCPrime[0] = dC[0]; dCPrime[1] = dC[1]; dCPrime[2] = dC[2]; 

    // Apply the variation in the rest of the coefficients
    for(size_t i = 4; i < dC.size(); i+=4)
      dCPrime.update(A * dC.extract(3, i), i);

    return dCPrime;
    }
};

#endif // __AffineTransformDescriptor_h_
