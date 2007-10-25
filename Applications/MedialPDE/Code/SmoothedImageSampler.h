#ifndef __SmoothedImageSampler_h_
#define __SmoothedImageSampler_h_

#include <math.h>

// ERF is defined in GNU C, but not on windows. Since this software
// is not windows-oriented, we do it this way. Alternative is to use
// vnl_erf, which is going to be slower.
double erf(double);

/**
 * This class samples the function F = G * I where G is a Gaussian
 * and I is an image. It also samples the derivative of the function
 *
 * The sampling of the function is exact (up to floating point error)
 * as opposed to convolving the image with a Gaussian and trilerping
 * it. The function is also smooth. The idea is to use this function
 * in image-based optimization algorithms, where we want the objective
 * function to be as smooth as possible. This function is coded for
 * optimal speed, but of course, it is going to be considerably slower
 * than simply sampling images. 
 */
class SmoothedImageSampler
{
public:
  
  SmoothedImageSampler(
    double sigma,           // Standard deviation of gaussian
    double alpha,           // The number of sigmas at which to cut off the Gaussian
    float *img,             // Pointer to raw image values
    double *bb_start,       // Start of the bounding box
    double *bb_end);

  ~SmoothedImageSampler();

  double Sample(const double *X, double *grad_f);
  double Sample(const double *X);

private:

  int ntx, nty, ntz;
  int stride_z, stride_y;
  double sfx, sfy, sfz;
  double cutx, cuty, cutz;

  float *img;
  double *dx, *dy, *dz, *ddx, *ddy, *ddz;
  double sigma, alpha;
  double bb_start[3], bb_end[3];

  bool compute_erf_array(
    double *dx_erf,         // The output array of erf(p+i+1) - erf(p+i)
    int &k0, int &k1,       // The range of integration 0 <= k0 < k1 <= n
    double b,               // Lower bound of the bounding box
    int n,                  // Size of the bounding box in steps
    double cut,             // The distance at which to cut off
    double p,               // the value p
    double sfac)            // scaling factor 1 / (Sqrt[2] sigma)
      {
      // Determine the range of voxels along the line where to evaluate erf
      k0 = (int) floor(p - b - cut);
      k1 = (int) ceil(p - b + cut);
      if(k0 < 0) k0 = 0;
      if(k1 > n) k1 = n;
      if(k0 >= k1)
        return false;

      // Start at the first voxel
      double t = (b - p + k0) * sfac;
      double e_last = erf(t);
      for(int i = k0; i < k1; i++)
        {
        t += sfac;
        double e_now = erf(t);
        dx_erf[i] = e_now - e_last;
        e_last = e_now;
        }

      return true;
      }


  bool compute_erf_array_with_deriv(
    double *dx_erf,         // The output array of erf(p+i+1) - erf(p+i)
    double *dx_erf_deriv,   // The derivative of the output array wrt x
    int &k0, int &k1,       // The range of integration 0 <= k0 < k1 <= n
    double b,               // Lower bound of the bounding box
    int n,                  // Size of the bounding box in steps
    double cut,             // The distance at which to cut off
    double p,               // the value p
    double sfac)            // scaling factor 1 / (Sqrt[2] sigma)
      {
      // Determine the range of voxels along the line where to evaluate erf
      k0 = (int) floor(p - b - cut);
      k1 = (int) ceil(p - b + cut);
      if(k0 < 0) k0 = 0;
      if(k1 > n) k1 = n;
      if(k0 >= k1)
        return false;

      // Scaling factor for derivative computations
      double dscale = - sfac * 1.128379167;

      // Start at the first voxel
      double t = (b - p + k0) * sfac;
      double e_last = erf(t);
      double d_last = dscale * exp(-t * t);
      for(int i = k0; i < k1; i++)
        {
        t += sfac;
        double e_now = erf(t);
        double d_now = dscale * exp(-t * t);
        dx_erf[i] = e_now - e_last;
        dx_erf_deriv[i] = d_now - d_last;
        e_last = e_now;
        d_last = d_now;
        }

      return true;
      }
};

#endif
