#include "SmoothedImageSampler.h"

SmoothedImageSampler
::SmoothedImageSampler(
  double sigma,           // Standard deviation of gaussian
  double alpha,           // The number of sigmas at which to cut off the Gaussian
  float *img,             // Pointer to raw image values
  double *bb_start,       // Start of the bounding box
  double *bb_end)
{
  // Store the parameters
  this->sigma = sigma;
  this->alpha = alpha;
  this->img = img;
  for(unsigned int d = 0; d < 3; d++)
    {
    this->bb_start[d] = bb_start[d];
    this->bb_end[d] = bb_end[d];
    }

  // Number of interpixel ticks in the bounding box
  ntx = (int)(bb_end[0] - bb_start[0] + 0.5);
  nty = (int)(bb_end[1] - bb_start[1] + 0.5);
  ntz = (int)(bb_end[2] - bb_start[2] + 0.5);

  // Arrays to store the erf values in x, y and z
  dx = new double[ntx]; dy = new double[nty]; dz = new double[ntz];

  // Derivative arrays to store the erf values in x, y and z
  ddx = new double[ntx]; ddy = new double[nty]; ddz = new double[ntz];

  // Set up the stides
  stride_z = nty * ntx; 
  stride_y = ntx;

  // Get the sigma scale factors
  sfx = 1.0 / (sqrt(2.0) * sigma);
  sfy = 1.0 / (sqrt(2.0) * sigma);
  sfz = 1.0 / (sqrt(2.0) * sigma);

  // Compute the cutoffs alpha * sigma
  cutx = sigma * alpha;
  cuty = sigma * alpha;
  cutz = sigma * alpha;
}

SmoothedImageSampler
::~SmoothedImageSampler()
{
  delete dx; delete dy; delete dz;
  delete ddx; delete ddy; delete ddz;
}

void 
SmoothedImageSampler
::Sample(double *X, double *f, double *grad_f)
{
  // The bound variables for x, y, z
  int ix0, ix1, iy0, iy1, iz0, iz1;

  // Compute the ERF difference arrays
  compute_erf_array_with_deriv(dx, ddx, ix0, ix1, bb_start[0], ntx, cutx, X[0], sfx);
  compute_erf_array_with_deriv(dy, ddy, iy0, iy1, bb_start[1], nty, cuty, X[1], sfy);
  compute_erf_array_with_deriv(dz, ddz, iz0, iz1, bb_start[2], ntz, cutz, X[2], sfz);

  // Get a pointer to the output value
  double sum_wf = 0.0, sum_w = 0.0;
  double sum_wfx = 0.0, sum_wx = 0.0;
  double sum_wfy = 0.0, sum_wy = 0.0;
  double sum_wfz = 0.0, sum_wz = 0.0;

  // Loop over the voxels in the region identified
  for(int iz = iz0; iz < iz1; iz++)
    {
    double wz = dz[iz];
    double dwdz = ddz[iz];
    int oz = stride_z * iz;
    for(int iy = iy0; iy < iy1; iy++)
      {
      double wyz = wz * dy[iy];
      double dwdy = ddy[iy];
      int oyz = oz + stride_y * iy;
      for(int ix = ix0; ix < ix1; ix++)
        {
        int oxyz = oyz + ix;
        double w = wyz * dx[ix];
        double dwdx = ddx[ix];

        // Sample the image at this location
        double ival = img[oxyz];

        // Accumulate the function
        sum_wf += w * ival;
        sum_w += w;

        // Compute the derivatives 
        sum_wfx += dwdx * ival;
        sum_wx += dwdx;
        
        // Compute the derivatives 
        sum_wfy += dwdy * ival;
        sum_wy += dwdy;

        // Compute the derivatives 
        sum_wfz += dwdz * ival;
        sum_wz += dwdz;
        }
      }
    }

  // Set the output value
  *f = sum_wf / sum_w;

}
