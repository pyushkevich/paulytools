/**
 * Program Objective: to compute the shape operator on the implicit surface
 * defined by the given floating point image. The shape operator can be used to
 * derive principal curvatures and principal directions on the surface. 
 */

#include "ReadWriteImage.h"
#include <iostream>
#include <string>
#include <cmath>
#include <ctime>


#include <itkConstNeighborhoodIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIterator.h>

#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>

#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkSTLReader.h>
#include <vtkCellArray.h>

using namespace std;
  
// Define the matrix pixel type
typedef vnl_vector_fixed<float, 3> Vec;
typedef vnl_matrix_fixed<float, 3, 3> Mat;

// PI because WINDOWS SUCKS!
const double MY_PI = acos(-1.0);

/**
 * Fast one dimensional kernel interface. Used as a component by the three-dimensional
 * kernels.
 */
class FastAbstractKernel1D
{
public:
  // The order to which the kernel may be evaluated. Set to 3.
  const static unsigned int DEPTH;
  
  /** Constructor - sets the radius */
  FastAbstractKernel1D(unsigned int radius);
  virtual ~FastAbstractKernel1D();

  /** Compute the kernel values for some t between 0 and 1 */
  virtual void ComputeKernel(double t) = 0;
  
  /** Gaussian kernel derivatives, first index is the depth, second is 
   * the kernel offset, where the index is related */
  double **G;  
  
protected:
  // Radius and width of the kernel
  unsigned int r,w;
};

FastAbstractKernel1D::FastAbstractKernel1D(unsigned int radius)
{
  // Store the parameters
  r = radius;
  w = radius * 2;
  
  // Allocate the Gaussian array
  G = new double *[DEPTH+1];
  for(unsigned int d = 0;d <= DEPTH; d++)
    G[d] = new double[w];
}

FastAbstractKernel1D::~FastAbstractKernel1D()
{
  for(unsigned int d = 0;d < DEPTH; d++)
    delete G[d];
  delete G;
}

const unsigned int FastAbstractKernel1D::DEPTH = 3;

/** 
 * For some t between 0 and 1, this class is used to compute the
 * sequence G(t-r) .. G(t) .. G(t+r-1), where r is the 'radius' of 
 * the sequence, and G() is the 1-D Gaussian kernel and its derivatives
 * 
 * This class has been verified for correctness.
 */
class FastGaussianKernel1D : public FastAbstractKernel1D
{
public:
  
  // Constructor - destructor
  FastGaussianKernel1D(unsigned int radius, double sigma);
  ~FastGaussianKernel1D();

  /** Compute the kernel values for some t between 0 and 1 */
  virtual void ComputeKernel(double t);
  
protected:
  // Standard deviation of the kernel
  double s;
  
  // Normalizing term on the 1D Gaussian
  double Z1;

  // Factors of sigma, Z2 =  - 1 / 2 s^2
  double s2inv, s2inv3, Z2;

  // An array of exponents of n^2, n = 1..r
  double *expN2;

  // Temporary computation space 
  double *Gtemp;
};

FastGaussianKernel1D::FastGaussianKernel1D(unsigned int radius, double sigma)
  : FastAbstractKernel1D(radius)
{
  // Store the parameters
  s = sigma;
  
  // Compute constant terms
  s2inv = 1.0 / (s * s);
  s2inv3 = 3.0 * s2inv;

  Z1 = 1.0 / ( sqrt(2 * acos(-1.0)) * s );
  Z2 = - 0.5 *s2inv;

  // Compute the terms of the form exp(-n^2 / 2 s^2), where n is an integer between 1 and r
  expN2 = new double[r];
  for(unsigned int n = 1; n <= r; n++)
    expN2[n-1] = exp( -0.5 * s2inv * n * n ); 

  Gtemp = new double [w];
}

FastGaussianKernel1D::~FastGaussianKernel1D()
{
  delete expN2;
  delete Gtemp;
}

void FastGaussianKernel1D::ComputeKernel(double t)
{
  // Compute G(t), i.e., the kernel at position r. This is the first 
  // exponentiation that we require
  double *G0 = G[0];
  G0[r] = Gtemp[r] = Z1 * exp( Z2 * t * t );

  // Compute the similar term that is first order in t
  double M1 = exp( - s2inv * t );
  double M2 = 1.0 / M1;

  // Compute the terms G(t+n) as G(t) * M1^n * expN2[n]
  unsigned int i, j = 0; 
  for(i = 1; i < r; i++, j++)
    {
    Gtemp[r+i] = Gtemp[r+j] * M1;
    Gtemp[r-i] = Gtemp[r-j] * M2;
    G0[r+i] = Gtemp[r+i] * expN2[j];
    G0[r-i] = Gtemp[r-i] * expN2[j];
    }

  // Compute the zeroth term
  G0[0] = Gtemp[1] * M2 * expN2[r-1];

  // At this point, we should have the 1st dimension Gaussian kernel
  // fully computed. Now, we have to compute the derivatives of the 
  // Gaussian kernel. These derivatives can be obtained by multiplying
  // G[t] by polynomials. 
  for(unsigned int k = 0; k < w; k++) 
    {
    float u = (t + k - r) * s2inv;
    float u2 = u * u;
    G[1][k] = -u * G0[k];
    G[2][k] = (u2 - s2inv) * G0[k];
    G[3][k] = - u * (u2 - s2inv3) * G0[k];
    }

  // At this point, all the Gaussian derivatives should be computed
}


/** 
 * For some t between 0 and 1, this class is used to compute the
 * sequence G(t-r) .. G(t) .. G(t+r-1), where r is the 'radius' of 
 * the sequence, and G() is the 1-D Gaussian kernel and its derivatives
 * 
 * This class has been verified for correctness.
 */
class FastWindowedGaussianKernel1D : public FastGaussianKernel1D
{
public:
  
  // Constructor - destructor
  FastWindowedGaussianKernel1D(unsigned int radius, double sigma);
  virtual ~FastWindowedGaussianKernel1D() {};

  /** Compute the kernel values for some t between 0 and 1 */
  virtual void ComputeKernel(double t);
  
private:
  // Coefficients of the added polynomials in each derivative
  double a0, b0, c0, d0;
  double a1, b1, c1, d1;
  double a2, b2, c2, d2;
  double a3, b3, c3, d3;  
};

FastWindowedGaussianKernel1D
::FastWindowedGaussianKernel1D(unsigned int radius, double sigma)
  : FastGaussianKernel1D(radius, sigma)
{
  // Precompute some terms
  double u = (double) radius; 
  double u2 = u * u;
  double u4 = u2 * u2;
  double u6 = u2 * u4;
  double s2 = s * s;
  double s4 = s2 * s2;
  double s6 = s2 * s4;
  
  // This is a term common to all the coefficients of the polynomial. I use
  // the log here to reduce numerical error of taking sigma to the 7th power
  double Z3 = Z1 * exp( u2 * Z2 ) / (48 * s6);
  
  // A sixth order polynomial term is added to the gaussian function
  // resulting in zero derivatives at the ends of the interval -rad, rad.
  // The coefficients of this polynomial are a, b, c, and d, computed here
  a0 = - Z3;
  b0 = 3*(u2 + 2*s2) * Z3;
  c0 = -3*(u4 + 4*u2*s2 + 8*s4) * Z3;
  d0 = (u6 + 6*u4*s2 + 24*u2*s4 + 48*s6) * Z3; 

  // Now the coefficients of the derivatives
  a1 = 6.0 * a0; b1 = 4.0 * b0; c1 = 2.0 * c0; d1 = 0.0;
  a2 = 5.0 * a1; b2 = 3.0 * b1; c2 = c1;       d2 = 0.0;
  a3 = 4.0 * a2; b3 = 2.0 * b2; c3 = 0.0;      d3 = 0.0;
}

void FastWindowedGaussianKernel1D::ComputeKernel(double t)
{
  // Compute the Gaussian proper
  FastGaussianKernel1D::ComputeKernel(t);

  // Add the polynomial to the gaussian
  for(unsigned int k = 0; k < w; k++) 
    {
    float x = (t + k - r);
    float x2 = x * x;
    G[0][k] -= ( ( a0 * x2 + b0 ) * x2 + c0 ) * x2 + d0; 
    G[1][k] -= ( ( a1 * x2 + b1 ) * x2 + c1 ) * x; 
    G[2][k] -=   ( a2 * x2 + b2 ) * x2 + c2; 
    G[3][k] -=   ( a3 * x2 + b3 ) * x; 
    } 
}

void TestFastGaussianKernel1D()
{
  // Define constants
  unsigned int radius = 5;
  double sigma = 2.0, tTest = 0.25;
  
  // Create the object
  FastGaussianKernel1D K(radius, sigma);
  K.ComputeKernel(0.25);

  // Dump out the results
  cout << "FastGaussianKernel1D test results" << endl;
  for(unsigned int d=0;d<=3;d++)    
    {
    cout << "Order " << d << endl;
    for(unsigned int i=0;i<2*radius;i++)
      {
      cout << K.G[d][i] << "\t";
      }
    cout << endl << endl;
    }
}

void TestFastWindowedGaussianKernel1D()
{
  // Test the kernel above
  unsigned int radius = 5;
  double sigma = 1.2, tTest = 0.25;
  
  // Create the object
  FastWindowedGaussianKernel1D K(radius, sigma);
  K.ComputeKernel(0.25);

  // Dump out the results
  cout << "FastWindowedGaussianKernel1D test results" << endl;
  for(unsigned int d=0;d<=3;d++)    
    {
    cout << "Order " << d << endl;
    for(unsigned int i=0;i<2*radius;i++)
      {
      cout << K.G[d][i] << "\t";
      }
    cout << endl << endl;
    }
}



/** A base class for evaluating 3D implicit functions */
class ImplicitFunction3D
{
public:
  // Evaluate the function up to a given order
  virtual bool EvaluateAtPoint(double x, double y, double z, int maxOrder = 3) = 0;

  // The results of the evaluation
  double G0, G1[3], G2[3][3], G3[3][3][3];
};


/** 
 * Now, a class that uses the fast 1D kernel to compute the 
 * derivatives of the smoothed image function at any point in 
 * space. 
 *
 * Bounds checking: for the sake of speed, this class only evaluates at positions
 * that are inside the image and whose distance to the image box is at least equal
 * to the radius of the filter. This means that the input image may need to be 
 * padded by the radius.
 */
template <class TPixel>
class FastAbstractKernel3D : public ImplicitFunction3D
{
public:
  typedef itk::Image<TPixel,3> ImageType;
  
  // Constructor. Does not initialize the K pointers - that's up to the child class
  FastAbstractKernel3D(ImageType *img, unsigned int radius, bool spherical);
  virtual ~FastAbstractKernel3D();

  // Evaluate the Gaussian kernel at a point. Returns false if the input index
  // was out of range
  virtual bool EvaluateAtPoint(double x, double y, double z, int maxOrder = 3);

protected:
  // Three independent 1D Gaussians
  FastAbstractKernel1D *K[3];  

  // The radius of the kernel
  int r, w;
  
  // The image
  typename ImageType::Pointer I;

  // Size of the image
  int sz[3];

  // An array of flags describing whether kernel positions are used
  // for computation or not. This is used to form a circular kernel
  bool *inRange;

  // Number of locations that are in range
  unsigned int n;

  // A list of image offsets of in-range locations, of size n
  int *offset;

  // A list of 1D kernel offsets for the in-range locations, used for quick
  // lookup. The first index is the direction, the second is the derivative
  // order.
  unsigned int *Goff;
};

template<class TPixel>
FastAbstractKernel3D<TPixel>
::FastAbstractKernel3D(ImageType *img, unsigned int radius, bool spherical)
{
  // Store the image
  I = img;
  
  // Store the radius
  r = radius; w = r + r;

  // Compute the largest possible distance to the kernel's center
  double d2CenterMax = (r + 0.5) * (r + 0.5);
  
  // Compute the image strides
  typename ImageType::SizeType size = img->GetBufferedRegion().GetSize();
  sz[0] = (int) size[0];
  sz[1] = (int) size[1];
  sz[2] = (int) size[2];
  int stride2 = sz[0] * sz[1];
  int stride1 = sz[0];

  // Initialize the n counter
  n = 0;
  
  // Precompute the image offsets
  offset = new int[w * w * w]; 
  Goff = new unsigned int[3 * w * w * w];
  int *p = offset;
  for(int z = 1 - r; z <= r; z++)
    {
    for(int y = 1 - r; y <= r; y++)
      {
      for(int x = 1 - r; x <= r; x++)
        {
        if(spherical)
          {
          // Check whether this position is farther away from kernel center
          // than the radius
          double d2Center = 
            (0.5 - z) * (0.5 - z) + (0.5 - y) * (0.5 - y) + (0.5 - x) * (0.5 - x);

          // Only take the offset of points that are in the range
          if(d2Center > d2CenterMax)
            {
            continue;
            }
          }

        // Compute the offset
        offset[n] = stride2 * z + stride1 * y + x;

        // Compute the lookup positions for each of the 1D Gaussian kernels
        Goff[3 * n]     = r - x; 
        Goff[3 * n + 1] = r - y; 
        Goff[3 * n + 2] = r - z;

        n++;          
        }
      }
    }  
}

template<class TPixel>
bool FastAbstractKernel3D<TPixel>
::EvaluateAtPoint(double x, double y, double z, int maxOrder)
{
  // Compute the image offset of the coordinate
  int ix = (int) x;
  int iy = (int) y; 
  int iz = (int) z;

  // Check that the coordinates are in range
  if(ix < r - 1 || iy < r - 1 || iz < r - 1)
    return false;
  if(ix >= sz[0] - r || iy >= sz[1] - r || iz >= sz[2] - r)
    return false;

  // Get a pointer to the image data at the pixel location
  TPixel *p = I->GetBufferPointer();
  p += sz[0] * (sz[1] * iz + iy) + ix; 

  // Get the in-pixel offsets
  double tx = x - ix;
  double ty = y - iy;
  double tz = z - iz;
  
  // Compute the Gaussian kernels
  K[0]->ComputeKernel(tx);
  K[1]->ComputeKernel(ty);
  K[2]->ComputeKernel(tz);

  // Clear all the quantities of interest
  G0 = 0.0;
  
  G1[0] = 0.0; G1[1] = 0.0; G1[2] = 0.0;
  
  G2[0][0] = 0.0; G2[0][1] = 0.0; G2[0][2] = 0.0;
  G2[1][1] = 0.0; G2[1][2] = 0.0; G2[2][2] = 0.0;
  
  G3[0][0][0] = 0.0; G3[0][0][1] = 0.0; G3[0][0][2] = 0.0;
  G3[0][1][1] = 0.0; G3[0][1][2] = 0.0; G3[0][2][2] = 0.0;
  G3[1][1][1] = 0.0; G3[1][1][2] = 0.0; G3[1][2][2] = 0.0;
  G3[2][2][2] = 0.0;
  
  // Go through the image, computing all the necessary quantities
  int *pOffset = offset;
  unsigned int *pGoff = Goff;
  
  for(unsigned int i = 0; i < n; i++)
    {
    // Get the image intensity
    double q = (double) *(p + *pOffset++);

    // Get the x, y and z offsets
    unsigned int ox = *pGoff++;
    unsigned int oy = *pGoff++;
    unsigned int oz = *pGoff++;

    // Get the individual derivatives
    double gx0 = K[0]->G[0][ox], gy0 = K[1]->G[0][oy], gz0 = K[2]->G[0][oz];

    // Compute some common products to avoid repetition (34 mults instead of 60)
    double qx0 = q * gx0;
    double qx0y0 = qx0 * gy0;

    // Compute the Gaussian 
    G0          += qx0y0 * gz0;                 // q * gx0 * gy0 * gz0;

    // Only continue for higher orders
    if(maxOrder > 0) 
      {
      double gx1 = K[0]->G[1][ox], gy1 = K[1]->G[1][oy], gz1 = K[2]->G[1][oz];
      double gx2 = K[0]->G[2][ox], gy2 = K[1]->G[2][oy], gz2 = K[2]->G[2][oz];
      double gx3 = K[0]->G[3][ox], gy3 = K[1]->G[3][oy], gz3 = K[2]->G[3][oz];

      // Compute some common products to avoid repetition (34 mults instead of 60)
      double qy0 = q * gy0;
      double qz0 = q * gz0;    
      
      double qx0z0 = qx0 * gz0;
      double qy0z0 = gy0 * qz0;

      double qx0y1 = qx0 * gy1;
      double qz0x1 = qz0 * gx1;
      double qy0z1 = qy0 * gz1;
      
      G1[0]       += gx1 * qy0z0;                 // q * gx1 * gy0 * gz0;
      G1[1]       += gy1 * qx0z0;                 // q * gx0 * gy1 * gz0;
      G1[2]       += gz1 * qx0y0;                 // q * gx0 * gy0 * gz1;

      G2[0][0]    += qy0z0 * gx2;                 // q * gx2 * gy0 * gz0;
      G2[0][1]    += qz0x1 * gy1;                 // q * gx1 * gy1 * gz0;
      G2[0][2]    += qy0z1 * gx1;                 // q * gx1 * gy0 * gz1;
      G2[1][1]    += qx0z0 * gy2;                 // q * gx0 * gy2 * gz0;
      G2[1][2]    += qx0y1 * gz1;                 // q * gx0 * gy1 * gz1;
      G2[2][2]    += qx0y0 * gz2;                 // q * gx0 * gy0 * gz2;

      G3[0][0][0] += qy0z0 * gx3;                 // q * gx3 * gy0 * gz0;
      G3[0][0][1] += qz0 * gx2 * gy1;             // q * gx2 * gy1 * gz0;
      G3[0][0][2] += qy0z1 * gx2;                 // q * gx2 * gy0 * gz1;
      G3[0][1][1] += qz0x1 * gy2;                 // q * gx1 * gy2 * gz0;
      G3[0][1][2] += q * gx1 * gy1 * gz1;         // q * gx1 * gy1 * gz1;
      G3[0][2][2] += qy0 * gx1 * gz2;             // q * gx1 * gy0 * gz2;
      G3[1][1][1] += qx0z0 * gy3;                 // q * gx0 * gy3 * gz0;
      G3[1][1][2] += qx0 * gy2 * gz1;             // q * gx0 * gy2 * gz1;
      G3[1][2][2] += qx0y1 * gz2;                 // q * gx0 * gy1 * gz2;
      G3[2][2][2] += qx0y0 * gz3;                 // q * gx0 * gy0 * gz3;
      }
    // We are done!
    }
  
  return true;
}


template<class TPixel>
FastAbstractKernel3D<TPixel>
::~FastAbstractKernel3D()
{
}


/** 
 * Now, a class that uses the fast Gaussian kernel to compute the 
 * derivatives of the smoothed image function at any point in 
 * space. 
 *
 * Bounds checking: for the sake of speed, this class only evaluates at positions
 * that are inside the image and whose distance to the image box is at least equal
 * to the radius of the filter. This means that the input image may need to be 
 * padded by the radius.
 */
template <class TPixel>
class FastGaussianKernel3D : public FastAbstractKernel3D<TPixel>
{
public:
  typedef FastAbstractKernel3D<TPixel> Superclass;
  typedef itk::Image<TPixel,3> ImageType;
  
  FastGaussianKernel3D(ImageType *img, unsigned int radius, double sigma[3]);
  ~FastGaussianKernel3D();
};

template<class TPixel>
FastGaussianKernel3D<TPixel>
::FastGaussianKernel3D(ImageType *img, unsigned int radius, double sigma[3])
: Superclass(img, radius, true)
{
  // Initialize the 1D kernels
  K[0] = new FastGaussianKernel1D(radius, sigma[0]);
  K[1] = new FastGaussianKernel1D(radius, sigma[1]);
  K[2] = new FastGaussianKernel1D(radius, sigma[2]);
}

template<class TPixel>
FastGaussianKernel3D<TPixel>
::~FastGaussianKernel3D()
{
  // UnInitialize the 1D kernels
  delete K[0];
  delete K[1];
  delete K[2];
}


/** A kernel that uses the windowed gaussian */
template <class TPixel>
class FastWindowedGaussianKernel3D : public FastAbstractKernel3D<TPixel>
{
public:
  typedef FastAbstractKernel3D<TPixel> Superclass;
  typedef itk::Image<TPixel,3> ImageType;
  
  FastWindowedGaussianKernel3D(ImageType *img, unsigned int radius, double sigma[3]);
  ~FastWindowedGaussianKernel3D();
};

template<class TPixel>
FastWindowedGaussianKernel3D<TPixel>
::FastWindowedGaussianKernel3D(ImageType *img, unsigned int radius, double sigma[3])
: Superclass(img, radius, false)
{
  // Initialize the 1D kernels
  K[0] = new FastWindowedGaussianKernel1D(radius, sigma[0]);
  K[1] = new FastWindowedGaussianKernel1D(radius, sigma[1]);
  K[2] = new FastWindowedGaussianKernel1D(radius, sigma[2]);
}

template<class TPixel>
FastWindowedGaussianKernel3D<TPixel>
::~FastWindowedGaussianKernel3D()
{
  // UnInitialize the 1D kernels
  delete K[0];
  delete K[1];
  delete K[2];
}

void TestFastGaussianKernel3D()
{
  // Create a dummy image
  typedef itk::Image<float, 3> ImageType;
  ImageType::Pointer img = ImageType::New();
  
  ImageType::SizeType sz = {{40,40,40}};
  ImageType::RegionType region(sz);
  img->SetRegions(region);
  img->Allocate();
  img->FillBuffer(0.0f);

  // Set one pixel in the image
  ImageType::IndexType idx = {{3,2,3}};
  img->SetPixel(idx, 3.0f);

  // Create a Gaussian kernel of radius 4
  double sigma[] = {2.0, 2.0, 2.0};
  FastGaussianKernel3D<float> K(img, 6, sigma);

  // Evaluate the kernel at 4.25
  K.EvaluateAtPoint(4.25, 4.15, 4.6);

  // Print the result
  cout << "G: " << K.G0 << endl;

  cout << "G1[0] : " << K.G1[0] << endl;
  cout << "G1[1] : " << K.G1[1] << endl;
  cout << "G1[2] : " << K.G1[2] << endl;

  cout << "G2[0][0] : " << K.G2[0][0] << endl;
  cout << "G2[0][1] : " << K.G2[0][1] << endl;
  cout << "G2[0][2] : " << K.G2[0][2] << endl;
  cout << "G2[1][1] : " << K.G2[1][1] << endl;
  cout << "G2[1][2] : " << K.G2[1][2] << endl;
  cout << "G2[2][2] : " << K.G2[2][2] << endl;

  cout << "G3[0][0][0] : " << K.G3[0][0][0] << endl;
  cout << "G3[0][0][1] : " << K.G3[0][0][1] << endl;
  cout << "G3[0][0][2] : " << K.G3[0][0][2] << endl;
  cout << "G3[0][1][1] : " << K.G3[0][1][1] << endl;
  cout << "G3[0][1][2] : " << K.G3[0][1][2] << endl;
  cout << "G3[0][2][2] : " << K.G3[0][2][2] << endl;
  cout << "G3[1][1][1] : " << K.G3[1][1][1] << endl;
  cout << "G3[1][1][2] : " << K.G3[1][1][2] << endl;
  cout << "G3[1][2][2] : " << K.G3[1][2][2] << endl;
  cout << "G3[2][2][2] : " << K.G3[2][2][2] << endl;

  // Now test how many evaluations can be performed in a short time
  vector<double> v;
  for(unsigned int i = 0; i < 30000; i++)
    v.push_back(6.0 + rand() * 28.0 / RAND_MAX);
  
  bool flag = true;
  double c0 = clock();
  for(unsigned int j = 0; j < 30000; j+=3)
    {
    flag &= K.EvaluateAtPoint(v[j], v[j+1], v[j+2]);
    }
  double c1 = clock();
  cout << "10000 evaluations in " << (c1 - c0) / CLOCKS_PER_SEC << " seconds" << endl;
  cout << "Flag : " << flag << endl;
  cout << (int)(10000 * CLOCKS_PER_SEC / (c1 - c0)) << " evals per second" << endl;
}

inline void CrossProduct(double *x1,double *x2,double *x)
{
  x[0] = x1[1] * x2[2] - x1[2] * x2[1];
  x[1] = x1[2] * x2[0] - x1[0] * x2[2];
  x[2] = x1[0] * x2[1] - x1[1] * x2[0];
}

template<class Real = double>
class SpacePoint 
{
public:
  // Raw third derivative data, without duplicates
  Real Fxxx, Fxxy, Fxyy, Fyyy, Fxxz, Fxyz, Fyyz, Fxzz, Fyzz, Fzzz;

  // Raw second derivative data, without duplicates
  Real Fxx, Fxy, Fxz, Fyy, Fyz, Fzz;

  // Raw first derivative data
  Real Fx, Fy, Fz;

  // The value of the function at the point
  Real F;

  // The gradient vector, set to &Fx
  Real *G;

  // The gradient magnitude and its powers
  Real GNorm, GNorm2, GNorm3, GNorm_1, GNorm_3;

  // The normal vector to the surface (normalized G), N points to &Nx;
  Real Nx, Ny, Nz, *N; 

  // Curvatures
  Real kMin, kMax, kGauss, kMean;

  // Principal directions
  Real tMin[3], tMax[3];

  // Value of the ridge function, zero on ridges an valleys
  Real fRidgeMin, fRidgeMax, fRidgeProduct, fRidgeDominant, fRidge;

  // Whether the curvatures and the ridge have been computed
  bool flagCurvatures, flagRidge;

  /** Constructor */
  SpacePoint();
  SpacePoint(Real F0, Real F1[3], Real F2[3][3], Real F3[3][3][3]);

  /** Set the Jet values - identical to the second constructor */
  void SetFunctionJet(Real F0, Real F1[3], Real F2[3][3], Real F3[3][3][3]);

  /** Compute the curvature information at the point */
  void ComputePrincipalCurvatures();

  /** Compute the ridge function at the point */
  void ComputeRidgeFunction();

private:
  // Hidden quanitities, used to speed up computation
  Real Fx2, Fy2, Fz2;

  // Compute the change in principal curvature with respect to a given principal 
  // direction
  Real ComputeDkDtTerm1(Real *t);
};

template<class Real>
SpacePoint<Real>
::SpacePoint()
{
  double F1[3] = {0.0, 0.0, 0.0};
  double F2[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  double F3[3][3][3] = {
    {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
    {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
    {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}} };

  SetFunctionJet(0.0, F1, F2, F3);
}

template<class Real>
SpacePoint<Real>
::SpacePoint(Real F0, Real F1[3], Real F2[3][3], Real F3[3][3][3])
{
  // Use the data
  SetFunctionJet(F0, F1, F2, F3);
}

template<class Real>
void
SpacePoint<Real>
::SetFunctionJet(Real F0, Real F1[3], Real F2[3][3], Real F3[3][3][3])
{
  // Pull out the derivative data
  Fxxx = F3[0][0][0];
  Fxxy = F3[0][0][1];
  Fxyy = F3[0][1][1];
  Fyyy = F3[1][1][1];
  Fxxz = F3[0][0][2];
  Fxyz = F3[0][1][2];
  Fyyz = F3[1][1][2];
  Fxzz = F3[0][2][2];
  Fyzz = F3[1][2][2];
  Fzzz = F3[2][2][2];

  Fxx  = F2[0][0]; 
  Fxy  = F2[0][1]; 
  Fxz  = F2[0][2]; 
  Fyy  = F2[1][1]; 
  Fyz  = F2[1][2]; 
  Fzz  = F2[2][2]; 
  
  Fx   = F1[0];
  Fy   = F1[1];
  Fz   = F1[2];

  F    = F0;

  // Assign pointers
  G = &Fx;
  N = &Nx;

  // Compute squared Fx, Fy, Fz
  Fx2 = Fx * Fx;
  Fy2 = Fy * Fy;
  Fz2 = Fz * Fz;

  // Compute the grad magnitude and related quantities
  GNorm2 = Fx2 + Fy2 + Fz2;
  GNorm = sqrt(GNorm2);
  GNorm3 = GNorm2 * GNorm;
  GNorm_1 = 1.0 / GNorm;
  GNorm_3 = 1.0 / GNorm3;

  // Compute the normal vector
  Nx = - Fx * GNorm_1;
  Ny = - Fy * GNorm_1;
  Nz = - Fz * GNorm_1;

  // Clear the flags
  flagCurvatures = flagRidge = false;
}

template<class Real>
void 
SpacePoint<Real>
::ComputePrincipalCurvatures()
{
  // Compute the Shape operator, given by 
  // (I - (G G^T) / |G|^2) Hess[F] / |G| where G is grad F
  // this formula taken from Eberly's paper
  // http://www.magic-software.com/Documentation/PrincipalCurvature.pdf
  // For simplicity, I got all the computations out of Mathematica
  
  // Here are some terms that get repeated in the computation, so I save them
  Real Fxz2 = Fxz * Fxz, Fxy2 = Fxy * Fxy, Fyz2 = Fyz * Fyz;

  // Compute the shape operator (up to the factor GNorm3)
  Real S[3][3] = 
    {
      {
      Fx*(Fxy*Fy + Fxz*Fz) - Fxx*(Fy2 + Fz2), 
      Fx*(Fy*Fyy + Fyz*Fz) - Fxy*(Fy2 + Fz2), 
      -(Fxz*(Fy2 + Fz2)) + Fx*(Fy*Fyz + Fz*Fzz)
      }, 
      {
      -(Fx2*Fxy) + Fx*Fxx*Fy + Fz*(Fxz*Fy - Fxy*Fz), 
      Fx*Fxy*Fy - Fx2*Fyy + Fz*(Fy*Fyz - Fyy*Fz), 
      Fx*Fxz*Fy - Fx2*Fyz + Fz*(-(Fyz*Fz) + Fy*Fzz)
      }, 
      {
      -(Fxz*(Fx2 + Fy2)) + (Fx*Fxx + Fxy*Fy)* Fz,
      -((Fx2 + Fy2)*Fyz) + (Fx*Fxy + Fy*Fyy)*Fz, 
      (Fx*Fxz + Fy*Fyz)*Fz - (Fx2 + Fy2)*Fzz
      }
    };

  // The Gaussian curvature
  kGauss = (-(Fxz2*Fy2) + 2*Fxz*(Fx*Fy*Fyz + Fxy*Fy*Fz - Fx*Fyy*Fz)
    + Fz*(-2*Fxx*Fy*Fyz - Fxy2*Fz + Fxx*Fyy*Fz) + Fxx*Fy2*Fzz
    + 2*Fx*Fxy*(Fyz*Fz - Fy*Fzz) + Fx2*(-Fyz2 + Fyy*Fzz)) /
    (GNorm2 * GNorm2);

  // The mean curvature
  kMean = (-2*Fy*Fyz*Fz + Fyy*Fz2 - 
    2*Fx*(Fxy*Fy + Fxz*Fz) + 
    Fxx*(Fy2 + Fz2) + Fy2*Fzz + 
    Fx2*(Fyy + Fzz))/
    ( 2 * GNorm3 );

  // Compute the principal curvatures
  Real M2K = sqrt(kMean * kMean-kGauss);
  kMax = kMean + M2K; kMin = kMean - M2K;

  // Compute the principal vectors of the modified shape operator
  Real kAdd = kMax * GNorm3;
  S[0][0] += kAdd;
  S[1][1] += kAdd;
  S[2][2] += kAdd;

  Real v[3][3];
  CrossProduct(S[0],S[1],v[0]);
  CrossProduct(S[1],S[2],v[1]);
  CrossProduct(S[2],S[0],v[2]);

  // The principal direction has the largest norm
  Real nv[3] = 
    { v[0][0]*v[0][0] + v[0][1]*v[0][1] + v[0][2]*v[0][2],
      v[1][0]*v[1][0] + v[1][1]*v[1][1] + v[1][2]*v[1][2],
      v[2][0]*v[2][0] + v[2][1]*v[2][1] + v[2][2]*v[2][2] };

  int iMax = (nv[0] > nv[1]) ? (nv[0] > nv[2] ? 0 : 2) : (nv[1] > nv[2] ? 1 : 2);
  Real factor1 = 1.0 / sqrt(nv[iMax]);
  tMax[0] = v[iMax][0] * factor1;
  tMax[1] = v[iMax][1] * factor1;
  tMax[2] = v[iMax][2] * factor1;

  // The second principal direction is the cross-product of the first with the normal
  CrossProduct(tMax, G, tMin);
  Real factor2 = 1.0 / GNorm;
  tMin[0] *= factor2;
  tMin[1] *= factor2;
  tMin[2] *= factor2;

  // Set the flag
  flagCurvatures = true;
}

template<class Real>
Real
SpacePoint<Real>
::ComputeDkDtTerm1(Real *t)
{
  Real dkds = 0.0;
  Real tx = t[0], ty = t[1], tz = t[2];
  Real tx2 = tx * tx, ty2 = ty * ty, tz2 = tz * tz;

  dkds +=       Fxxx * tx2 * tx;
  dkds += 3.0 * Fxxy * tx2 * ty;
  dkds += 3.0 * Fxyy * tx * ty2;
  dkds +=       Fyyy * ty2 * ty;
  dkds += 3.0 * Fxxz * tx2 * tz;
  dkds += 6.0 * Fxyz * tx * ty * tz;
  dkds += 3.0 * Fyyz * ty2 * tz;
  dkds += 3.0 * Fxzz * tx * tz2;
  dkds += 3.0 * Fyzz * ty * tz2;
  dkds +=       Fzzz * tz * tz2;

  return dkds;
}


/**
 * Compute the ridge function, equation 10 from the paper Belyaev et al.,
 * "Ridges and ravines on implicit surfaces "
 */
template<class Real>
void
SpacePoint<Real>
::ComputeRidgeFunction()
{
  // Compute curvatures if needed 
  if(!flagCurvatures)
    ComputePrincipalCurvatures();
  
  // Compute the third derivative terms
  Real T1max = ComputeDkDtTerm1(tMax);
  Real T1min = ComputeDkDtTerm1(tMin);

  // Compute the common parts in the second term
  Real T2x = Fxx * Nx + Fxy * Ny + Fxz * Nz;
  Real T2y = Fxy * Nx + Fyy * Ny + Fyz * Nz;
  Real T2z = Fxz * Nx + Fyz * Ny + Fzz * Nz;

  // Compute the second terms
  Real T2max = T2x * tMax[0] + T2y * tMax[1] + T2z * tMax[2];
  Real T2min = T2x * tMin[0] + T2y * tMin[1] + T2z * tMin[2];

  // Compute the ridge implicit functions
  fRidgeMax = (T1max + 3 * kMax * T2max) * GNorm_1;
  fRidgeMin = (T1min + 3 * kMin * T2min) * GNorm_1;
  
  // Compute the final answer
  fRidgeProduct = fRidgeMin * fRidgeMax;
  fRidgeDominant = fabs(kMin) < fabs(kMax)
    ? fRidgeMax : fRidgeMin;

  // Choose which ridge to use in tracking
  fRidge = fRidgeDominant;

  // Set the ridge flag
  flagRidge = true;
}

void TestCurvatureComputation()
{
  // Get some data from Mathematica. The function is 
  // F = Sin[x*y]Cos[z] - Cos[y*z]Sin[x] + Sin[z*x]Cos[y]
  // and the point is {x -> 0.45, y -> 0.34, z -> 0.23};
  double F0 = -0.187845;
  
  double F1[3] = 
    {-0.354843, 0.406393, 0.398778};
  
  double F2[3][3] = 
    {
      {0.41133, 0.879476, 0.874938}, 
      {0.879476, -0.104511, -0.182768}, 
      {0.874938, -0.182768, -0.117986}
    };

  double F3[3][3][3] = 
    {
      {
        {0.848464, -0.156957, -0.0746632}, 
        {-0.156957, -0.36799, -0.407592}, 
        {-0.0746632, -0.407592, -0.354741}
      }, 
      {
        {-0.156957, -0.36799, -0.407592}, 
        {-0.36799, -0.053648, -0.216072}, 
        {-0.407592, -0.216072, -0.132085}
      }, 
      {
        {-0.0746632, -0.407592, -0.354741}, 
        {-0.407592, -0.216072, -0.132085}, 
        {-0.354741, -0.132085, -0.0520397}
      }
    };

  // Create a space point
  SpacePoint<double> P(F0,F1,F2,F3);

  // Compute the curvature information
  P.ComputePrincipalCurvatures();

  // Dump out the principal curvature information
  // Expect 1.97909 and 0.106361 for principal curvatures
  // Expect {0.848674, 0.376491, 0.371492} and {-0.00124454, -0.700941, 0.713218}
  // as the principal directions.
  cout << "Curvatures are: " << P.kMax << " and " << P.kMin << endl;
  cout << "Princ. Vec. 1 : {" << P.tMax[0] << "," << P.tMax[1] << "," << P.tMax[2] << endl;
  cout << "Princ. Vec. 2 : {" << P.tMin[0] << "," << P.tMin[1] << "," << P.tMin[2] << endl;

  // Compute the ridge function at the point
  P.ComputeRidgeFunction();

  // Dump out the result (expect -0.219344)
  cout << "Ridge function: " << P.fRidge << endl;
}

void TestFastWindowedGaussianKernel3D()
{
  // Create a dummy image of random pixels
  typedef itk::Image<float, 3> ImageType;
  ImageType::Pointer img = ImageType::New();
  
  ImageType::SizeType sz = {{40,40,40}};
  ImageType::RegionType region(sz);
  img->SetRegions(region);
  img->Allocate();
  
  // Fill the image with randomness
  itk::ImageRegionIterator<ImageType> it(img, img->GetBufferedRegion());
  while(!it.IsAtEnd())
    {
    it.Set(rand() * 1.0 / RAND_MAX);
    ++it;
    }

  // Create a Gaussian kernel of radius 4
  double sigma[] = {1.0, 1.0, 1.0};
  FastGaussianKernel3D<float> KBad(img, 5, sigma);
  FastWindowedGaussianKernel3D<float> K(img, 5, sigma);

  // Transition accross a pixel boundary, making sure that everything is
  // legit
  double xStep = 0.0001, yStep = 0.0001, zStep = 0.0001;
  double xCenter = 20.0, yCenter = 17.0, zCenter = 22.0;
  for(int i=-5;i<=5;i++)
    {
    K.EvaluateAtPoint(xStep * i + xCenter, yStep * i + yCenter, zStep * i + zCenter);
    SpacePoint<double> p(K.G0, K.G1, K.G2, K.G3 );
    p.ComputePrincipalCurvatures();
    p.ComputeRidgeFunction();

    KBad.EvaluateAtPoint(xStep * i + xCenter, yStep * i + yCenter, zStep * i + zCenter);
    SpacePoint<double> pBad(KBad.G0, KBad.G1, KBad.G2, KBad.G3 );
    pBad.ComputePrincipalCurvatures();
    pBad.ComputeRidgeFunction();

    cout << "kGood[ " << i << "] = \t" << p.fRidge << "\t\t" << "kBad = " << pBad.fRidge << endl;
    }

  // Now test how many evaluations can be performed in a short time
  vector<double> v;
  for(unsigned int i = 0; i < 30000; i++)
    v.push_back(6.0 + rand() * 28.0 / RAND_MAX);
  
  bool flag = true;
  double c0 = clock();
  for(unsigned int j = 0; j < 30000; j+=3)
    {
    flag &= K.EvaluateAtPoint(v[j], v[j+1], v[j+2]);
    }
  double c1 = clock();
  cout << "10000 evaluations in " << (c1 - c0) / CLOCKS_PER_SEC << " seconds" << endl;
  cout << "Flag : " << flag << endl;
  cout << (int)(10000 * CLOCKS_PER_SEC / (c1 - c0)) << " evals per second" << endl;
}


class TestImplicitFunction01 : public ImplicitFunction3D {
public:
  bool EvaluateAtPoint(double x, double y, double z, int maxOrder = 3);
};

bool 
TestImplicitFunction01
::EvaluateAtPoint(double x, double y, double z, int maxOrder)
{
  G0 = -1. + 0.25*x*x + 0.1111111111111111*y*y + z*z;
  
  G1[0] = 0.5*x;
  G1[1] = 0.2222222222222222*y;
  G1[2] = 2*z;

  G2[0][0] = 0.5;
  G2[1][1] = 0.2222222222222222;
  G2[2][2] = 2.0;
  G2[0][1] = G2[1][2] = G2[0][2] = 0.0;

  G3[0][0][0] = G3[0][0][1] = G3[0][0][2] = G3[0][1][1] = G3[0][1][2] = 
    G3[0][2][2] = G3[1][1][1] = G3[1][1][2] = G3[1][2][2] = G3[2][2][2] = 0.0;

  return true;
}

class TestImplicitFunction02 : public ImplicitFunction3D {
public:
  bool EvaluateAtPoint(double x, double y, double z, int maxOrder = 3);
};

bool 
TestImplicitFunction02
::EvaluateAtPoint(double x, double y, double z, int maxOrder)
{
  G0 = -1. + sqrt(sqrt(pow(x,4) + pow(y,4)) + pow(z,2));
  
  G1[0] = pow(x,3)/(sqrt(pow(x,4) + pow(y,4))*sqrt(sqrt(pow(x,4) + pow(y,4)) + pow(z,2)));
  G1[1] = pow(y,3)/(sqrt(pow(x,4) + pow(y,4))*sqrt(sqrt(pow(x,4) + pow(y,4)) + pow(z,2)));
  G1[2] = z/sqrt(sqrt(pow(x,4) + pow(y,4)) + pow(z,2));

  G2[0][0] = (pow(x,2)*(pow(x,4)*pow(z,2) + 3*pow(y,4)*(sqrt(pow(x,4) + pow(y,4)) + pow(z,2))))/
   (pow(pow(x,4) + pow(y,4),1.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),1.5));
  
  G2[0][1] = -((pow(x,3)*pow(y,3)*(3*sqrt(pow(x,4) + pow(y,4)) + 2*pow(z,2)))/
     (pow(pow(x,4) + pow(y,4),1.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),1.5)));
  G2[0][2] = -((pow(x,3)*z)/(sqrt(pow(x,4) + pow(y,4))*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),1.5)));  
  G2[1][1] = (pow(y,2)*(pow(y,4)*pow(z,2) + 3*pow(x,4)*(sqrt(pow(x,4) + pow(y,4)) + pow(z,2))))/
   (pow(pow(x,4) + pow(y,4),1.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),1.5));
  G2[1][2] = -((pow(y,3)*z)/(sqrt(pow(x,4) + pow(y,4))*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),1.5)));
  G2[2][2] = sqrt(pow(x,4) + pow(y,4))/pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),1.5);

  G3[0][0][0] = (-3*x*(pow(x,8)*(5*pow(y,4) + sqrt(pow(x,4) + pow(y,4))*pow(z,2)) - 
      2*pow(y,8)*(pow(y,4) + 2*sqrt(pow(x,4) + pow(y,4))*pow(z,2) + pow(z,4)) + 
      pow(x,4)*pow(y,4)*(3*pow(y,4) + 7*sqrt(pow(x,4) + pow(y,4))*pow(z,2) + 2*pow(z,4))))/
    (pow(pow(x,4) + pow(y,4),2.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5));
  G3[0][0][1] =  (3*pow(x,2)*pow(y,3)*(4*pow(x,8) + pow(x,4)*(pow(y,4) + 5*sqrt(pow(x,4) + pow(y,4))*pow(z,2) + 2*pow(z,4)) - 
      pow(y,4)*(3*pow(y,4) + 5*sqrt(pow(x,4) + pow(y,4))*pow(z,2) + 2*pow(z,4))))/
    (pow(pow(x,4) + pow(y,4),2.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5));
  G3[0][0][2] =  -((pow(x,2)*z*(pow(x,4)*(-2*sqrt(pow(x,4) + pow(y,4)) + pow(z,2)) + 
        3*pow(y,4)*(sqrt(pow(x,4) + pow(y,4)) + pow(z,2))))/
    (pow(pow(x,4) + pow(y,4),1.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5)));
  G3[0][1][0] =  (3*pow(x,2)*pow(y,3)*(4*pow(x,8) + pow(x,4)*(pow(y,4) + 5*sqrt(pow(x,4) + pow(y,4))*pow(z,2) + 2*pow(z,4)) - 
      pow(y,4)*(3*pow(y,4) + 5*sqrt(pow(x,4) + pow(y,4))*pow(z,2) + 2*pow(z,4))))/
    (pow(pow(x,4) + pow(y,4),2.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5));
  G3[0][1][1] =  (3*pow(x,3)*pow(y,2)*(-3*pow(x,8) + pow(x,4)*(pow(y,4) - 5*sqrt(pow(x,4) + pow(y,4))*pow(z,2) - 2*pow(z,4)) + 
      pow(y,4)*(4*pow(y,4) + 5*sqrt(pow(x,4) + pow(y,4))*pow(z,2) + 2*pow(z,4))))/
    (pow(pow(x,4) + pow(y,4),2.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5));
  G3[0][1][2] =  (pow(x,3)*pow(y,3)*z*(5*sqrt(pow(x,4) + pow(y,4)) + 2*pow(z,2)))/
    (pow(pow(x,4) + pow(y,4),1.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5));
  G3[0][2][0] =  -((pow(x,2)*z*(pow(x,4)*(-2*sqrt(pow(x,4) + pow(y,4)) + pow(z,2)) + 
        3*pow(y,4)*(sqrt(pow(x,4) + pow(y,4)) + pow(z,2))))/
    (pow(pow(x,4) + pow(y,4),1.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5)));
  G3[0][2][1] = (pow(x,3)*pow(y,3)*z*(5*sqrt(pow(x,4) + pow(y,4)) + 2*pow(z,2)))/
    (pow(pow(x,4) + pow(y,4),1.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5));
  G3[0][2][2] = -((pow(x,3)*(sqrt(pow(x,4) + pow(y,4)) - 2*pow(z,2)))/
    (sqrt(pow(x,4) + pow(y,4))*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5)));
  G3[1][0][0] = (3*pow(x,2)*pow(y,3)*(4*pow(x,8) + pow(x,4)*(pow(y,4) + 5*sqrt(pow(x,4) + pow(y,4))*pow(z,2) + 2*pow(z,4)) - 
      pow(y,4)*(3*pow(y,4) + 5*sqrt(pow(x,4) + pow(y,4))*pow(z,2) + 2*pow(z,4))))/
    (pow(pow(x,4) + pow(y,4),2.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5));
  G3[1][0][1] = (3*pow(x,3)*pow(y,2)*(-3*pow(x,8) + pow(x,4)*(pow(y,4) - 5*sqrt(pow(x,4) + pow(y,4))*pow(z,2) - 2*pow(z,4)) + 
      pow(y,4)*(4*pow(y,4) + 5*sqrt(pow(x,4) + pow(y,4))*pow(z,2) + 2*pow(z,4))))/
    (pow(pow(x,4) + pow(y,4),2.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5));
  G3[1][0][2] =  (pow(x,3)*pow(y,3)*z*(5*sqrt(pow(x,4) + pow(y,4)) + 2*pow(z,2)))/
    (pow(pow(x,4) + pow(y,4),1.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5));
  G3[1][1][0] =  (3*pow(x,3)*pow(y,2)*(-3*pow(x,8) + pow(x,4)*(pow(y,4) - 5*sqrt(pow(x,4) + pow(y,4))*pow(z,2) - 2*pow(z,4)) + 
      pow(y,4)*(4*pow(y,4) + 5*sqrt(pow(x,4) + pow(y,4))*pow(z,2) + 2*pow(z,4))))/
    (pow(pow(x,4) + pow(y,4),2.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5));
  G3[1][1][1] = (-3*y*(-2*pow(x,12) + pow(y,8)*sqrt(pow(x,4) + pow(y,4))*pow(z,2) + 
      pow(x,4)*pow(y,4)*(5*pow(y,4) + 7*sqrt(pow(x,4) + pow(y,4))*pow(z,2) + 2*pow(z,4)) + 
      pow(x,8)*(3*pow(y,4) - 2*(2*sqrt(pow(x,4) + pow(y,4))*pow(z,2) + pow(z,4)))))/
    (pow(pow(x,4) + pow(y,4),2.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5));
  G3[1][1][2] = -((pow(y,2)*z*(pow(y,4)*(-2*sqrt(pow(x,4) + pow(y,4)) + pow(z,2)) + 
        3*pow(x,4)*(sqrt(pow(x,4) + pow(y,4)) + pow(z,2))))/
    (pow(pow(x,4) + pow(y,4),1.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5)));
  G3[1][2][0] = (pow(x,3)*pow(y,3)*z*(5*sqrt(pow(x,4) + pow(y,4)) + 2*pow(z,2)))/
    (pow(pow(x,4) + pow(y,4),1.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5));
  G3[1][2][1] = -((pow(y,2)*z*(pow(y,4)*(-2*sqrt(pow(x,4) + pow(y,4)) + pow(z,2)) + 
        3*pow(x,4)*(sqrt(pow(x,4) + pow(y,4)) + pow(z,2))))/
    (pow(pow(x,4) + pow(y,4),1.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5)));
  G3[1][2][2] = -((pow(y,3)*(sqrt(pow(x,4) + pow(y,4)) - 2*pow(z,2)))/
    (sqrt(pow(x,4) + pow(y,4))*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5)));
  G3[2][0][0] =  -((pow(x,2)*z*(pow(x,4)*(-2*sqrt(pow(x,4) + pow(y,4)) + pow(z,2)) + 
        3*pow(y,4)*(sqrt(pow(x,4) + pow(y,4)) + pow(z,2))))/
    (pow(pow(x,4) + pow(y,4),1.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5)));
  G3[2][0][1] = (pow(x,3)*pow(y,3)*z*(5*sqrt(pow(x,4) + pow(y,4)) + 2*pow(z,2)))/
    (pow(pow(x,4) + pow(y,4),1.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5));
  G3[2][0][2] =  -((pow(x,3)*(sqrt(pow(x,4) + pow(y,4)) - 2*pow(z,2)))/
    (sqrt(pow(x,4) + pow(y,4))*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5)));
  G3[2][1][0] = (pow(x,3)*pow(y,3)*z*(5*sqrt(pow(x,4) + pow(y,4)) + 2*pow(z,2)))/
    (pow(pow(x,4) + pow(y,4),1.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5));
  G3[2][1][1] =  -((pow(y,2)*z*(pow(y,4)*(-2*sqrt(pow(x,4) + pow(y,4)) + pow(z,2)) + 
        3*pow(x,4)*(sqrt(pow(x,4) + pow(y,4)) + pow(z,2))))/
    (pow(pow(x,4) + pow(y,4),1.5)*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5)));
  G3[2][1][2] =  -((pow(y,3)*(sqrt(pow(x,4) + pow(y,4)) - 2*pow(z,2)))/
    (sqrt(pow(x,4) + pow(y,4))*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5)));
  G3[2][2][0] =  -((pow(x,3)*(sqrt(pow(x,4) + pow(y,4)) - 2*pow(z,2)))/
    (sqrt(pow(x,4) + pow(y,4))*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5)));
  G3[2][2][1] =  -((pow(y,3)*(sqrt(pow(x,4) + pow(y,4)) - 2*pow(z,2)))/
    (sqrt(pow(x,4) + pow(y,4))*pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5)));
  G3[2][2][2] = (-3*sqrt(pow(x,4) + pow(y,4))*z)/pow(sqrt(pow(x,4) + pow(y,4)) + pow(z,2),2.5);

  return true;
}



void TestRidgeFunctions()
{
  TestImplicitFunction02 f;
  double x[9][3] = {
    {0.34, 0.55, -0.77},
    {-0.14, 0.568, -0.44},
    {0.87, -0.05, 2.19},
    {0.11, 0.11, 0.11},
    {-0.11, 0.11, 0.11},
    {0.11, -0.11, 0.11},
    {0.11, 0.11, -0.11},
    {0.0, 3.0, 0.0},
    {0.1, 3.0, 0.0}};
  
  for(unsigned int i = 0; i < 9; i++)
    {
    f.EvaluateAtPoint(x[i][0],x[i][1],x[i][2]);
    SpacePoint<double> p(f.G0, f.G1, f.G2, f.G3);
    p.ComputeRidgeFunction();
    cout << " k1 = " << p.kMax << ", k2 = " << p.kMin 
      << ", fMax = " << p.fRidgeMax << ", fMin = " << p.fRidgeMin << endl;
    }
}

int usage()
{
  cout << "PROGRAM: sulcalfilter" << endl;
  cout << "  This program refines a mesh based on the level set of a binary or" << endl;
  cout << "  a floating point image. It can either generate the mesh using a " << endl;
  cout << "  marching cubes variant, or it can refine an existing mesh by finding " << endl;
  cout << "  the ridges in the isosurface and including them in the mesh. A third " << endl;
  cout << "  mode of operation is to input a VTK mesh and to augment it with " << endl;
  cout << "  curvature and ridge information at the vertices" << endl;
  cout << "usage : " << endl;
  cout << "  sulcalfilter [options] image.hdr" << endl;
  cout << "options (mode selection) : " << endl;
  cout << "  -cubes out.vtk         Extract the isosurface and save it in VTK fmt. " << endl;
  cout << "  -ridge in.vtk out.vtk  Find ridges in input mesh " << endl;
  cout << "  -vgeo  in.vtk out.vtk  Augment mesh with vertex geometry " << endl;
  cout << "options (common to both modes): " << endl;
  cout << "  -s X.XX                Sigma (scale) of the Gaussian kernel in space units" << endl;
  cout << "  -r N                   Radius of the Gaussian kernel in pixels " << endl;
  cout << "  -i                     Invert the image before processing. [Out:+, In:-]" << endl;
  cout << "  -l X.XX                Level set function value (def. 0.0) " << endl;
  cout << "options (with -cubes only)" << endl;
  cout << "  -c X.XX                Size of the cube used for implicit surface computation" << endl;
  return -1;
}

// ===============================================================================
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/Types/PolyMesh_ArrayKernelT.hh>

// Define the mesh traits
class MyMeshTraits : public OpenMesh::DefaultTraits
{
public:
  // Use double vectors
  typedef OpenMesh::Vec3d Point;
  typedef OpenMesh::Vec3d Normal;

  // Allow storing normal along with the data (although can we trust this normal
  // as well as the one we get from the implicit function?)
  VertexAttributes ( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Status );

  // Store status for edges and faces because we need deletion
  FaceAttributes    ( OpenMesh::Attributes::Status );
  EdgeAttributes    ( OpenMesh::Attributes::Status );
  
  // Store the differential geometry along with the point
  VertexTraits {
  private:
    SpacePoint <double> p;
  public:
    VertexT() { }
    SpacePoint<double> & GetGeometry() { return p; }
    const SpacePoint<double>& GetGeometry() const { return p; }
  };
};

// Define the mesh type
typedef OpenMesh::PolyMesh_ArrayKernelT<MyMeshTraits>  MyMesh;

#include "implicit.h"

bool operator < (const itk::ContinuousIndex<double, 3> &c1, const itk::ContinuousIndex<double, 3> &c2)
{
  return c1[0] == c2[0] ? 
    (c1[1] == c2[1] ? c1[2] < c2[2] : c1[1] < c2[1]) : c1[0] < c2[0];
}

class ImplicitSurfaceExtractor {
public:
  typedef itk::Image<float, 3> ImageType;
  
  static void ExtractSurface(
    ImageType *img, ImplicitFunction3D *K, MyMesh *mesh, double cubeSize);

private:
  static ImplicitFunction3D *K;
  static MyMesh *mesh;
  static ImageType *image;

  static double ImplicitFunction(double x, double y, double z);

  static int TriangleCallback(int i1, int i2, int i3, IMP_VERTICES vInfo);

  // Convert from implicit's coordinates to image coordinates (why?)
  typedef itk::ContinuousIndex<double, 3> CIndex;
  static CIndex MapToImageSpace(double x, double y, double z) 
    {
    // Convert to pixel coorinates 
    itk::Point<double, 3> p; p[0] = x; p[1] = y; p[2] = z;
    CIndex index;
    image->TransformPhysicalPointToContinuousIndex(p, index);
    return index;
    }

  // Map from implicit's vertex indices to mesh handles
  typedef std::map<int, MyMesh::VertexHandle> HandleMapType; 
  static HandleMapType m_MapIndexToHandle;

  // A cache for holding previously computed differential geometry
  // typedef std::map< CIndex, SpacePoint<double> > CacheType;
  // static CacheType m_GeometryCache;

  // Counters for debugging
  static int nFunctionCalls;
  static int nFaceCalls;
  
};


double 
ImplicitSurfaceExtractor
::ImplicitFunction(double x, double y, double z)
{
  // Map the point to image space
  CIndex index = MapToImageSpace(x,y,z);
  
  // Evaluate the function
  K->EvaluateAtPoint(index[0], index[1], index[2], 0);

  // Track the number of callbacks
  // if(++nFunctionCalls % 10000 == 0)
  //   cout << "   ... " << nFunctionCalls << " function calls " << endl;

  // Return the result
  return K->G0;
}


/** This method creates a mesh from an implicit function */
void 
ImplicitSurfaceExtractor
::ExtractSurface(ImageType *img, ImplicitFunction3D *K, MyMesh *mesh, double cubeSize)
{
  // Store the kernel and the mesh
  ImplicitSurfaceExtractor::K = K;
  ImplicitSurfaceExtractor::mesh = mesh;
  ImplicitSurfaceExtractor::image = img;

  // We need to fit the image into a box {-1,1}^3 and choose a cube size
  // In order to do this we need to compute the extents of the image data. In 
  // addition, we will find the closest-to-the-level-set point
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> Iterator;
  Iterator it(img, img->GetBufferedRegion());

  // Get the image size
  ImageType::SizeType size = img->GetBufferedRegion().GetSize();

  // Record the corners of the bounding box as well as the point with smallest
  // function value. Inside is negative, outside is positive.
  typedef ImageType::IndexType Index;
  Index iLower = {{ size[0], size[1], size[2] }};
  Index iUpper = {{ 0, 0, 0 }};
  Index iNearest = {{ 0, 0, 0 }};
  float xNearest = 1.0e10f;

  while(!it.IsAtEnd()) 
    {
    // Check if the point is inside the surface
    if(it.Value() < 0.0f) 
      {
      Index i = it.GetIndex();
      for(unsigned int j = 0; j < 3; j++)
        {
        if(iLower[j] > i[j]) iLower[j] = i[j];
        if(iUpper[j] < i[j]) iUpper[j] = i[j];
        }

      // Check if this point is closest to nil
      if(- it.Value() < xNearest)
        {
        xNearest = -it.Value();
        iNearest = i;
        }
      }
    ++it;
    }

  // Print out the results of this computation
  cout << "  implicit surface is bounded by " << iLower << " and " << iUpper << endl;
  cout << "  good starting point is at " << iNearest << " with intensity " << xNearest << endl;

  // Find the position of the image center in space coordinates
  // Index iImageCenter = {{ size[0] / 2, size[1] / 2, size[2] / 2 }};
  
  // Map the bounding box to spatial coordinates
  itk::Point<double,3> pLower, pUpper, pNearest;
  img->TransformIndexToPhysicalPoint(iLower, pLower);
  img->TransformIndexToPhysicalPoint(iUpper, pUpper);
  img->TransformIndexToPhysicalPoint(iNearest, pNearest);

  // Compute the distance to each side of the bounding box, and take the max of
  // those six quantities.
  double q[6] = {
    pUpper[0] - pNearest[0],
    pUpper[1] - pNearest[1],
    pUpper[2] - pNearest[2],
    pNearest[0] - pLower[0],
    pNearest[1] - pLower[1],
    pNearest[2] - pLower[2]
  };

  double q0 = q[0] > q[1] ? q[0] : q[1];
  double q1 = q[2] > q[3] ? q[2] : q[3];
  double q2 = q[4] > q[5] ? q[4] : q[5];
  double qMax = q0 > q1 ? (q0 > q2 ? q0 : q2) : (q1 > q2 ? q1 : q2);  
  
  int qSize = (int) (qMax / cubeSize) + 4;

  cout << "  the cube range will be from " << -qSize << " to " << qSize << endl;
  cout << "  tracking progress ( each '*' represents 1000 triangles )" << endl << "  ";

  // Clear the index map
  m_MapIndexToHandle.clear();

  // Clear the geometry cache
  //m_GeometryCache.clear();
  
  // Clear the mesh
  mesh->clear();

  // Clear counters
  nFunctionCalls++;
  nFaceCalls++;
  
  // Run the implicit surface code from Blumenthal
  polygonize(
    &ImplicitSurfaceExtractor::ImplicitFunction,
    cubeSize,
    qSize,
    pNearest[0],
    pNearest[1],
    pNearest[2],
    &ImplicitSurfaceExtractor::TriangleCallback,
    TET);

  // Clear global variables
  ImplicitSurfaceExtractor::K = NULL;
  ImplicitSurfaceExtractor::mesh = NULL;
  ImplicitSurfaceExtractor::image = NULL;


  // Clear the progress
  cout << endl;
}

int ImplicitSurfaceExtractor
::TriangleCallback(int i1, int i2, int i3, IMP_VERTICES vInfo)
{
  // Use an array for easier access
  int idx[3] = {i1, i2, i3};

  // An array of handles for creating a face entry
  vector<MyMesh::VertexHandle> lHandles(3);

  // Loop over the vertices
  for(int i=0;i<3;i++)
    {
    // Check if each of the points already exists in the mesh
    map<int, MyMesh::VertexHandle>::const_iterator it = m_MapIndexToHandle.find(idx[i]);
    if(it != m_MapIndexToHandle.end())
      {
      // Just get the handle
      lHandles[i] = it->second;
      }
    else
      {
      // Create and add the vertex
      IMP_POINT *vertex = &vInfo.ptr[idx[i]].position;
      IMP_POINT *normal = &vInfo.ptr[idx[i]].normal;
      
      // Map to image space
      CIndex index = MapToImageSpace(vertex->x, vertex->y, vertex->z);
      
      // Get the handle
      lHandles[i] = mesh->add_vertex(MyMesh::Point(index[0], index[1], index[2]));

      // Evaluate the function
      K->EvaluateAtPoint(index[0], index[1], index[2]);

      // Construct the space point and save it
      SpacePoint<double> &P = mesh->vertex(lHandles[i]).GetGeometry();
      P.SetFunctionJet( K->G0, K->G1, K->G2, K->G3 );         

      // Set the normal       
      mesh->set_normal(lHandles[i], MyMesh::Normal( P.N ) );
      // mesh->set_normal(lHandles[i], MyMesh::Normal(normal->x, normal->y, normal->z));

      // Stick the handle in the map
      m_MapIndexToHandle[idx[i]] = lHandles[i];
      }
    }

  // Add the face to the mesh
  mesh->add_face(lHandles);

  // Track the number of callbacks
  if(++nFaceCalls % 1000 == 0)
    if(nFaceCalls % 70000 == 0) cout << endl << "  *";
    else cout << "*" << flush;

  // Continue
  return 1;
}

// Static stuff
MyMesh* ImplicitSurfaceExtractor::mesh = NULL;
ImplicitFunction3D* ImplicitSurfaceExtractor::K = NULL;
itk::Image<float, 3>* ImplicitSurfaceExtractor::image = NULL;
ImplicitSurfaceExtractor::HandleMapType ImplicitSurfaceExtractor::m_MapIndexToHandle;
int ImplicitSurfaceExtractor::nFunctionCalls = 0;
int ImplicitSurfaceExtractor::nFaceCalls = 0;

void ConvertOpenMeshToVTK(MyMesh *source, vtkPolyData *target)
{
  // Get the number of faces and vertices
  unsigned int nVertices = 0;
  MyMesh::VertexIter itVertex;
  for(itVertex = source->vertices_begin(); itVertex != source->vertices_end(); ++itVertex)
    nVertices++;

  unsigned int nFaces = 0;
  MyMesh::FaceIter itFace;
  for(itFace = source->faces_begin(); itFace != source->faces_end(); ++itFace)
    nFaces++;

  // Allocate points for the mesh
  vtkPoints *lPoints = vtkPoints::New();
  lPoints->Allocate(nFaces);
  
  // Create and array for storing the normal vectors
  vtkFloatArray *lNormals = vtkFloatArray::New();
  lNormals->SetNumberOfComponents(3);
  lNormals->SetNumberOfTuples(nVertices);

  // Allocate the points and cells for the VTK data
  target->Reset();
  target->Allocate(nFaces);

  // Put the points and normals into the Vtk mesh
  target->SetPoints(lPoints);
  target->GetPointData()->SetNormals(lNormals);

  // A map from handles to vertex ID's, just in case of an inconsistency
  typedef std::map< MyMesh::VertexHandle, vtkIdType> MyMap;
  MyMap handleMap;

  // Insert points one at a time
  for(itVertex = source->vertices_begin(); itVertex != source->vertices_end(); ++itVertex)
    {
    // Get point and normal info
    MyMesh::Point p = source->point(itVertex.handle());
    MyMesh::Normal n = source->normal(itVertex.handle());

    // Insert the point
    vtkIdType id = lPoints->InsertNextPoint(p[0], p[1], p[2]);

    // Insert the normal
    lNormals->SetTuple3(id, n[0], n[1], n[2]);

    // Store the ID
    handleMap[itVertex.handle()] = id;
    }

  // Insert faces one at a time
  for(itFace = source->faces_begin(); itFace != source->faces_end(); ++itFace)
    {
    vtkIdType vList[256], *vPtr = vList; 
    
    for(MyMesh::FaceVertexIter itLoop = source->fv_iter(itFace.handle());itLoop;++itLoop)
      {
      *vPtr++ = handleMap[itLoop.handle()];
      }

    if(vPtr - vList == 3)
      target->InsertNextCell(VTK_TRIANGLE, 3, vList);
    else
      target->InsertNextCell(VTK_POLYGON, vPtr - vList, vList);
    }
}

void ConvertVTKToOpenMesh(vtkPolyData *source, MyMesh *target)
{
  // Get the elements of the VTK mesh
  vtkPoints *lPoints = source->GetPoints();
  vtkDataArray *lNormals = source->GetPointData()->GetNormals();
  vtkCellArray *lPolys = source->GetPolys();

  // Create a point-id to handle map
  MyMesh::VertexHandle *idMap = new MyMesh::VertexHandle[lPoints->GetNumberOfPoints()];
  
  // Clear the target
  target->clear();

  // Add each of the vertices
  for(vtkIdType iVertex = 0; iVertex < lPoints->GetNumberOfPoints(); iVertex++)
    {
    // Get the point
    double *xPoint = lPoints->GetPoint(iVertex); 

    // Store the point in the mesh
    idMap[iVertex] = target->add_vertex( MyMesh::Point( xPoint[0], xPoint[1], xPoint[2] ));

    // Store the normal 
    if(lNormals && lNormals->GetNumberOfTuples() > iVertex) 
      {
      double *xNormal = lNormals->GetTuple(iVertex);
      target->set_normal( idMap[iVertex], 
        MyMesh::Normal( xNormal[0], xNormal[1], xNormal[2] )); 
      }
    }

  // Add each of the faces
  vtkIdType nCellPoints, *lCellPoints;
  for ( lPolys->InitTraversal(); lPolys->GetNextCell(nCellPoints,lCellPoints); ) 
    {
    // Vector to hold the faces
    vector<MyMesh::VertexHandle> lHandles(nCellPoints);

    // Add each point
    for( unsigned int iLoop = 0; iLoop < nCellPoints; iLoop++)
      {
      lHandles[iLoop] = idMap[ lCellPoints[iLoop] ];
      }

    // Add the face
    target->add_face(lHandles);
    }

  // Clean up
  delete idMap;
}

template<class TFunction>
void BrentRootSearch(
  TFunction &f, 
  double x1, double x2, 
  double f1, double f2,
  double &xRoot, double &fRoot,
  double tol = 1.0e-5)
{
  int iter, ITMAX = 100;
  const double EPS = 1.0e-8;

  double a=x1,b=x2,c=x2,d,e,min1,min2;
  double fa=f1;
  double fb=f2,fc,p,q,r,s,tol1,xm;

  fc=fb;
  for (iter=1;iter<=ITMAX;iter++) 
    {

    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) 
      {
      c=a;                           
      fc=fa;                       
      e=d=b-a;
      }

    if (fabs(fc) < fabs(fb)) 
      {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
      }

    tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol;     
    xm= 0.5 * ( c - b );

    if (fabs(xm) <= tol1 || fb == 0.0)
      { 
      xRoot = b;
      fRoot = fb;
      return;
      }

    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
      {
      s=fb/fa;                 
      if (a == c) 
        {
        p=2.0*xm*s;
        q=1.0-s;
        } 
      else 
        {
        q=fa/fc;
        r=fb/fc;
        p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
        }

      if (p > 0.0) q = -q;      
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);

      if (2.0*p < (min1 < min2 ? min1 : min2))
        {
        e=d;                      
        d=p/q;
        }
      else 
        {
        d=xm;                     
        e=d;
        }
      }
    else 
      {
      d=xm;                 
      e=d;
      }
    a=b;                   
    fa=fb;
    if (fabs(d) > tol1)      
      b += d;
    else 
      b += xm >= 0.0 ? fabs(tol1) : -fabs(tol1);
    fb=f(b);
    }
  
  assert(0);
  //   nrerror("Maximum number of iterations exceeded in zbrent");										    
  return; //0.0;            
}

class BrentTestFunction {
public:
  double operator() (double x) 
    {
    cout << "Testing location " << x << endl;
    return x * x * x - 3 * x * x + 2 * x - 1;
    }
};

void TestBrentSearch() 
{
  BrentTestFunction f;

  double x1 = 2.0, x2 = 3.0;
  double f1 = f(x1), f2 = f(x2);
  double xRoot, fRoot;

  BrentRootSearch(f,x1,x2,f1,f2,xRoot,fRoot);
  cout << "Root found at " << xRoot << " value " << fRoot << endl;
}

class MeshRidgeExtractor {
public:
  typedef MyMesh::Point MPoint;
  typedef SpacePoint<double> Geometry;

  void ExtractRidges(MyMesh *mesh, ImplicitFunction3D *K);

private:

  /** Class used in conjuction with Brent's search */
  class EdgeSearchFunction {
  public:
    EdgeSearchFunction(ImplicitFunction3D *K_) 
      : K(K_), nEvals(0) {};

    void SetPoints(MPoint &p1, MPoint &p2)
      { this->p1 = p1; this->p2 = p2; }

    // Used in ridge search
    double operator () ( double t )
      {
      pTest[0] = t * p2[0] + (1 - t) * p1[0];
      pTest[1] = t * p2[1] + (1 - t) * p1[1];
      pTest[2] = t * p2[2] + (1 - t) * p1[2];
      
      K->EvaluateAtPoint( pTest[0], pTest[1], pTest[2] );
      xGeometry.SetFunctionJet(K->G0, K->G1, K->G2, K->G3);
      xGeometry.ComputeRidgeFunction();

      // Print a star for every 10,000 evaluations
      if(++nEvals % 10000 == 0) 
        if(nEvals % 700000 == 0) cout << endl << "  *";
        else cout << "*" << flush;

      return xGeometry.fRidge;
      }  

    // Get the last geometry computation
    Geometry &GetLastGeometry() { return xGeometry; }

    // Get the last point tested
    MPoint &GetLastPoint() { return pTest; }

  private:
    // Initialization data
    ImplicitFunction3D *K;
    MPoint p1, p2, pTest;
    Geometry xGeometry;

    // Number of evaluations
    unsigned int nEvals;
  };

};

// This is the main method here. Find ridges in the implicit function
// defined by K and embed them as edges in the mesh. 
void 
MeshRidgeExtractor
::ExtractRidges(MyMesh *mesh, ImplicitFunction3D *K)
{
  // First of all, evaluate the ridge function at all points
  MyMesh::VertexIter itVertex = mesh->vertices_begin();
  while( itVertex != mesh->vertices_end() )
    {
    // Evaluate the kernel at the point
    MyMesh::Point p = mesh->point( itVertex );
    K->EvaluateAtPoint(p[0], p[1], p[2]);

    // Assign the geometry
    Geometry &g = mesh->vertex( itVertex ).GetGeometry();
    g.SetFunctionJet( K->G0, K->G1, K->G2, K->G3 );
    g.ComputePrincipalCurvatures();
    g.ComputeRidgeFunction();

    // On to the next vertex
    ++itVertex;
    }

  // Function used to search the edges for ridge crossings
  EdgeSearchFunction f( K );

  // We begin by finding all edges in the mesh that happen to cross a ridge, and 
  // creating the vertex that lies on the ridge. However, we do not yet split the
  // edges. Instead, we record edge subdivision rules, which we apply later on. 
  vector<MyMesh::VertexHandle> lToDoStart, lToDoNew, lToDoEnd;
  MyMesh::EdgeIter itEdge = mesh->edges_begin();
  while(itEdge != mesh->edges_end())
    {
    // Check whether the two vertices on the opposite sides of the edge have
    // different signs
    MyMesh::HalfedgeHandle heh0 = mesh->halfedge_handle( itEdge, 0 );
    MyMesh::VertexHandle vh0 = mesh->to_vertex_handle( heh0 );
    MyMesh::VertexHandle vh1 = mesh->from_vertex_handle( heh0 );
    MyMesh::Vertex &v0 = mesh->vertex( vh0 );
    MyMesh::Vertex &v1 = mesh->vertex( vh1 );

    // Check the signs
    double r0 = v0.GetGeometry().fRidge;
    double r1 = v1.GetGeometry().fRidge;
    if(r0 < 0.0 && r1 > 0.0 || r0 > 0.0 && r1 < 0.0) 
      {
      // Set up a functor to pass to the Brent's algorithm
      f.SetPoints(mesh->point(vh0), mesh->point(vh1));

      // Perform the Brent search, pass dummy variables for results
      double xRoot = 0.0, fRoot = 0.0;
      BrentRootSearch( f, 0.0, 1.0, r0, r1, xRoot, fRoot, 1.0e-6 );

      // Add a new vertex to the mesh
      MyMesh::VertexHandle vhNew = mesh->add_vertex( f.GetLastPoint() );

      // Copy the geometry to the new point, and set its normal
      mesh->vertex( vhNew ).GetGeometry() = f.GetLastGeometry();
      mesh->set_normal( vhNew, MyMesh::Normal(f.GetLastGeometry().N) ); 

      // Store the start, end and middle of the insertion that we must make
      lToDoStart.push_back(vh0);
      lToDoNew.push_back(vhNew);
      lToDoEnd.push_back(vh1);
      }

    ++itEdge;
    }

  // Now, we apply all the edge subdivision rules by deleting edges and inserting new
  // faces into the mesh. This is probably not the most efficient approach.
  for(unsigned int iInsert = 0; iInsert < lToDoStart.size(); iInsert++)
    {
    // Get the handles
    MyMesh::VertexHandle vh0 = lToDoStart[iInsert];
    MyMesh::VertexHandle vh1 = lToDoEnd[iInsert];
    MyMesh::VertexHandle vhNew = lToDoNew[iInsert];

    // Find the edge leading from one point to the other
    MyMesh::HalfedgeHandle heh0 = mesh->find_halfedge( vh0, vh1 );
    assert(heh0.is_valid());

    // Find out which two faces are adjacent to the edge
    MyMesh::HalfedgeHandle heh1 = mesh->opposite_halfedge_handle(heh0);
    MyMesh::FaceHandle fh0 = mesh->face_handle( heh0 );
    MyMesh::FaceHandle fh1 = mesh->face_handle( heh1 );

    // Create vectors that would hold the vertices in the new, expanded faces
    vector<MyMesh::VertexHandle> vl0, vl1;
    vl0.push_back( vhNew );
    vl1.push_back( vhNew );

    // The vertex list vl0 needs to contain vertices vNew, v0, ..., v1
    MyMesh::FaceVertexIter it = mesh->fv_iter(fh0);
    while( it.handle() != vh1 ) ++it;
    do 
      {
      vl0.push_back( it.handle() );
      ++it;
      } 
    while ( it.handle() != vh1 );
      
    // The vertex list vl1 needs to contain vertices vNew, v1, ..., v0
    it = mesh->fv_iter( fh1 );
    while( it.handle() != vh0 ) ++it;
    do 
      {
      vl1.push_back( it.handle() );
      ++it;
      } 
    while ( it.handle() != vh0 );

    // Delete the split edge
    mesh->delete_edge( mesh->edge_handle( heh0 ) );

    // Add the newly created faces to the mesh
    mesh->add_face( vl0 );
    mesh->add_face( vl1 );
    }

  // Garbage collect the mesh
  mesh->garbage_collection();

}

void AugmentVTKMeshVertices( vtkPolyData *meshVTK, ImplicitFunction3D *K )
{
  // Create data arrays for curvatures
  const unsigned int nArrays = 10;
  vtkFloatArray *xCrvArray[nArrays];
  const char *names[] = 
    { "MAJOR", "MINOR", "GAUSS", "MEAN", 
      "KOEN_C", "KOEN_S", "RIDGE_MAJOR", "RIDGE_MINOR",
      "RIDGE_MAIN", "RIDGE_BOTH"};

  // Initialize the arrays 
  unsigned int i;
  for(i=0; i < nArrays; i++) 
    {
    xCrvArray[i] = vtkFloatArray::New();
    xCrvArray[i]->SetName(names[i]);
    }

  // Iterate through the meshVTK
  for(vtkIdType iVertex = 0; iVertex < meshVTK->GetNumberOfPoints(); iVertex++)
    {
    // Compute the jet at each vertex
    double *x = meshVTK->GetPoint(iVertex);

    // A space point to store the results of all evaluations
    SpacePoint<double> P;

    if( K->EvaluateAtPoint(x[0],x[1],x[2]) )
      { P.SetFunctionJet(K->G0, K->G1, K->G2, K->G3); }

    P.ComputePrincipalCurvatures();
    P.ComputeRidgeFunction();

    // Assign curvatures to data arrays 
    xCrvArray[0]->InsertNextValue(P.kMax);
    xCrvArray[1]->InsertNextValue(P.kMin);
    xCrvArray[2]->InsertNextValue(P.kGauss);
    xCrvArray[3]->InsertNextValue(P.kMean);

    // Koenderink coloring. Metric C measures how curved a surface is, with zero representing
    // 'unit' curvature -inf being very flat and +inf being very curvy. Metric S describes the 
    // local shape, with -1 being a concave umbillic, to saddle point to convex umbillic + 1
    const float xLogBoundLow = exp(-4.0f);
    const float xLogBoundHi = exp(4.0f);
    float xKoenderinkC = P.kMax * P.kMax + P.kMin * P.kMin;
    float xKoenderinkCLog = 
      xKoenderinkC > xLogBoundHi ? 4.0f : 
      ( xKoenderinkC < xLogBoundLow ? -4.0f : log(xKoenderinkC) );

    float xKoenderinkS = atan2(P.kMax + P.kMin, P.kMax - P.kMin) / MY_PI;

    // Store the Koenderink values
    xCrvArray[4]->InsertNextValue( 0.25 * xKoenderinkCLog );
    xCrvArray[5]->InsertNextValue( xKoenderinkS );

    // Store the ridge and the implicit function values
    xCrvArray[6]->InsertNextValue( P.fRidgeMax );
    xCrvArray[7]->InsertNextValue( P.fRidgeMin );
    xCrvArray[8]->InsertNextValue( P.fRidgeDominant );
    xCrvArray[9]->InsertNextValue( P.fRidgeProduct );
    }

  // Store the arrays with the data
  for(i=0; i < nArrays; i++) 
    meshVTK->GetPointData()->AddArray(xCrvArray[i]);
}

void SaveVTKMesh(vtkPolyData *meshVTK, const char *file)
{
  // Dump out mesh information
  cout << "Saving mesh " << file << endl;
  cout << "  Mesh Summary " << endl;
  cout << "    Vertices:    \t\t" << meshVTK->GetNumberOfPoints() << endl;
  cout << "    Cells:       \t\t" << meshVTK->GetNumberOfCells() << endl;

  // Count the number of triangles, quadrangles, etc.
  int nPolys[] = {0,0,0,0,0,0,0,0}, nPolysMax = 7;
  for(vtkIdType id = 0; id < meshVTK->GetNumberOfCells(); id++)
    {
    vtkCell *cell = meshVTK->GetCell(id);
    int np = cell->GetNumberOfPoints();
    if(np < nPolysMax)
      nPolys[np]++;
    else
      nPolys[nPolysMax]++;
    }

  for(unsigned int i = 0; i < nPolysMax; i++)
    {
    if(nPolys[i] > 0)
      cout << "    Cells w/ " << i << " verts: \t\t" << nPolys[i] << endl;
    }
  if(nPolys[nPolysMax] > 0)
    cout << "    Cells w/ many verts: \t" << nPolys[nPolysMax] << endl;

  // Save the mesh
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(meshVTK);
  writer->SetFileName(file);
  writer->SetFileTypeToBinary();
  writer->Update();
  writer->Delete();
}

int main(int argc, char *argv[])
{
  // TestFastGaussianKernel1D();
  // TestFastWindowedGaussianKernel1D();
  // TestFastWindowedGaussianKernel3D();
  // TestFastGaussianKernel3D();
  // TestCurvatureComputation();
  // TestBrentSearch();
  // TestRidgeFunctions();

  // Usage
  if(argc < 3) return usage();

  // Process the options
  const char *fnInput = argv[argc-1];

  // Mesh input and output files
  const char *fnRidgeModeInputMesh = "";
  const char *fnRidgeModeOutputMesh = "";
  const char *fnCubesModeOutputMesh = "";
  const char *fnVGModeInputMesh = "";
  const char *fnVGModeOutputMesh = "";

  // The mode in which to run the program
  bool flagModeCubes = false, flagModeRidge = false, flagModeVertexGeometry = false;

  // Parameters of the algorithm
  double sigma = 1.0, cubeSize = 1.0, levelSet = 0.0, flipScale = 1.0;
  int radius = 6;

  // Parse the command line arguments
  for(unsigned int iArg = 1; iArg < argc - 1; iArg++)
    {
    if(!strcmp(argv[iArg], "-s"))
      {
      sigma = atof(argv[++iArg]);
      }
    else if(!strcmp(argv[iArg], "-r"))
      {
      radius = atoi(argv[++iArg]);
      }
    else if(!strcmp(argv[iArg],"-i"))
      {
      flipScale = -1.0;
      }
    else if(!strcmp(argv[iArg],"-c"))
      {
      cubeSize = atof(argv[++iArg]);
      }
    else if(!strcmp(argv[iArg],"-l"))
      {
      levelSet = atof(argv[++iArg]);
      }
    else if(!strcmp(argv[iArg],"-cubes"))
      {
      flagModeCubes = true;
      fnCubesModeOutputMesh = argv[++iArg];
      }
    else if(!strcmp(argv[iArg],"-vgeo"))
      {
      flagModeVertexGeometry = true;
      fnVGModeInputMesh = argv[++iArg];
      fnVGModeOutputMesh = argv[++iArg];
      }
    else if(!strcmp(argv[iArg],"-ridge"))
      {
      flagModeRidge = true;
      fnRidgeModeInputMesh = argv[++iArg];
      fnRidgeModeOutputMesh = argv[++iArg];
      }
    }

  // Update the user as to what we are doing
  cout << "==============================" << endl;
  cout << "Running SulcalFilter " << endl;
  cout << "  input image:             " << fnInput << endl;
  cout << "  kernel sigma:            " << sigma << endl;
  cout << "  kernel radius:           " << radius << endl;
  cout << "  flip enabled:            " << (flipScale < 0) << endl;
  cout << "  level set value:         " << levelSet << endl;
  if( flagModeCubes )
    {
    cout << "  iso-surface extraction:  ON" << endl;
    cout << "    output mesh:           " << fnCubesModeOutputMesh << endl;
    cout << "    marching cube size:    " << cubeSize << endl;
    }
  if( flagModeRidge )
    {
    cout << "  ridge extraction:        ON" << endl;
    cout << "    input mesh:            " << fnRidgeModeInputMesh << endl;
    cout << "    output mesh:           " << fnRidgeModeOutputMesh << endl;
    }
  if( flagModeVertexGeometry )  
    {
    cout << "  vertex geom. augment-n:  ON" << endl;
    cout << "    input mesh:            " << fnVGModeInputMesh << endl;
    cout << "    output mesh:           " << fnVGModeOutputMesh << endl;
    }
  cout << "------------------------------" << endl;

  // Define the float based image
  typedef itk::Image<float,3> ImageType;
  ImageType::Pointer imgInput;

  // Read the image
  ReadImage(imgInput,fnInput);

  // Flip the image if necessary, to make the inside negative and the outside positive
  if(flipScale < 0.0 || levelSet != 0.0) 
    {
    itk::ImageRegionIterator<ImageType> it(imgInput,imgInput->GetBufferedRegion());
    while(!it.IsAtEnd())
      {
      it.Set( - (it.Get() - levelSet) );
      ++it;
      }
    }

  // Create the xyz sigmas
  double s3[3];
  s3[0] = sigma / imgInput->GetSpacing()[0];
  s3[1] = sigma / imgInput->GetSpacing()[1];
  s3[2] = sigma / imgInput->GetSpacing()[2];

  // Create the Gaussian filter
  FastWindowedGaussianKernel3D<float> K(imgInput, radius, s3);

  // From here, we split, depending on what the user wants to do
  if(flagModeCubes) 
    {
    // Keep the user informed
    cout << "Extracting the iso-surface. " << endl;

    // Create a meshVTK
    MyMesh mymesh;

    // Compute implicit surface
    ImplicitSurfaceExtractor::ExtractSurface(imgInput, &K, &mymesh, cubeSize);

    // Convert to a VTK meshVTK
    vtkPolyData *meshVTK = vtkPolyData::New();
    ConvertOpenMeshToVTK(&mymesh,meshVTK);

    // Save the array
    SaveVTKMesh( meshVTK, fnCubesModeOutputMesh);
    }

  // The ridge mode
  if(flagModeRidge) 
    {
    // Keep the user informed
    cout 
      << "Augmenting mesh with ridges. Each '*' represents 10000 evaluations." 
      << endl << "  ";

    // Load the input mesh
    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    reader->SetFileName(fnRidgeModeInputMesh);
    reader->Update();

    // Get the polygon data
    vtkPolyData *meshVTK = reader->GetOutput();

    // Convert the mesh to OpenMesh
    MyMesh mymesh;
    ConvertVTKToOpenMesh(meshVTK, &mymesh);

    // Create the ridge extractor
    MeshRidgeExtractor mre;
    mre.ExtractRidges(&mymesh, &K);
    cout << endl;

    // Convert the mesh back to VTK 
    ConvertOpenMeshToVTK(&mymesh, meshVTK);

    // Save the mesh
    SaveVTKMesh(meshVTK, fnRidgeModeOutputMesh);
    }

  // The ridge mode
  if(flagModeVertexGeometry) 
    {
    // Keep the user informed
    cout << "Augmenting mesh vertices with curvature info" << endl;

    // Load the input mesh
    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    reader->SetFileName( fnVGModeInputMesh );
    reader->Update();

    // Get the polygon data
    vtkPolyData *meshVTK = reader->GetOutput();

    // Augment the mesh
    AugmentVTKMeshVertices(meshVTK, &K);

    // Save the mesh
    SaveVTKMesh( meshVTK, fnVGModeOutputMesh );
    }
}
