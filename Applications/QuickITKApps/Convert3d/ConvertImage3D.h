#ifndef __ConvertImage3D_h_
#define __ConvertImage3D_h_

// We define manual instantiation to speed up compilation
#define ITK_MANUAL_INSTANTIATION 1

// ITK includes
#include "itkAddImageFilter.h"
#include "itkAntiAliasBinaryImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkByteSwapper.h"
#include "itkCastImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkIOCommon.h"
#include "itkLinearInterpolateImageFunction.h" 
#include "itkMetaDataObject.h"
#include "itkMultiplyImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkPovRayDF3ImageIOFactory.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkResampleImageFilter.h"
#include <itkSegmentationLevelSetImageFilter.h>
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkSTAPLEImageFilter.h"
#include <itksys/SystemTools.hxx>
#include "itkVnlFFTRealToComplexConjugateImageFilter.h"
#include "itkVoxBoCUBImageIOFactory.h"

#include <iostream>
#include <cctype>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>

using namespace std;

template<class TPixel, unsigned int VDim>
class ImageConverter
{
public:
  ImageConverter();
  int ProcessCommandLine(int argc, char *argv[]);

private:
  // Image typedef
  typedef itk::Image<TPixel, VDim> ImageType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef typename ImageType::SizeType SizeType;
  typedef typename ImageType::IndexType IndexType;
  typedef typename ImageType::RegionType RegionType;
  typedef vnl_vector_fixed<double, VDim> RealVector;
  typedef vnl_vector_fixed<int, VDim> IntegerVector;

  // Complex stuff
  typedef std::complex<TPixel> ComplexPixel;
  typedef itk::Image<ComplexPixel, VDim> ComplexImageType;

  // Iterators
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> ConstIterator;
  
  // Internal functions
  void ReadImage(const char *file);
  void WriteImage(const char *file, bool force);
  void AntiAliasImage(double iso);
  void AddImages();
  void CopyImage();
  void ComputeFFT();
  void ComputeOverlaps(double value);
  void ConnectedComponents();
  void ImageERF(double thresh, double scale);
  void SignedDistanceTransform();
  void ExtractRegion(RegionType bbox);
  void LevelSetSegmentation(int nIter);
  void MultiplyImages();
  void SampleImage(const RealVector &x);
  void ScaleShiftImage(double a, double b);
  void StapleAlgorithm(double ival);
  void PrintImageInfo(bool flagFullInfo);
  void ReplaceIntensities(vector<double> &rules);
  void ThresholdImage(double u1, double u2, double vIn, double vOut);
  void TrimImage(const RealVector &margin);
  int ProcessCommand(int argc, char *argv[]);
  int ProcessResampleCommand(int argc, char *argv[]);
  int ProcessSmoothCommand(int argc, char *argv[]);

  // Read vectors, etc from command line
  SizeType ReadSizeVector(char *vec);
  IndexType ReadIndexVector(char *vec);
  RealVector ReadRealVector(char *vec);

  // Get bounding box of an image
  void GetBoundingBox(ImageType *image, RealVector &bb0, RealVector &bb1);

  // Templated write function
  template<class TOutPixel> void TemplatedWriteImage(const char *file, double xRoundFactor);

  // Stack of images from the command line
  vector<ImagePointer> m_ImageStack;

  // Typeid of the image to be saved
  string m_TypeId;

  // Interpolation type
  string m_Interpolation;

  // Background intensity
  double m_Background;

  // Rounding factor (not used for float and double) is 0.5 unless -noround was specified
  double m_RoundFactor;

  // Whether SPM extensions are used
  bool m_FlagSPM;

  // Number of iterations for various algorithms
  size_t m_Iterations;

  // Root mean square error for anti-aliasing algorithm
  double m_AntiAliasRMS;

  // Level set algorithm parameters
  double m_LevSetCurvature, m_LevSetAdvection;

  // Verbose output stream
  std::ostringstream devnull;
  std::ostream *verbose;
};

#endif

