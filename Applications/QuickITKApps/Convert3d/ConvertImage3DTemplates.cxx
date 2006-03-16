#undef ITK_MANUAL_INSTANTIATION

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkExtractImageFilter.h"
#include "itkByteSwapper.h"
#include "itkMetaDataObject.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h" 
#include "itkDiscreteGaussianImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkVoxBoCUBImageIOFactory.h"
#include "itkPovRayDF3ImageIOFactory.h"

#include <string>

typedef itk::Image<double, 3> DoubleImage;

// Double image stuff
template DoubleImage;
template itk::ImageFileReader<DoubleImage>;
template itk::ImageFileWriter<DoubleImage>;
template itk::ResampleImageFilter<DoubleImage,DoubleImage>;
template itk::DiscreteGaussianImageFilter<DoubleImage,DoubleImage>;
template itk::ImageToImageFilter<DoubleImage,DoubleImage>;
template itk::ShiftScaleImageFilter<DoubleImage,DoubleImage>;
template itk::FixedArray<double,3>;
template itk::Array2D<double>;
template itk::LinearInterpolateImageFunction<DoubleImage, double>;
template itk::NearestNeighborInterpolateImageFunction<DoubleImage>;
template itk::BSplineInterpolateImageFunction<DoubleImage>;
template itk::MetaDataObject<std::string>;
template itk::Transform<double, 3, 3>;
template itk::BinaryThresholdImageFilter<DoubleImage, DoubleImage>;


// Stuff associated with writing images
template itk::Image<unsigned char,3>;
template itk::Image<char,3>;
template itk::Image<unsigned short,3>;
template itk::Image<short,3>;
template itk::Image<unsigned int,3>;
template itk::Image<int,3>;
template itk::Image<float,3>;

template itk::ImageFileWriter<itk::Image<unsigned char,3> >;
template itk::ImageFileWriter<itk::Image<char,3> >;
template itk::ImageFileWriter<itk::Image<unsigned short,3> >;
template itk::ImageFileWriter<itk::Image<short,3> >;
template itk::ImageFileWriter<itk::Image<unsigned int,3> >;
template itk::ImageFileWriter<itk::Image<int,3> >;
template itk::ImageFileWriter<itk::Image<float,3> >;
