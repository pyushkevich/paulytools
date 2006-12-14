#undef ITK_MANUAL_INSTANTIATION

#include "itkImage.h"
#include "itkImageRegion.h"
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
#include "itkAntiAliasBinaryImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkVoxBoCUBImageIOFactory.h"
#include "itkPovRayDF3ImageIOFactory.h"
#include "itkArray.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkVnlFFTRealToComplexConjugateImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"

#include <string>

typedef itk::Image<double, 3> DoubleImage;
typedef itk::Image<std::complex<double>, 3> ComplexImage;

// Double image stuff
template class itk::Image<double, 3>;
template class itk::ImageRegion<3>;
template class itk::ImageFileReader<DoubleImage>;
template class itk::ImageFileWriter<DoubleImage>;
template class itk::ResampleImageFilter<DoubleImage,DoubleImage>;
template class itk::AntiAliasBinaryImageFilter<DoubleImage,DoubleImage>;
template class itk::DiscreteGaussianImageFilter<DoubleImage,DoubleImage>;
template class itk::ImageToImageFilter<DoubleImage,DoubleImage>;
template class itk::ShiftScaleImageFilter<DoubleImage,DoubleImage>;
template class itk::ExtractImageFilter<DoubleImage, DoubleImage>;
template class itk::FixedArray<double,3>;
template class itk::Array2D<double>;
template class itk::LinearInterpolateImageFunction<DoubleImage, double>;
template class itk::NearestNeighborInterpolateImageFunction<DoubleImage>;
template class itk::BSplineInterpolateImageFunction<DoubleImage>;
template class itk::MetaDataObject<std::string>;
template class itk::Transform<double, 3, 3>;
template class itk::BinaryThresholdImageFilter<DoubleImage, DoubleImage>;
template class itk::ImageRegionIteratorWithIndex<DoubleImage>;
template class itk::ImageRegionConstIteratorWithIndex<DoubleImage>;
template class itk::Array<double>;

// All the junk we need for complex
template class itk::VnlFFTRealToComplexConjugateImageFilter<double, 3>;
template class itk::ComplexToRealImageFilter<ComplexImage, DoubleImage>;
template class itk::ComplexToImaginaryImageFilter<ComplexImage, DoubleImage>;
template class itk::ImageToImageFilter<ComplexImage,DoubleImage>;
template class itk::ImageToImageFilter<DoubleImage,ComplexImage>;
template class itk::Image<std::complex<double>, 3>;
template class itk::ImageSource<ComplexImage>;
template class itk::InPlaceImageFilter<ComplexImage, DoubleImage>;

// Stuff associated with writing images
template class itk::Image<unsigned char,3>;
template class itk::Image<char,3>;
template class itk::Image<unsigned short,3>;
template class itk::Image<short,3>;
template class itk::Image<unsigned int,3>;
template class itk::Image<int,3>;
template class itk::Image<float,3>;

template class itk::ImageFileWriter<itk::Image<unsigned char,3> >;
template class itk::ImageFileWriter<itk::Image<char,3> >;
template class itk::ImageFileWriter<itk::Image<unsigned short,3> >;
template class itk::ImageFileWriter<itk::Image<short,3> >;
template class itk::ImageFileWriter<itk::Image<unsigned int,3> >;
template class itk::ImageFileWriter<itk::Image<int,3> >;
template class itk::ImageFileWriter<itk::Image<float,3> >;
