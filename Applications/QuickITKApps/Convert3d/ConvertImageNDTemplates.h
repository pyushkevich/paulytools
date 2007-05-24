#ifndef __ConvertImage3DTemplates_h_
#define __ConvertImage3DTemplates_h_

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
#include "itkVectorLinearInterpolateImageFunction.h" 
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
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkCastImageFilter.h"
#include <itkSegmentationLevelSetImageFilter.h>
#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"

#include <string>

typedef itk::Image<double, C3D_VDIM> DoubleImage;
typedef itk::Image<int, C3D_VDIM> IntImage;
typedef itk::Image<std::complex<double>, C3D_VDIM> ComplexImage;
typedef itk::Image<itk::FixedArray<double, C3D_VDIM>, C3D_VDIM> VectorImage;

// Double image stuff
template class itk::Image<double, C3D_VDIM>;
template class itk::ImageRegion<C3D_VDIM>;
template class itk::ImageFileReader<DoubleImage>;
template class itk::ImageFileWriter<DoubleImage>;
template class itk::ResampleImageFilter<DoubleImage,DoubleImage>;
template class itk::AntiAliasBinaryImageFilter<DoubleImage,DoubleImage>;
template class itk::DiscreteGaussianImageFilter<DoubleImage,DoubleImage>;
template class itk::ImageToImageFilter<DoubleImage,DoubleImage>;
template class itk::ShiftScaleImageFilter<DoubleImage,DoubleImage>;
template class itk::ExtractImageFilter<DoubleImage, DoubleImage>;
template class itk::FixedArray<double,C3D_VDIM>;
template class itk::Array2D<double>;
template class itk::LinearInterpolateImageFunction<DoubleImage, double>;
template class itk::LinearInterpolateImageFunction<DoubleImage, float>;
template class itk::VectorLinearInterpolateImageFunction<VectorImage, double>;
template class itk::VectorLinearInterpolateImageFunction<VectorImage, float>;
template class itk::NearestNeighborInterpolateImageFunction<DoubleImage>;
template class itk::BSplineInterpolateImageFunction<DoubleImage>;
template class itk::MetaDataObject<std::string>;
template class itk::Transform<double, C3D_VDIM, C3D_VDIM>;
template class itk::BinaryThresholdImageFilter<DoubleImage, DoubleImage>;
template class itk::ImageRegionIteratorWithIndex<DoubleImage>;
template class itk::ImageRegionConstIteratorWithIndex<DoubleImage>;
template class itk::Array<double>;
template class itk::SegmentationLevelSetImageFilter<DoubleImage, DoubleImage, double>;
template class itk::SegmentationLevelSetFunction<DoubleImage, DoubleImage>;
template class itk::ConnectedComponentImageFilter<DoubleImage, IntImage>;
template class itk::RelabelComponentImageFilter<IntImage, IntImage>;
template class itk::UnaryFunctorImageFilter<IntImage, DoubleImage, 
  itk::Functor::Cast<int, double> >;
template class itk::UnaryFunctorImageFilter<ComplexImage, DoubleImage, 
  itk::Function::ComplexToImaginary<std::complex<double>, double> >;
template class itk::UnaryFunctorImageFilter<ComplexImage, DoubleImage, 
  itk::Function::ComplexToReal<std::complex<double>, double> >;
template class itk::SparseFieldLevelSetImageFilter<DoubleImage, DoubleImage>;
template class itk::AddImageFilter<DoubleImage, DoubleImage, DoubleImage>;
template class itk::MultiplyImageFilter<DoubleImage, DoubleImage, DoubleImage>;
template class itk::BinaryFunctorImageFilter<DoubleImage, DoubleImage, DoubleImage,
  itk::Functor::Add2<double,double,double> >;
template class itk::BinaryFunctorImageFilter<DoubleImage, DoubleImage, DoubleImage,
  itk::Function::Mult<double,double,double> >;

// All the junk we need for complex
template class itk::VnlFFTRealToComplexConjugateImageFilter<double, C3D_VDIM>;
template class itk::ComplexToRealImageFilter<ComplexImage, DoubleImage>;
template class itk::ComplexToImaginaryImageFilter<ComplexImage, DoubleImage>;
template class itk::ImageToImageFilter<ComplexImage,DoubleImage>;
template class itk::ImageToImageFilter<DoubleImage,ComplexImage>;
template class itk::Image<std::complex<double>, C3D_VDIM>;
template class itk::ImageSource<ComplexImage>;
template class itk::InPlaceImageFilter<ComplexImage, DoubleImage>;

// All the junk we need for integer images
template class itk::CastImageFilter<IntImage, DoubleImage>;
template class itk::ImageSource<IntImage>;
template class itk::ImageToImageFilter<DoubleImage, IntImage>;
template class itk::ImageToImageFilter<IntImage, DoubleImage>;
template class itk::InPlaceImageFilter<IntImage, DoubleImage>;
template class itk::InPlaceImageFilter<IntImage, IntImage>;

// Stuff associated with writing images
template class itk::Image<unsigned char,C3D_VDIM>;
template class itk::Image<char,C3D_VDIM>;
template class itk::Image<unsigned short,C3D_VDIM>;
template class itk::Image<short,C3D_VDIM>;
template class itk::Image<unsigned int,C3D_VDIM>;
template class itk::Image<int,C3D_VDIM>;
template class itk::Image<float,C3D_VDIM>;

template class itk::ImageFileWriter<itk::Image<unsigned char,C3D_VDIM> >;
template class itk::ImageFileWriter<itk::Image<char,C3D_VDIM> >;
template class itk::ImageFileWriter<itk::Image<unsigned short,C3D_VDIM> >;
template class itk::ImageFileWriter<itk::Image<short,C3D_VDIM> >;
template class itk::ImageFileWriter<itk::Image<unsigned int,C3D_VDIM> >;
template class itk::ImageFileWriter<itk::Image<int,C3D_VDIM> >;
template class itk::ImageFileWriter<itk::Image<float,C3D_VDIM> >;

#endif
