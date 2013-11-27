#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkMetaDataObject.h"
#include <iostream>

int usage()
{
  std::cout << "shrink_warp : a lossy warp compression utility" << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "  shrink_warp dim input_warp.nii.gz outputwarp.nii.gz " << std::endl;
  return -1;
}

template <class TInput, class TOutput, unsigned int VDim>
class CompressFunctor 
{
public:
  typedef itk::Vector<double, VDim> SpacingType;
  typedef typename TInput::ComponentType InputComponentType;
  typedef typename TOutput::ComponentType OutputComponentType;
  typedef CompressFunctor<TInput, TOutput, VDim> Self;

  CompressFunctor() {}

  CompressFunctor(const SpacingType &spacing, double scale)
    {
    for(unsigned int d = 0; d < VDim; d++)
      m_MultFactor = scale * (1.0 / spacing[d]);
    }

  CompressFunctor(const Self &src)
    { m_MultFactor = src.m_MultFactor; }

  TOutput operator() (const TInput &in)
    {
    TOutput output(VDim);
    for(unsigned int d = 0; d < VDim; d++)
      {
      double x = in[d] * m_MultFactor[d];
      double rnd = std::floor(x + 0.5);
      output[d] = static_cast<OutputComponentType>(rnd / m_MultFactor[d]);
      }
    return output;
    }

  bool operator == (const Self &src)
    {
    return src.m_MultFactor == m_MultFactor;
    }

  bool operator != (const Self &src)
    {
    return src.m_MultFactor != m_MultFactor;
    }


private:
  SpacingType m_MultFactor;
};


template <int VDim> int mymain(int argc, char *argv[])
{
  if(argc != 4)
    return usage();

  char *fnInput = argv[argc-2];
  char *fnOutput = argv[argc-1];

  typedef itk::VectorImage<float, VDim> ImageType;
  typedef itk::VectorImage<float, VDim> CompressedImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fnInput);
  reader->Update();

  typedef CompressFunctor<
    typename ImageType::PixelType, typename CompressedImageType::PixelType, VDim> FunctorType;

  FunctorType functor(reader->GetOutput()->GetSpacing(), 10.0);
  typedef itk::UnaryFunctorImageFilter<ImageType, CompressedImageType, FunctorType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetFunctor(functor);
  filter->SetInput(reader->GetOutput());

  typedef itk::ImageFileWriter<CompressedImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(fnOutput);
  writer->SetInput(filter->GetOutput());
  writer->Update();

  return 0;
}




int main(int argc, char *argv[])
{
  if(argc != 4)
    return usage();

  int dim = atoi(argv[1]);

  if(dim == 2)
    return mymain<2>(argc, argv);
  else if(dim == 3)
    return mymain<3>(argc, argv);
  else
    {
    std::cerr << "The first parameter must be the dimensionality (2 or 3)" << std::endl;
    return usage();
    }
}
