#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryWellComposed3DImageFilter.h"

int main( int argc, char * argv[])
{
 
  const unsigned int      ImageDimension = 3;
  typedef int             PixelType;

  typedef itk::Image<PixelType,  ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[1]);
  reader->Update();  


/*  typedef itk::NeighborhoodIterator<ImageType>  NeighborhoodIteratorType;
  NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(1);
  
  NeighborhoodIteratorType It(radius, reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    if (!It.InBounds())
    {
      continue;
    } 
    It.SetCenterPixel(2);
  } 
  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholdFilterType;
  ThresholdFilterType::Pointer thresh = ThresholdFilterType::New();

  thresh->SetInput(reader->GetOutput());
  thresh->SetInsideValue(0);
  thresh->SetOutsideValue(1);
  thresh->SetLowerThreshold(0);
  thresh->SetUpperThreshold(120);
  thresh->Update();
*/
  typedef itk::BinaryWellComposed3DImageFilter<ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();

  filter->SetInput(reader->GetOutput());
  filter->SetBackgroundValue(static_cast<PixelType>(0));
  filter->SetForegroundValue(static_cast<PixelType>(1));
  filter->Update();
  filter->Print(std::cout);

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(filter->GetOutput());
  writer->SetFileName(argv[2]);
  writer->Update();


  return 0;
};















