#ifndef __ReadWriteImage_h_
#define __ReadWriteImage_h_

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

template <class TImageType> 
void ReadImage(itk::SmartPointer<TImageType> &target, const char *file)
{
  typename itk::ImageFileReader<TImageType>::Pointer reader = 
    itk::ImageFileReader<TImageType>::New();
  reader->SetFileName(file);
  reader->Update();
  target = reader->GetOutput();
}

template <class TImageType> 
void WriteImage(itk::SmartPointer<TImageType> image, const char *file)
{
  typename itk::ImageFileWriter<TImageType>::Pointer writer = 
    itk::ImageFileWriter<TImageType>::New();
  writer->SetFileName(file);
  writer->SetInput(image);
  writer->Update();
}

#endif
