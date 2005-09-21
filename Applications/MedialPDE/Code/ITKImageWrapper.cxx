#include "ITKImageWrapper.h"

#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

template<typename TPixel>
class ITKImageWrapperImpl : public virtual ITKImageWrapper<TPixel>
{
public:
  // The typedefs for this class
  typedef ITKImageWrapper<TPixel> Superclass;
  typedef typename Superclass::ImageType ImageType;
  typedef itk::LinearInterpolateImageFunction<ImageType, float> InterpolatorType;
  
  // Load the image from a file
  void LoadFromFile(const char *file);

  // Save the image to a file
  void SaveToFile(const char *file);

  // Set the internal image to some ITK image
  void SetInternalImage(ImageType *xImage) 
    {
    if(xImage == NULL)
      this->xImage = ImageType::New();
    else
      this->xImage = xImage;
    OnImageUpdate();
    }

  // Check whether an image is loaded into the wrapper
  virtual bool IsImageLoaded()
    { return xImage->GetBufferedRegion().GetSize()[0] > 0; }

  // Get the basic image properties
  unsigned int GetImageSize(unsigned int d)
    { return (unsigned int) xImage->GetBufferedRegion().GetSize()[d]; }
  
  double GetImageSpacing(unsigned int d)
    { return xImage->GetSpacing()[d]; }
  
  double GetImageOrigin(unsigned int d)
    { return xImage->GetOrigin()[d]; }

  // Interpolate the image at a continuous index, or return background if out of bounds
  float Interpolate(float x, float y, float z, float xBackground);

  // Get the internal ITK image pointer
  ImageType *GetInternalImage()
    { return xImage.GetPointer(); }

private:
  // The constructor
  ITKImageWrapperImpl();
  
  // The internaly stored image
  typename ImageType::Pointer xImage;

  // The linear interpolator for the currently stored image
  typename InterpolatorType::Pointer fnInterpolator;

  // Function called when the image is updated externally
  void OnImageUpdate();

  // Factory is our friend
  friend class ITKImageWrapperFactory<TPixel>;
};

template<typename TPixel>
ITKImageWrapperImpl<TPixel>
::ITKImageWrapperImpl()
{
  xImage = ImageType::New();
  fnInterpolator = InterpolatorType::New();
}


template<typename TPixel>
float
ITKImageWrapperImpl<TPixel>
::Interpolate(float x, float y, float z, float xBackground)
{
  // Create a point with the passed in coordinates
  itk::Point<float, 3> xPoint;
  xPoint[0] = x; xPoint[1] = y; xPoint[2] = z;

  // Create a continuous index from the point
  itk::ContinuousIndex<float, 3> xIndex;
  xImage->TransformPhysicalPointToContinuousIndex(xPoint, xIndex);

  // Interpolate at the index
  if(fnInterpolator->IsInsideBuffer(xIndex))
    return (float) fnInterpolator->EvaluateAtContinuousIndex(xIndex);
  else
    return xBackground;
}

template<typename TPixel>
void ITKImageWrapperImpl<TPixel>
::OnImageUpdate()
{
  fnInterpolator->SetInputImage(xImage);
}

template<typename TPixel>
void ITKImageWrapperImpl<TPixel>
::LoadFromFile(const char *file)
{
  // Load the image
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(file);
  fltReader->Update();
  xImage = fltReader->GetOutput();

  // Refresh other components
  OnImageUpdate();
}

template<typename TPixel>
void ITKImageWrapperImpl<TPixel>
::SaveToFile(const char *file)
{
  // Load the image
  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer fltWriter = WriterType::New();
  fltWriter->SetInput(xImage);
  fltWriter->SetFileName(file);
  fltWriter->Update();
}

template <typename TPixel>
ITKImageWrapper<TPixel> *
ITKImageWrapperFactory<TPixel>
::NewImageWrapper()
{
  // Create a new instance of image wrapper impl
  return new ITKImageWrapperImpl<TPixel>;
}

// Define image wrappers for necessary types
template class ITKImageWrapperFactory<float>;
template class ITKImageWrapperFactory<unsigned char>;

