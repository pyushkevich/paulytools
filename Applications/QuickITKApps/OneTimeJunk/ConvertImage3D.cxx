#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkVoxBoCUBImageIOFactory.h"
#include <itkImageRegionIterator.h>
#include <itkExtractImageFilter.h>
#include <itkByteSwapper.h>

#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

class ConverterBase
{
public:
  virtual void Convert(const char *in, const char *out, double *origin, double xScale, double xShift) = 0;
};

template<class TPixelIn, class TPixelOut, unsigned int Dim>
class Converter : public ConverterBase
{
public:
  typedef itk::Image<TPixelIn, Dim> InputImageType;
  typedef itk::Image<TPixelOut, Dim> OutputImageType;
  typedef itk::ImageFileReader<InputImageType> ReaderType;
  typedef itk::ImageFileWriter<OutputImageType> WriterType;

  void ScaleAndShift(InputImageType *in, OutputImageType *out, double xScale, double xShift)
    {
    typedef itk::ImageRegionIterator<InputImageType> IteratorIn;
    typedef itk::ImageRegionIterator<OutputImageType> IteratorOut;
    IteratorIn itIn(in, in->GetBufferedRegion());
    IteratorOut itOut(out, out->GetBufferedRegion());
    for( ; !itIn.IsAtEnd(); ++itIn, ++itOut)
      {
      TPixelIn z = itIn.Get();
      if(!vnl_math_isfinite(z))
        z = 0;
      
      if(xScale == 1.0 && xShift == 0.0)
        itOut.Set((TPixelOut) z);
      else
        {
        double val = xShift + xScale * z;
        if(typeid(TPixelOut) != typeid(float))
          itOut.Set((TPixelOut)round(val));
        else
          itOut.Set((TPixelOut)val);
        }
      }
    }
  
  void Convert(const char *in, const char *out, double *origin, 
    double xScale = 1.0, double xShift = 0.0)
    {
    // Read the input image
    typename ReaderType::Pointer fltReader = ReaderType::New();
    fltReader->SetFileName(in);
    try {
      fltReader->Update();
    } 
    catch (itk::ExceptionObject &exc) {
      cerr << "Exception loading image!" << endl;
      cerr << exc << endl;
      return;
    }

    // Input and output images
    typename InputImageType::Pointer imgInput = fltReader->GetOutput();
    typename OutputImageType::Pointer imgOutput = OutputImageType::New();
    imgOutput->SetRegions(imgInput->GetBufferedRegion());
    imgOutput->SetSpacing(imgInput->GetSpacing());
    imgOutput->Allocate();

    // Set the origin if necessary
    if(origin)
      imgOutput->SetOrigin(origin);
    else
      imgOutput->SetOrigin(imgInput->GetOrigin());

    // Scale and shift image
    ScaleAndShift(imgInput, imgOutput, xScale, xShift);

    // Check for custom writer
    if(strstr(out,".df3") != NULL)
      {
      WriteDF3(imgOutput, out);
      }
    else
      {
      // Save the output image
      typename WriterType::Pointer fltWriter = WriterType::New();
      fltWriter->SetInput(imgOutput);
      fltWriter->SetFileName(out);
      fltWriter->Update();
      }
    }

  void WriteDF3(OutputImageType *img, const char *file)
    {
    FILE *f = fopen(file, "wb");
    unsigned short xyz[3];
    xyz[0] = (unsigned short) img->GetBufferedRegion().GetSize(0);
    xyz[1] = (unsigned short) img->GetBufferedRegion().GetSize(1);
    xyz[2] = (unsigned short) img->GetBufferedRegion().GetSize(2);

    itk::ByteSwapper<unsigned short>::SwapRangeFromSystemToBigEndian(xyz, 3);
    itk::ByteSwapper<TPixelOut>::SwapRangeFromSystemToBigEndian(
      img->GetBufferPointer(), 
      img->GetBufferedRegion().GetNumberOfPixels());

    cout << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;
    cout << (int)(((char *)xyz)[0]) << " ";
    cout << (int)(((char *)xyz)[1]) << " ";
    cout << (int)(((char *)xyz)[2]) << " ";
    cout << (int)(((char *)xyz)[3]) << " ";
    cout << (int)(((char *)xyz)[4]) << " ";
    cout << (int)(((char *)xyz)[5]) << endl;

    fwrite(xyz, sizeof(char), 3 * sizeof(short), f);
    fwrite(img->GetBufferPointer(), sizeof(TPixelOut), img->GetBufferedRegion().GetNumberOfPixels(), f);
    fclose(f);
    }
};

int usage()
{
  cout << "USAGE: [flags] program input_image output_image" << endl;
  cout << "FLAGS: " << endl;
  cout << "  -f                 Convert to floating point" << endl;
  cout << "  -b                 Convert to byte" << endl;
  cout << "  -d                 Convert to double" << endl;
  cout << "  -us                Convert to unsigned short" << endl;
  cout << "  -o X Y Z           Set the origin in the output image" << endl;
  cout << "  -scale X.XX        Scale the image by factor x.xx" << endl;
  cout << "  -shift X.XX        Shift the image intensity by x.xx" << endl;
  cout << "  -2                 Perform operation on 2D images" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  // Enable support for VoxBo
  itk::ObjectFactoryBase::RegisterFactory(itk::VoxBoCUBImageIOFactory::New());
  
  if(argc < 3)
    return usage();

  // Create abstract converter
  ConverterBase *converter = NULL;

  // Origin vector
  double *origin = NULL;
  double xScale = 1.0, xShift = 0.0;
  int iSlice = -1;

  unsigned int k = 3;
  for(int j = 1; j < argc; j++)
    if(strcmp(argv[j],"-2") == 0)
      k = 2;

  cout << "DOING " << k << "-dimensional" << endl;

  // Parse arguments
  for(int i = 1; i < argc; i++)
    {
    std::string arg(argv[i]);
    if(arg == "-f")
      if(k==2) converter = new Converter<float, float, 2>; 
      else converter = new Converter<float, float, 3>;
    else if(arg == "-b")
      if(k==2) converter = new Converter<unsigned char, unsigned char, 2>; 
      else converter = new Converter<unsigned char, unsigned char, 3>;
    else if(arg == "-d")
      if(k==2) converter = new Converter<double, double, 2>; 
      else converter = new Converter<double, double, 3>;
    else if(arg == "-us")
      if(k==2) converter = new Converter<unsigned short, unsigned short, 2>; 
      else converter = new Converter<unsigned short, unsigned short, 3>;
    else if(arg == "-fus")
      if(k==2) converter = new Converter<float, unsigned short, 2>; 
      else converter = new Converter<float, unsigned short, 3>;
    else if(arg == "-scale")
      xScale = atof(argv[++i]);
    else if(arg == "-shift")
      xShift = atof(argv[++i]);
    else if(arg == "-o")
      {
      origin = new double[3];
      origin[0] = atof(argv[++i]);
      origin[1] = atof(argv[++i]);
      origin[2] = atof(argv[++i]);
      }
    }

  // Check the converter
  if(!converter)
    if(k==2) converter =  new Converter<short, short, 2>; else converter = new Converter<short, short, 3>;

  // Perform the conversion
  converter->Convert(argv[argc-2], argv[argc-1], origin, xScale, xShift);

  // Clean up
  if(origin)
    delete origin;
  delete converter;
}
