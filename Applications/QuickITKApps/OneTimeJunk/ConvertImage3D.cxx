#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkVoxBoCUBImageIOFactory.h"
#include <itkImageRegionIterator.h>
#include <itkExtractImageFilter.h>

#include <iostream>

using namespace std;
using namespace itk;

class ConverterBase
{
public:
  virtual void Convert(const char *in, const char *out, double *origin, double xScale, double xShift) = 0;
};

template<class TPixel, unsigned int Dim>
class Converter : public ConverterBase
{
public:
  typedef Image<TPixel, Dim> ImageType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef ImageFileReader<ImageType> ReaderType;

  void ScaleAndShift(ImageType *img, double xScale, double xShift)
    {
    typedef ImageRegionIterator<ImageType> Iterator;
    Iterator it(img, img->GetBufferedRegion());
    for( ; !it.IsAtEnd(); ++it)
      {
      double val = xShift + xScale * it.Get();
      if(typeid(TPixel) != typeid(float))
        it.Set((TPixel)round(val));
      else
        it.Set((TPixel)val);
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
    typename ImageType::Pointer img = fltReader->GetOutput();

    // Set the origin if necessary
    if(origin)
      img->SetOrigin(origin);

    // Scale and shift image
    if(xScale != 1.0 || xShift != 0.0)
      ScaleAndShift(img, xScale, xShift);

    // Save the output image
    typedef ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer fltWriter = WriterType::New();
    fltWriter->SetInput(fltReader->GetOutput());
    fltWriter->SetFileName(out);
    fltWriter->Update();
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
    string arg(argv[i]);
    if(arg == "-f")
      if(k==2) converter = new Converter<float,2>; else converter = new Converter<float,3>;
    else if(arg == "-b")
      if(k==2) converter = new Converter<unsigned char,2>; else converter = new Converter<unsigned char,3>;
    else if(arg == "-d")
      if(k==2) converter =  new Converter<double,2>; else converter = new Converter<double,3>;
    else if(arg == "-us")
      if(k==2) converter =  new Converter<unsigned short,2>; else converter = new Converter<unsigned short,3>;
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
    if(k==2) converter =  new Converter<short,2>; else converter = new Converter<short,3>;

  // Perform the conversion
  converter->Convert(argv[argc-2], argv[argc-1], origin, xScale, xShift);

  // Clean up
  if(origin)
    delete origin;
  delete converter;
}
