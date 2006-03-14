#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkVoxBoCUBImageIOFactory.h"
#include "itkMetaDataObject.h"
#include "itkIOCommon.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
  // Enable support for VoxBo
  itk::ObjectFactoryBase::RegisterFactory(itk::VoxBoCUBImageIOFactory::New());
  
  if(argc < 2)
    {
    cerr << "Usage: identimg image_file" << endl;
    return -1;
    }

  typedef itk::Image<float, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(argv[1]);
  fltReader->Update();

  ImageType::Pointer image = fltReader->GetOutput();

  cout << "Image: " << argv[1] << endl;
  cout << "Dimensions: " << image->GetBufferedRegion() << endl;
  cout << "Spacing: " << image->GetSpacing() << endl;
  cout << "Origin: " << image->GetOrigin() << endl;

  // Compute the image histogram
  double iMin, iMax, iSum = 0;
  bool flagFirst = true;
  itk::ImageRegionIterator<ImageType> it(image, image->GetBufferedRegion());
  for( ; !it.IsAtEnd(); ++it)
    {
    double val = it.Value();
    if(!vnl_math_isfinite(val))
      continue;
    
    if(flagFirst)
      { 
      iMin = iMax = iSum = val;
      flagFirst = false;
      }
    else
      {
      iMax = val > iMax ? val : iMax;
      iMin = val < iMin ? val : iMin;
      iSum += val;
      }
    }

  cout << "Intensity Range: [" << iMin << ", " << iMax << "]" << endl;
  cout << "Intensity Mean: " << 
    iSum / image->GetBufferedRegion().GetNumberOfPixels() << endl;

  // See if there is an SPM origin
  string temp;
  if(itk::ExposeMetaData<std::string>(
      image->GetMetaDataDictionary(), itk::ITK_FileOriginator, temp))
    {
    short ox = (temp[1] << 8) + temp[0];
    short oy = (temp[3] << 8) + temp[2];
    short oz = (temp[5] << 8) + temp[4];
    cout << "SPM Origin: " << ox << "," << oy << "," << oz << endl;
    }

  // Print metadata
  cout << "Meta Data Dictionary" << endl;
  itk::MetaDataDictionary &mdd = image->GetMetaDataDictionary();
  itk::MetaDataDictionary::ConstIterator itMeta = mdd.Begin();
  while(itMeta != mdd.End())
    {
    // Get the metadata as a generic object
    std::string key = itMeta->first;
    itk::MetaDataObjectBase *meta = itMeta->second;

    // Check if the meta data string is a
    if( typeid(std::string) == meta->GetMetaDataObjectTypeInfo() )
      {
      // Cast the value to a string and print it
      typedef itk::MetaDataObject<std::string> ObjectType;
      std::string value = ((ObjectType *)(meta))->GetMetaDataObjectValue();

      // For some weird reason, some of the strings returned by this method
      // contain '\0' characters. We will replace them by spaces
      std::ostringstream sout("");
      for(unsigned int i=0;i<value.length();i++)
        if(value[i] >= ' ') sout << value[i];
      value = sout.str();

      // Make sure the value has more than blanks
      if(value.find_first_not_of(" ") != value.npos)
        {
        cout << key << " : " << value << endl;
        }
      }
    ++itMeta;
    }



  return 0;
}
