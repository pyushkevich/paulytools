#include "ConvertImage3D.h"

int main(int argc, char *argv[])
{
  // Load the ITK factories
  itk::ObjectFactoryBase::RegisterFactory(itk::VoxBoCUBImageIOFactory::New());
  itk::ObjectFactoryBase::RegisterFactory(itk::PovRayDF3ImageIOFactory::New());

  ImageConverter<double, 3> convert;
  return convert.ProcessCommandLine(argc, argv);
}

template class ImageConverter<double, 3>;
