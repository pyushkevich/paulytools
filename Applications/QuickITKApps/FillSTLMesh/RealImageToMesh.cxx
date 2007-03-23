#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkVTKImageExport.h"
#include "vtkImageImport.h"
#include "vtkImageData.h"
#include "vtkImageMarchingCubes.h"
#include "ReadWriteVTK.h"
#include <iostream>

using namespace std;

int usage()
{
  cout << "Usage: real2mesh input.img output.vtk threshold" << endl;
  return -1;
}

template<class TImage>
void ConnectITKToVTK(itk::VTKImageExport<TImage> *fltExport,vtkImageImport *fltImport)
{
  fltImport->SetUpdateInformationCallback( fltExport->GetUpdateInformationCallback());
  fltImport->SetPipelineModifiedCallback( fltExport->GetPipelineModifiedCallback());
  fltImport->SetWholeExtentCallback( fltExport->GetWholeExtentCallback());
  fltImport->SetSpacingCallback( fltExport->GetSpacingCallback());
  fltImport->SetOriginCallback( fltExport->GetOriginCallback());
  fltImport->SetScalarTypeCallback( fltExport->GetScalarTypeCallback());
  fltImport->SetNumberOfComponentsCallback( fltExport->GetNumberOfComponentsCallback());
  fltImport->SetPropagateUpdateExtentCallback( fltExport->GetPropagateUpdateExtentCallback());
  fltImport->SetUpdateDataCallback( fltExport->GetUpdateDataCallback());
  fltImport->SetDataExtentCallback( fltExport->GetDataExtentCallback());
  fltImport->SetBufferPointerCallback( fltExport->GetBufferPointerCallback());
  fltImport->SetCallbackUserData( fltExport->GetCallbackUserData());
}


int main(int argc, char *argv[])
{
  if(argc != 4)
    return usage();

  // Read the input image
  typedef itk::Image<float, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(argv[1]);
  fltReader->Update();
  ImageType::Pointer imgInput = fltReader->GetOutput();

  // Get the range of the input image
  float imax = imgInput->GetBufferPointer()[0];
  float imin = imax;
  for(size_t i = 0; i < imgInput->GetBufferedRegion().GetNumberOfPixels(); i++)
    {
    float x = imgInput->GetBufferPointer()[i];
    imax = std::max(imax, x);
    imin = std::min(imin, x);
    }

  float cut = imin + (imax - imin) * atof(argv[3]);
  cout << "Image Range: [" << imin << ", " << imax << "]" << endl;
  cout << "Taking level set at " << cut << endl;

  // Create an importer and an exporter in VTK
  typedef itk::VTKImageExport<ImageType> ExporterType;
  ExporterType::Pointer fltExport = ExporterType::New();
  fltExport->SetInput(imgInput);
  vtkImageImport *fltImport = vtkImageImport::New();
  ConnectITKToVTK(fltExport.GetPointer(), fltImport);

  // Run marching cubes on the input image
  vtkImageMarchingCubes *fltMarching = vtkImageMarchingCubes::New();
  fltMarching->SetInput(fltImport->GetOutput());
  fltMarching->ComputeScalarsOff();
  fltMarching->ComputeGradientsOff();
  fltMarching->ComputeNormalsOn();
  fltMarching->SetNumberOfContours(1);
  fltMarching->SetValue(0,cut);
  fltMarching->Update();
  vtkPolyData *meshCubes = fltMarching->GetOutput();

  // Write the output
  WriteVTKData(meshCubes, argv[2]);
}
