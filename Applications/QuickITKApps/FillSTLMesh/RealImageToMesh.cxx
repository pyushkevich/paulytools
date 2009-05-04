#include "itkImageFileReader.h"
#include "itkOrientedRASImage.h"
#include "itkVTKImageExport.h"
#include "vtkImageImport.h"
#include "vtkImageData.h"
#include "vtkImageMarchingCubes.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"
#include <vtkMatrix4x4.h>
#include "vnl/vnl_matrix_fixed.h"
#include "ReadWriteVTK.h"
#include <iostream>

using namespace std;

int usage()
{
  cout << "Usage: vtklevelset input.img output.vtk threshold" << endl;
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
  typedef itk::OrientedRASImage<float, 3> ImageType;
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

  float cut = atof(argv[3]);
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

  // Map to the RAS coordinates
  for(size_t i = 0; i < meshCubes->GetNumberOfPoints(); i++)
    {
    double *pvtk = meshCubes->GetPoint(i);
    itk::ContinuousIndex<double, 3> idx;
    idx[0] = (pvtk[0] - imgInput->GetOrigin()[0]) / imgInput->GetSpacing()[0];
    idx[1] = (pvtk[1] - imgInput->GetOrigin()[1]) / imgInput->GetSpacing()[1];
    idx[2] = (pvtk[2] - imgInput->GetOrigin()[2]) / imgInput->GetSpacing()[2];
    itk::Point<double, 3> pitk;
    imgInput->TransformContinuousIndexToRASPhysicalPoint(idx, pitk);
    pvtk[0] = pitk[0];
    pvtk[1] = pitk[1];
    pvtk[2] = pitk[2];
    meshCubes->GetPoints()->SetPoint(i, pitk[0], pitk[1], pitk[2]);
    }

  // Write the output
  WriteVTKData(meshCubes, argv[2]);
}
