/********************************************************************************
 * Create a VTK mesh from the white matter image
 *******************************************************************************/
#ifndef __BinaryImageToMeshFilter_h_
#define __BinaryImageToMeshFilter_h_

#include "itkImage.h"
#include "itkAntiAliasBinaryImageFilter.h"
#include "itkCommand.h"

#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkImageMarchingCubes.h>
#include <vtkImageImport.h>
#include <itkVTKImageExport.h>
#include <vtkDecimate.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkStripper.h>

#include <iostream>

using namespace std;

                                                                                                                                                                                   
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

template<class TImage>
class BinaryImageToMeshFilter : public itk::ProcessObject
{
public:
  typedef BinaryImageToMeshFilter Self;
  typedef itk::ProcessObject Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;
  typedef itk::Image<float,3> FloatImageType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int, TImage::ImageDimension);

  /** Get the result mesh */
  vtkPolyData *GetMesh()
    {
    return fltStripper->GetOutput();
    }

  /** Whether to invert the binary image */
  itkSetMacro(InvertInput,bool);
  itkGetMacro(InvertInput,bool);

  /** Set the input */
  void SetInput(TImage *image)
    {
    this->SetNthInput(0,image);
    }

  /** Update method (why?) */
  void Update()
    {
    this->GenerateData();
    }

  /** Set the anti-aliasing quality parameter */
  void SetAntiAliasMaxRMSError(double value) 
    {
    fltAlias->SetMaximumRMSError(value);
    }

  /** Get the 'distance image' based on anti-aliasing the binary image */
  FloatImageType *GetDistanceImage()
    {
    return fltAlias->GetOutput();
    }

protected:
  BinaryImageToMeshFilter() 
    {
    this->SetNumberOfInputs(1);
    this->SetNumberOfOutputs(1);
  
    // Create an anti-aliasing image filter
    fltAlias = AAFilter::New();
    fltAlias->SetMaximumRMSError(0.024);

    // Cast the image to VTK
    fltExport = ExportFilter::New();
    fltImport = vtkImageImport::New();
    fltExport->SetInput(fltAlias->GetOutput());
    ConnectITKToVTK(fltExport.GetPointer(),fltImport);

    // Compute marching cubes
    fltMarching = vtkImageMarchingCubes::New();
    fltMarching->SetInput(fltImport->GetOutput());
    fltMarching->ComputeScalarsOff();
    fltMarching->ComputeGradientsOff();
    fltMarching->SetNumberOfContours(1);
    fltMarching->SetValue(0,0.0f);

    // Keep the largest connected component
    fltConnect = vtkPolyDataConnectivityFilter::New();
    fltConnect->SetInput(fltMarching->GetOutput());
    fltConnect->SetExtractionModeToLargestRegion();

    // Compute triangle strips for faster display
    fltStripper = vtkStripper::New();
    fltStripper->SetInput(fltConnect->GetOutput());

    // Set up progress
    typedef itk::SimpleMemberCommand<Self> CommandType;
    typename CommandType::Pointer cmd = CommandType::New();
    cmd->SetCallbackFunction(this, &Self::ProgressCommand);
    fltAlias->AddObserver(itk::ProgressEvent(),cmd);

    // Invert - NO
    m_InvertInput = false;
    }

  ~BinaryImageToMeshFilter()
    {
    // CLean up
    fltMarching->Delete();
    fltConnect->Delete();
    fltImport->Delete();
    fltStripper->Delete();
    }

  /** Generate Data */
  virtual void GenerateData( void )
    {
    // Get the input and output pointers
    typename TImage::ConstPointer inputImage = 
      reinterpret_cast<TImage *>(this->GetInput(0));
    fltAlias->SetInput(inputImage);

    // Run the computation
    cout << "Computing white matter mesh" << endl;

    // Make a negative image
    if(m_InvertInput)
      {
      typename TImage::Pointer imgInverse = TImage::New();
      imgInverse->SetRegions(inputImage->GetBufferedRegion());
      imgInverse->Allocate();

      typename TImage::PixelType iMax = 
        std::numeric_limits< typename TImage::PixelType >::max();
      
      ConstIteratorType itSrc(inputImage, inputImage->GetBufferedRegion());
      IteratorType itTrg(imgInverse, imgInverse->GetBufferedRegion());
      while(!itSrc.IsAtEnd())
        {
        itTrg.Set(iMax - itSrc.Get());
        ++itTrg; ++itSrc;
        }

      inputImage = imgInverse;
      }

    // Run the filters
    cout << "   anti-aliasing the image " << endl;
    fltAlias->Update();

    cout << "   converting image to VTK" << endl;
    fltImport->Update();

    cout << "   running marching cubes algorithm" << endl;
    fltMarching->Update();

    cout << "      mesh has " << fltMarching->GetOutput()->GetNumberOfCells() << " cells." << endl;

    cout << "   extracting the largest component" << endl;
    fltConnect->Update();

    cout << "      mesh has " << fltConnect->GetOutput()->GetNumberOfCells() << " cells." << endl;

    cout << "   converting mesh to triangle strips" << endl;
    fltStripper->Update();
    }

private:
  typedef itk::ImageRegionIterator<TImage> IteratorType;
  typedef itk::ImageRegionConstIterator<TImage> ConstIteratorType;
  typedef itk::AntiAliasBinaryImageFilter<TImage,FloatImageType> AAFilter;
  typedef itk::VTKImageExport<FloatImageType> ExportFilter;

  typename AAFilter::Pointer fltAlias;
  typename ExportFilter::Pointer fltExport;
  vtkImageImport *fltImport;
  vtkImageMarchingCubes *fltMarching;
  vtkPolyDataConnectivityFilter *fltConnect;
  vtkStripper *fltStripper;

  bool m_InvertInput;

  void ProgressCommand() 
    {
    cout << "." << flush;
    }
};

#endif

