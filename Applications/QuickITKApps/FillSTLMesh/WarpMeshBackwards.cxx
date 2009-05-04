#include <iostream>
#include <string>
#include <sstream>
#include <set>

#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <itkImageFileReader.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkOrientedRASImage.h>
#include <itkImage.h>
#include "ReadWriteVTK.h"

using namespace std;
using namespace itk;

int usage()
{
  cout << "WarpMeshBackward - Applies a warp field (Brian's format) to a VTK mesh" << endl;
  cout << "usage: " << endl;
  cout << "   WarpMeshBackward [options] mesh.vtk warpname out.vtk " << endl;
  cout << "options: " << endl;
  cout << "   -e ext : Specify the extension of the warp files (default: img)" << endl;
  return -1;
}

int main(int argc, char **argv)
{
  // Check the parameters
  if(argc < 4) return usage();

  // Check the extension
  string ext = "img";
  for(size_t iarg = 1; iarg < argc - 3; iarg++)
    {
    if(0 == strcmp(argv[iarg], "-e"))
      {
      ext = string(argv[++iarg]);
      }
    else return usage();
    }

  // Read the mesh
  vtkPolyData *p = ReadVTKData(argv[1]);

  // Read the warp field
  typedef itk::OrientedRASImage<double,3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::LinearInterpolateImageFunction<ImageType,double> FuncType;
  FuncType::Pointer warp[3];
  ImageType::Pointer wimg[3];
  for(size_t i = 0; i < 3; i++)
    {
    ReaderType::Pointer fltReader = ReaderType::New();
    std::string xyz = "xyz";
    std::ostringstream oss;
    oss << argv[2] << xyz[i] << "vec.nii.gz";
    cout << "Reading " << oss.str() << endl;
    fltReader->SetFileName(oss.str().c_str());
    fltReader->Update();
    warp[i] = FuncType::New();
    wimg[i] = fltReader->GetOutput();
    warp[i]->SetInputImage(wimg[i]);
    }

  // Update the coordinates
  for(size_t k = 0; k < p->GetNumberOfPoints(); k++)
    {
    // Get the point (in RAS coords)
    double *pt = p->GetPoint(k);
    ImageType::PointType itk_point;
    itk_point[0] = pt[0]; itk_point[1] = pt[1]; itk_point[2] = pt[2];
    
    // Map the point to a continuous index
    FuncType::ContinuousIndexType idx;
    wimg[0]->TransformRASPhysicalPointToContinuousIndex(itk_point, idx); 
    cout << "Evaluate at index " << idx[0] << " " << idx[1] << " " << idx[2] << endl;

    // Compute the transformation. We assume the transformation is in ITK
    // physical space. 
    vnl_vector_fixed<double, 4> w;
    w[0] = warp[0]->EvaluateAtContinuousIndex(idx);
    w[1] = warp[1]->EvaluateAtContinuousIndex(idx);
    w[2] = warp[2]->EvaluateAtContinuousIndex(idx);
    w[3] = 0.0;

    // Map the transformation to RAS
    // vnl_vector_fixed<double, 4> w_ras = 
    //  wimg[0]->GetSpacingOriginPhysicalSpaceToRASPhysicalSpaceMatrix() * w;

    // vnl_vector_fixed<double, 4> w_ras = w;
      
    // Assume transformation is in spacing units
    itk::FixedArray<double, 3> aw, awras;
    aw[0] = w[0];
    aw[1] = w[1];
    aw[2] = w[2];
    wimg[0]->TransformLocalVectorToPhysicalVector(aw, awras);

    p->GetPoints()->SetPoint(k, pt[0] + awras[0], pt[1] + awras[1], pt[2] + awras[2]);
    // p->GetPoints()->SetPoint(k, pt[0] + w_ras[0], pt[1] + w_ras[1], pt[2] + w_ras[2]);
    }
    
  // Write the mesh
  WriteVTKData(p, argv[3]);
  return 0;
}
