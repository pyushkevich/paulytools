#include <iostream>
#include <string>
#include <sstream>
#include <set>

#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <itkImageFileReader.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkImage.h>
#include "ReadWriteVTK.h"

using namespace std;
using namespace itk;

int usage()
{
  cout << "WarpMeshBackward - Applies a warp field (Brian's format) to a VTK mesh" << endl;
  cout << "usage: " << endl;
  cout << "   WarpMeshBackward mesh.vtk warpname out.vtk " << endl;
  return -1;
}

int main(int argc, char **argv)
{
  // Check the parameters
  if(argc < 4) return usage();

  // Read the mesh
  vtkPolyData *p = ReadVTKData(argv[1]);

  // Read the warp field
  typedef itk::Image<double,3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::LinearInterpolateImageFunction<ImageType,double> FuncType;
  FuncType::Pointer warp[3];
  for(size_t i = 0; i < 3; i++)
    {
    ReaderType::Pointer fltReader = ReaderType::New();
    std::string xyz = "xyz";
    std::ostringstream oss;
    oss << argv[2] << xyz[i] << "vec.img";
    cout << "Reading " << oss.str() << endl;
    fltReader->SetFileName(oss.str().c_str());
    fltReader->Update();
    warp[i] = FuncType::New();
    warp[i]->SetInputImage(fltReader->GetOutput());
    }

  // Update the coordinates
  for(size_t k = 0; k < p->GetNumberOfPoints(); k++)
    {
    double *pt = p->GetPoint(k);
    FuncType::ContinuousIndexType idx;
    idx[0] = pt[0]; idx[1] = pt[1]; idx[2] = pt[2];
    double newpt[3];
    double x = pt[0] + warp[0]->EvaluateAtContinuousIndex(idx);
    double y = pt[1] + warp[1]->EvaluateAtContinuousIndex(idx);
    double z = pt[2] + warp[2]->EvaluateAtContinuousIndex(idx);
    p->GetPoints()->SetPoint(k, x, y, z);
    }
    
  // Write the mesh
  WriteVTKData(p, argv[3]);
  return 0;
}
