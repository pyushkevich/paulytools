#include <iostream>
#include <string>
#include <sstream>

#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkLODActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkBYUReader.h>
#include <vtkSTLReader.h>
#include <vtkTriangleFilter.h>

#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>
#include <itkSignedDanielssonDistanceMapImageFilter.h>
#include <itkImage.h>

#include "DrawTriangles.h"

using namespace std;
using namespace itk;

int usage()
{
  cout << "meshdiff - compare two VTK meshes" << endl;
  cout << "usage: " << endl;
  cout << "   meshdiff [options] mesh1.vtk mesh2.vtk" << endl;
  cout << "options: " << endl;
  cout << "   -s X.XX             Size of the voxel for overlap/distance measurements" << endl;
  return -1;
}

typedef itk::Image<unsigned char, 3> ByteImageType;
typedef itk::Image<float, 3> FloatImageType;

vtkPolyData *ReadVTKData(string fn)
{
  vtkPolyData *p1 = NULL;

  // Choose the reader based on extension

  if(fn.rfind(".byu") == fn.length() - 4)
    {
    vtkBYUReader *reader = vtkBYUReader::New();
    reader->SetFileName(fn.c_str());
    reader->Update();
    p1 = reader->GetOutput();
    }
  else if(fn.rfind(".vtk") == fn.length() - 4)
    {
    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    reader->SetFileName(fn.c_str());
    reader->Update();
    p1 = reader->GetOutput();
    }
  else
    {
    cout << "Could not find a reader for " << fn << endl;
    cout << fn.rfind(".byu") << endl;
    cout << fn.rfind(".vtk") << endl;
    return NULL;
    }

  // Convert the model to triangles
  vtkTriangleFilter *tri = vtkTriangleFilter::New();
  tri->SetInput(p1);
  cout << "Converting to triangles ..." << endl;
  tri->Update();
  vtkPolyData *pd = tri->GetOutput();

  return pd;
}

void ComputeDistanceTransform(ByteImageType *binary, FloatImageType::Pointer &dist)
{
  typedef SignedDanielssonDistanceMapImageFilter<
    ByteImageType, FloatImageType> DistanceMapType;
  DistanceMapType::Pointer fltDistance = DistanceMapType::New();
  fltDistance->SetInput(binary);
  fltDistance->SetSquaredDistance(float);
  fltDistance->SetInputIsBinary(true);
  fltDistance->SetUseImageSpacing(true);
  dist = fltDistance->Update();
}

void IntegrateDistance(vtkPolyData *pd, FloatImageType *dist, 
  double &xOneNorm, double &xTwoNorm, double &xInfNorm)
{
  // Set up an interpolator
  typedef itk::LinearInterpolateImageFunction<FloatImageType, double> FuncType;
  FuncType::Pointer fnInterp = FuncType::New();
  fnInterp->SetInputImage(dist);

  // Clear the norms
  xOneNorm = xTwoNorm = xInfNorm = 0.0;
  
  // Go over all the points in the image
  for(unsigned int i = 0; i < pd->GetNumberOfPoints(); i++)
    {
    // Get the point
    FuncType::PointType P(pd->GetPoint(i));
    double d = fabs(fnInterp->Evaluate(P));
    
    if(d > xInfNorm) xInfNorm = d;
    xOneNorm += d;
    xTwoNorm += d * d;
    }

  xOneNorm /= pd->GetNumberOfPoints();
  xTwoNorm = sqrt( xTwoNorm / pd->GetNumberOfPoints() );
}

void ScanConvertPolyData(vtkPolyData *pd, double *b0, double *b1, int *res, ByteImageType *img)
{
  // Create a vertex table from the polydata
  unsigned int nt = pd->GetNumberOfPolys(), it = 0;
  double **vtx = new double*[nt*3];  

  vtkCellArray *poly = pd->GetPolys();
  vtkIdType npts;
  vtkIdType *pts;
  for(poly->InitTraversal();poly->GetNextCell(npts,pts);)
    {
    for(unsigned int i=0;i<3;i++)
      {
      double *x = pd->GetPoints()->GetPoint(pts[i]);
      vtx[it] = (double *) malloc(3*sizeof(double));
      for(unsigned int j=0;j<3;j++)
        {
        vtx[it][j] = res[j] * (x[j] - b0[j]) / b1[j];
        }

      ++it;
      }      
    }

  // Create a ITK image to store the results
  ByteImageType::RegionType region;
  region.SetSize(0,res[0]); region.SetSize(1,res[1]); region.SetSize(2,res[2]);
  img->SetRegions(region);
  img->Allocate();
  img->FillBuffer(0);

  // Convert the polydata to an image
  drawBinaryTrianglesFilled(img->GetBufferPointer(), res, vtx, nt);

  // Set the origin and spacing of the image
  ByteImageType::SpacingType xSpacing;
  ByteImageType::PointType xOrigin;
  for(unsigned int d = 0; d < 3; d++)
    {
    xOrigin[d] = b0[d];
    xSpacing[d] = b1[d] / res[d];
    }
  img->SetOrigin(xOrigin);
  img->SetSpacing(xSpacing);
}


int main(int argc, char **argv)
{
  // Parameter values
  double bb[2][3] = {{0.,0.,0.},{128.,128.,128.}};
  double xVox = 1.0;
  int res[3];

  // Input specifications
  string fn1, fn2;

  // Check the parameters
  if(argc < 3) return usage();

  // Get the file names
  fn1 = argv[argc-2];
  fn2 = argv[argc-1];

  // Parse the optional parameters
  try 
    {
    for(unsigned int i=1;i<argc-2;i++)
      {
      string arg = argv[i];
      if(arg == "-s")
        {
        xVox = atof(argv[++i]);
        }
      }
    }
  catch(...) 
    {
    cout << "Error parsing command line options!" << endl;
    return usage();
    }

  // The real program begins here

  // Read the appropriate meshes
  vtkPolyData *p1 = ReadVTKData(fn1);
  vtkPolyData *p2 = ReadVTKData(fn2);

  // Get the extents of the data
  double *b1 = p1->GetBounds();  
  double *b2 = p2->GetBounds();  

  // Compute the bounding box
  for(size_t i = 0; i < 3; i++) 
    {
    bb[0][i] = (b1[2*i] < b2[2*i] ? b1[2*i] : b2[2*i]) - 4;
    bb[1][i] = (b1[2*i+1] < b2[2*i+1] ? b1[2*i+1] : b2[2*i+1]) + 4;
    cout << "Bounds[" << i << "]: " << bb[0][i] << " to " << bb[1][i] << endl;
    res[i] = (int)(ceil((bb[1][i] - bb[0][i]) / xVox));
    }

  // Scan convert the mesh to images
  ByteImageType::Pointer i1 = ByteImageType::New();
  ByteImageType::Pointer i2 = ByteImageType::New();
  ScanConvertPolyData(p1, bb[0], bb[1], res, i1);
  ScanConvertPolyData(p2, bb[0], bb[1], res, i2);

  // Compute the volume overlap between images
  unsigned long nUnion = 0, nIntersect = 0, nReference = 0, nModel = 0;
  typedef ImageRegionConstIterator<ByteImageType> IterType;
  IterType it1(i1, i1->GetBufferedRegion());
  IterType it2(i2, i2->GetBufferedRegion());
  for( ; !it1.IsAtEnd() ; ++it1, ++it2)
    {
    if(it1.Value() > 0 && it2.Value() > 0)
      nIntersect++;
    if(it1.Value() > 0 || it2.Value() > 0)
      nUnion++;
    if(it2.Value() > 0)
      nReference++;
    if(it1.Value() > 0)
      nModel++;
    }

  // Compute the distance transform for each image
  FloatImageType::Pointer d1 = FloatImageType::New();
  FloatImageType::Pointer d2 = FloatImageType::New();
  ComputeDistanceTransform(i1, d1);
  ComputeDistanceTransform(i2, d2);

  // Integrate the distance transform over each mesh
  double xOneNorm[2], xTwoNorm[2], xInfNorm[2];
  IntegrateDistance(p1, d2, xOneNorm[0], xTwoNorm[0], xInfNorm[0]);
  IntegrateDistance(p2, d1, xOneNorm[1], xTwoNorm[1], xInfNorm[1]);

  // Report our findings
  cout << "Unbiased overlap between meshes is " << nIntersect * 1.0 / nUnion << endl;
  cout << "Biased overlap between meshes is " << nIntersect * 1.0 / nReference << endl;
  cout << "Reverse biased overlap between meshes is " << nIntersect * 1.0 / nModel << endl;
  cout << "Distance from mesh 1 to mesh 2 is " << xOneNorm[0] << "; " << xTwoNorm[0] << "; " << xInfNorm[0] << endl;
  cout << "Distance from mesh 2 to mesh 1 is " << xOneNorm[1] << "; " << xTwoNorm[1] << "; " << xInfNorm[1] << endl;
}
