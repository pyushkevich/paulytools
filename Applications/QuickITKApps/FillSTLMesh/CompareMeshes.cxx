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
#include <vtkBYUReader.h>
#include <vtkSTLReader.h>
#include <vtkTriangleFilter.h>

#include <itkImageFileWriter.h>
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

  // Read the model
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(fn.c_str());
  cout << "Reading the BYU file ..." << endl;
  reader->Update();
  p1 = reader->GetOutput();

  // Convert the model to triangles
  vtkTriangleFilter *tri = vtkTriangleFilter::New();
  tri->SetInput(p1);
  cout << "Converting to triangles ..." << endl;
  tri->Update();
  vtkPolyData *pd = tri->GetOutput();

  return pd;
}

void ScanConvertPolyData(vtkPolyData *pd, double **bb, int *res, ByteImageType *img)
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
        vtx[it][j] = res[j] * (x[j] - bb[0][j]) / bb[1][j];
        }

      ++it;
      }      
    }

  // Clean up
  fltGenericReader->Delete();

  // Create a ITK image to store the results
  ByteImageType::RegionType region;
  region.SetSize(0,res[0]); region.SetSize(1,res[1]); region.SetSize(2,res[2]);
  img->SetRegions(region);
  img->Allocate();
  img->FillBuffer(0);

  // Convert the polydata to an image
  if(doFloodFill)
    {
    cout << "Scan converting triangles and filling the interior ..." << endl;
    drawBinaryTrianglesFilled(img->GetBufferPointer(), res, vtx, nt);
    }
  else
    {
    cout << "Scan converting triangles ..." << endl;
    drawBinaryTrianglesSheetFilled(img->GetBufferPointer(), res, vtx, nt);
    }  

  // Set the origin and spacing of the image
  ByteImageType::SpacingType xSpacing;
  ByteImageType::PointType xOrigin;
  for(unsigned int d = 0; d < 3; d++)
    {
    xOrigin[d] = bb[0][d];
    xSpacing[d] = bb[1][d] / res[d];
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
    bb[i][0] = (b1[2*i] < b2[2*i] ? b1[2*i] : b2[2*i]) - 4;
    bb[i][1] = (b1[2*i+1] < b2[2*i+1] ? b1[2*i+1] : b2[2*i+1]) + 4;
    cout << "Bounds[" << i << "]: " << bb[i][0] << " to " << bb[i][1] << endl;
    res[i] = ceil((bb[i][1] - bb[i][0]) / xVox);
    }

  // Save the image to disk
  typedef ImageFileWriter<ByteImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(img);
  writer->SetFileName(fnOutput.c_str());

  cout << "Writing the output image ..." << endl;
  writer->Update();
}
