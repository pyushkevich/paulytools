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
  cout << "meshtoimg - fill in a 3D mesh, producing a binary image" << endl;
  cout << "description:" << endl;
  cout << "        Given an 3D mesh, this command fill scan-convert the faces in " << endl;
  cout << "   the mesh to a 3D image, optionally filling the interior (-f option)"  << endl;
  cout << "   The user must specify the corners of the image in the coordinate"  << endl;
  cout << "   space of the mesh using the -o and -s options, as well as the resolution"  << endl;
  cout << "   of the image using the -r option. " << endl;
  cout << "        To find out the spatial extent of the mesh, call the program without" << endl;
  cout << "   the -o and -s parameters. The extent will be printed out." << endl;
  cout << "usage: " << endl;
  cout << "   stltoimg [options] output.img" << endl;
  cout << "options: " << endl;
  cout << "   -stl file                   Specify input as an STL mesh" << endl;
  cout << "   -byu file                   Specify input as a BYU mesh" << endl;
  cout << "   -f                          Fill the interior of the mesh also" << endl;
  cout << "   -o NN NN NN                 Image origin in x,y,z in mesh space" << endl;
  cout << "   -s NN NN NN                 Image size in x,y,z in mesh space" << endl;
  cout << "   -r NN NN NN                 Image resolution in x,y,z in pixels" << endl;
  cout << "   -v                          Visualize (render) the STL mesh on the screen" << endl;
  return -1;
}

void drawPolyData(vtkPolyData *poly)
{
  vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
  mapper->SetInput(poly);

  vtkLODActor *actor = vtkLODActor::New();
  actor->SetMapper(mapper);

  vtkRenderer *ren = vtkRenderer::New();
  ren->AddActor(actor);
  ren->SetBackground(0.1,0.2,0.4);

  vtkRenderWindow *renWin = vtkRenderWindow::New();
  renWin->AddRenderer(ren);
  renWin->SetSize(500,500);

  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);
  iren->Initialize();

  renWin->Render();
  iren->Start();

  iren->Delete();
  renWin->Delete();
  ren->Delete();
  actor->Delete();
  mapper->Delete();
}

int main(int argc, char **argv)
{
  // Parameter values
  double bb[2][3] = {{0.,0.,0.},{128.,128.,128.}};
  int res[3] = {128,128,128};
  bool doRender = false, doFloodFill = false, flagInputGiven = false;

  // Input specifications
  enum MeshType {NONE, BYU, STL} iMeshType = NONE;
  string fnOutput, fnInput;

  // Check the parameters
  if(argc < 2) return usage();

  // Get the file names
  fnOutput = argv[argc-1];

  // Parse the optional parameters
  try 
    {
    for(unsigned int i=1;i<argc-2;i++)
      {
      string arg = argv[i];
      if(arg == "-stl")
        {
        fnInput = argv[++i];
        iMeshType = STL;
        }
      else if(arg == "-byu")
        {
        fnInput = argv[++i];
        iMeshType = BYU;
        }
      else if(arg == "-o")
        {
        bb[0][0] = atof(argv[++i]);
        bb[0][1] = atof(argv[++i]);
        bb[0][2] = atof(argv[++i]);
        }
      else if(arg == "-s")
        {
        bb[1][0] = atof(argv[++i]);
        bb[1][1] = atof(argv[++i]);
        bb[1][2] = atof(argv[++i]);
        }
      else if(arg == "-r")
        {
        res[0] = atoi(argv[++i]);
        res[1] = atoi(argv[++i]);
        res[2] = atoi(argv[++i]);
        }      
      else if(arg == "-v")
        {
        doRender = true;
        }
      else if(arg == "-f")
        {
        doFloodFill = true;
        }
      else
        {
        cout << "Unknown argument " << arg << endl;
        return usage();
        }
      }
    }
  catch(...) 
    {
    cout << "Error parsing command line options!" << endl;
    return usage();
    }

  // If no file given, return
  if(iMeshType == NONE)
    {
    cerr << "No input mesh specified!" << endl;
    return usage();
    }    	

  // Print the parameters
  cout << endl << "Parameters:" << endl;
  cout << "   Image resolution  : " << res[0] << " by " << res[1] << " by " << res[2] << endl;
  cout << "   Image box origin  : {" << bb[0][0] << ", " << bb[0][1] << ", " << bb[0][2] << "}" 
    << endl;
  cout << "   Image box size    : {" << bb[1][0] << ", " << bb[1][1] << ", " << bb[1][2] << "}" 
    << endl << endl;

  // The real program begins here

  // Read the appropriate mesh type
  vtkPolyData *polyInput = NULL;
  vtkObject *fltGenericReader = NULL;

  if(iMeshType == BYU)
    {
    vtkBYUReader *reader = vtkBYUReader::New();
    reader->SetFileName(fnInput.c_str());
    cout << "Reading the BYU file ..." << endl;
    reader->Update();

    fltGenericReader = reader;
    polyInput = reader->GetOutput();
    }
  else if(iMeshType == STL)
    {
    vtkSTLReader *reader = vtkSTLReader::New();
    reader->SetFileName(fnInput.c_str());
    cout << "Reading the STL file ..." << endl;
    reader->Update();

    fltGenericReader = reader;
    polyInput = reader->GetOutput();
    }

  // Convert the model to triangles
  vtkTriangleFilter *tri = vtkTriangleFilter::New();
  tri->SetInput(polyInput);
  cout << "Converting to triangles ..." << endl;
  tri->Update();
  vtkPolyData *pd = tri->GetOutput();

  // Display the polydata that we loaded
  if(doRender)
    drawPolyData(pd);

  // Get the extents of the data
  double *bounds = pd->GetBounds();  

  cout << "STL mesh bounds: " << endl;
  cout << "   X : " << bounds[0] << " to " << bounds[1] << endl;
  cout << "   Y : " << bounds[2] << " to " << bounds[3] << endl;
  cout << "   Z : " << bounds[4] << " to " << bounds[5] << endl;

  if(bounds[0] < bb[0][0] || bounds[1] > bb[0][0] + bb[1][0] ||
    bounds[2] < bb[0][1] || bounds[3] > bb[0][1] + bb[1][1] ||
    bounds[4] < bb[0][2] || bounds[5] > bb[0][2] + bb[1][2])
    {
    cout << "User specified bounds (-o -s) are out of range! Can't continue!" << endl;
    return -1;
    }

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
  typedef Image<unsigned char,3> ImageType;
  ImageType::Pointer img = ImageType::New();
  ImageType::RegionType region;
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
  ImageType::SpacingType xSpacing;
  ImageType::PointType xOrigin;
  for(unsigned int d = 0; d < 3; d++)
    {
    xOrigin[d] = bb[0][d];
    xSpacing[d] = bb[1][d] / res[d];
    }
  img->SetOrigin(xOrigin);
  img->SetSpacing(xSpacing);

  // Save the image to disk
  typedef ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(img);
  writer->SetFileName(fnOutput.c_str());

  cout << "Writing the output image ..." << endl;
  writer->Update();
}
