/**
 * Program Objective: to compute the shape operator on the implicit surface
 * defined by the given floating point image. The shape operator can be used to
 * derive principal curvatures and principal directions on the surface. 
 */

#include "ReadWriteImage.h"
#include <iostream>
#include <string>

#include <itkConstNeighborhoodIterator.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>

#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>

#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>

using namespace std;

int usage()
{
  cout << "usage sulcalfilter [options] input.img output_base_name \n" << endl;
  cout << "options : " << endl;
  cout << "  -m in.vtk out.vtk           Interpolate the curvature at vertices in a mesh " << endl;
  return -1;  
}
  
// Define the matrix pixel type
typedef vnl_vector_fixed<float, 3> Vec;
typedef vnl_matrix_fixed<float, 3, 3> Mat;
  
bool ComputePrincipalCurvatures(
  const Mat &xHessian, const Vec &xGradient, 
  float &xCurvMinor, float &xCurvMajor)
{
  // Compute the Shape operator, given by 
  // (I - (G G^T) / |G|^2) Hess[F] / |G| where G is grad F
  // this formula taken from Eberly's paper
  // http://www.magic-software.com/Documentation/PrincipalCurvature.pdf
  double xGradMagSqr = xGradient.squared_magnitude();

  // The computation is possible only if this is a non-zero quantity
  // because we have to divide by the gradient magnitude
  if(xGradMagSqr == 0.0f) return false;

  // The shape operator matrix
  Mat xShape;

  // The factor by which everything is multiplied (grad mag to -3/2 power)
  double xFactor = 1.0  / (xGradMagSqr * sqrt(xGradMagSqr));

  // The actual gradient magnitude
  double xGradMag = sqrt(xGradMagSqr);

  // Compute the elements of the shape operator
  for(unsigned int i=0; i < 3; i++)
    {
    xShape[i][i] = static_cast<float>(
      xFactor * xHessian[i][i] * 
      (xGradMagSqr - xGradient[i] * xGradient[i]));

    for(unsigned int j=0; j < i; j++)
      {
      xShape[i][j] = xShape[j][i] = static_cast<float>(
        xFactor * xHessian[i][j] *
        (xGradMagSqr - xGradient[i] * xGradient[j]));
      }
    }

  // Now, compute the principal curvatures by solving the 3x3 eigensystem.
  vnl_symmetric_eigensystem<float> eig(xShape);

  // One of the eigenvalues is necessarily zero and should be ignored
  if(eig.get_eigenvalue(0) == 0.0f)
    {
    xCurvMinor = eig.get_eigenvalue(1); xCurvMajor = eig.get_eigenvalue(2);
    }
  else if(eig.get_eigenvalue(1) == 0.0f)
    {
    xCurvMinor = eig.get_eigenvalue(0); xCurvMajor = eig.get_eigenvalue(2);
    }
  else
    {
    xCurvMinor = eig.get_eigenvalue(0); xCurvMajor = eig.get_eigenvalue(1);
    }

  return true;
}

int main(int argc, char *argv[])
{
  // Usage
  if(argc < 3) return usage();

  // Process the options
  const char *fnInput = argv[argc-2];
  string sOutBase = argv[argc-1];
  bool flagProcessMesh = false;
  const char *fnMeshInput, *fnMeshOutput;

  // Parse the command line arguments
  for(unsigned int iArg = 1; iArg < argc - 2; iArg++)
    {
    if(!strcmp(argv[iArg], "-m"))
      {
      flagProcessMesh = true;
      fnMeshInput = argv[++iArg];
      fnMeshOutput = argv[++iArg];
      }
    }

  // Define the float based image
  typedef itk::Image<float,3> ImageType;
  ImageType::Pointer imgInput;

  // Read the image
  ReadImage(imgInput,fnInput);

  // Create the necessary image types
  typedef itk::Image< Mat, 3> MatImage;
  typedef itk::Image< Vec, 3> VecImage;
  typedef itk::Image< float, 3> FloatImage;

  // Make an array of curvature images
  enum CurvatureType { MINOR = 0, MAJOR, GAUSS, MEAN, COUNT };
  FloatImage::Pointer imgCrv[4];
  for(unsigned int iCurv = 0; iCurv < COUNT; iCurv++)
    {
    imgCrv[iCurv] = FloatImage::New();
    imgCrv[iCurv]->SetRegions(imgInput->GetBufferedRegion());
    imgCrv[iCurv]->Allocate();
    imgCrv[iCurv]->FillBuffer(0.0f);
    }

  // Initialize the 9-dimensional geometry image that will hold
  // the gradient and Hessian information
  typedef itk::FixedArray<float, 9> GeomArray;
  typedef itk::Image< GeomArray, 3> GeomImage;
  
  GeomImage::Pointer imgGeom = GeomImage::New();
  imgGeom->SetRegions(imgInput->GetBufferedRegion());
  imgGeom->Allocate();

  // Initialize the necessary images

  // Set all images to dimensions of the input
  // MatImage::Pointer imgHessian = MatImage::New();
  // imgHessian->SetRegions(imgInput->GetBufferedRegion());
  // imgHessian->Allocate();
  // imgHessian->FillBuffer(Mat(0.0f));

  VecImage::Pointer imgGradient = VecImage::New();
  imgGradient->SetRegions(imgInput->GetBufferedRegion());
  imgGradient->Allocate();
  imgGradient->FillBuffer(Vec(0.0f));

  // Create a neighborhood iterator through the input image
  typedef itk::ConstNeighborhoodIterator<FloatImage> SourceIterator;
  itk::Size<3> radius = {{1,1,1}};
  SourceIterator it(radius, imgInput, imgInput->GetBufferedRegion());

  // Precompute the 'strides' for accessing neighbors
  unsigned int stride[3], center;
  center = it.Size() >> 1;
  stride[0] = it.GetStride(0);
  stride[1] = it.GetStride(1);
  stride[2] = it.GetStride(2);

  // Iterate through the image, computing the gradient image. This is the 
  // first step, after which we will compute the second derivatives. This is
  // faster than using the radius 2 neighborhood.
  while(!it.IsAtEnd())
    {
    // Define the gradient vector and the hessian matrix
    Vec xGradient(0.0f);

    // Get the current pixel value
    float xCenter = it.GetPixel( center );

    // Compute the Hessian matrix
    for(unsigned int i = 0; i < 3; i++)
      {
      // Get the neighboring values in the current direction
      float xNext = it.GetPixel( center + stride[i] );
      float xPrev = it.GetPixel( center - stride[i] );

      // Average the values to compute directional derivative using
      // central differences. (Argh!)
      xGradient[i] = 0.5f * ( xNext - xPrev);
      }

    // Store the gradient
    imgGradient->SetPixel(it.GetIndex(), xGradient);

    // Go to the next pixel
    ++it;
    }

  cout << " Computed the gradient image " << endl;

  // Create a neighborhood iterator for the gradient image. 
  typedef itk::ConstNeighborhoodIterator<VecImage> GradientIterator;
  GradientIterator itGrad(radius, imgGradient, imgGradient->GetBufferedRegion());

  // Compute the Hessian matrix using central differences in the gradient image
  while(!itGrad.IsAtEnd())
    {
    // Matrix to store the Hessian
    Mat xHessian(0.0f);
    
    // Get the gradient at the center position
    const Vec &xCenterGrad = itGrad.GetPixel( center );

    // Compute each of the Hessian components
    for(unsigned int i = 0; i < 3 ; i++)
      {
      // Get the previous and the next gradients
      const Vec &xNextGrad = itGrad.GetPixel( center + stride[i] );
      const Vec &xPrevGrad = itGrad.GetPixel( center - stride[i] );

      // Compute the second derivative in the i-th direction
      xHessian[i][i] = 0.5f * (xNextGrad[i] - xPrevGrad[i]);

      // Compute the mixed derivatives
      for(unsigned int j = 0; j < i; j++)
        {
        // Get the adjacent gradients in the J direction
        const Vec &xNextJGrad = itGrad.GetPixel( center + stride[j] );
        const Vec &xPrevJGrad = itGrad.GetPixel( center - stride[j] );

        // Compute the second derivatives
        xHessian[i][j] = xHessian[j][i] = 
          0.25 * (
            xNextJGrad[i] - xPrevJGrad[i] + 
            xNextGrad[j] - xPrevGrad[j]);
        }
      }

    // Create the Geometry information vector
    GeomArray G;
    G[0] = xCenterGrad[0];
    G[1] = xCenterGrad[1];
    G[2] = xCenterGrad[2];
    G[3] = xHessian[0][0];
    G[4] = xHessian[1][1];
    G[5] = xHessian[2][2];
    G[6] = xHessian[0][1];
    G[7] = xHessian[0][2];
    G[8] = xHessian[1][2];

    // Store the Hessian for later use
    imgGeom->SetPixel(itGrad.GetIndex(), G);

    // Compute the curvatures
    float xCurvMajor = 0.0, xCurvMinor = 0.0;
    if(ComputePrincipalCurvatures(xHessian,xCenterGrad,xCurvMinor,xCurvMajor))
      {
      // Save the curvatures
      imgCrv[MINOR]->SetPixel(itGrad.GetIndex(),xCurvMinor);
      imgCrv[MAJOR]->SetPixel(itGrad.GetIndex(),xCurvMajor);

      // Compute mean and gaussian curvatures
      imgCrv[GAUSS]->SetPixel(itGrad.GetIndex(),xCurvMinor * xCurvMajor);
      imgCrv[MEAN]->SetPixel(itGrad.GetIndex(),0.5 * (xCurvMinor + xCurvMajor));
      }
    
    // Go to the next pixel
    ++itGrad;
    }

  // Store the curvature files
  WriteImage(imgCrv[MAJOR], string( sOutBase + ".pc1.img.gz").c_str() );
  WriteImage(imgCrv[MINOR], string( sOutBase + ".pc2.img.gz").c_str() );
  WriteImage(imgCrv[GAUSS], string( sOutBase + ".gauss.img.gz").c_str() );
  WriteImage(imgCrv[MEAN], string( sOutBase + ".mean.img.gz").c_str() );

  // Process the mesh if requested
  if(flagProcessMesh)
    {
    unsigned int i;

    // Load the mesh from VTK file
    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    reader->SetFileName(fnMeshInput);
    reader->Update();

    // Get the polygon data
    vtkPolyData *mesh = reader->GetOutput();

    // Create a linear interpolator for the geometry data
    typedef itk::VectorLinearInterpolateImageFunction<GeomImage> Vecolator;
    Vecolator::Pointer fVecInterp = Vecolator::New();
    fVecInterp->SetInputImage(imgGeom);

    // Create data arrays for curvatures
    vtkFloatArray *xCrvArray[4];
    const char *names[] = { "MAJOR", "MINOR", "GAUSS", "MEAN" };
    for(i=0; i < 4; i++) 
      {
      xCrvArray[i] = vtkFloatArray::New();
      xCrvArray[i]->SetName(names[i]);
      // xCrvArray[i]->SetNumberOfComponents(1);
      // xCrvArray[i]->SetNumberOfTuples(mesh->GetNumberOfPoints());
      // xCrvArray[i]->Allocate(mesh->GetNumberOfPoints());
      }

    // Iterate through the mesh
    for(vtkIdType iVertex = 0; iVertex < mesh->GetNumberOfPoints(); iVertex++)
      {
      // Create an index structure
      double *x = mesh->GetPoint(iVertex);
      Vecolator::ContinuousIndexType idx;
      idx[0] = (double) x[0]; 
      idx[1] = (double) x[1];
      idx[2] = (double) x[2];
      
      // Interpolate the geometry data
      Vecolator::OutputType G = fVecInterp->EvaluateAtContinuousIndex(idx);

      // Assign the gradient info
      Vec xGradient; 
      xGradient[0] = G[0];
      xGradient[1] = G[1];
      xGradient[2] = G[2];

      // Assign the Hessian info
      Mat xHessian;
      xHessian[0][0] = G[3];
      xHessian[1][1] = G[4];
      xHessian[2][2] = G[5];
      xHessian[0][1] = xHessian[1][0] = G[6];
      xHessian[0][2] = xHessian[2][0] = G[7];
      xHessian[1][2] = xHessian[2][1] = G[8];

      // Compute the curvatures 
      float xCurvMajor = 0.0f, xCurvMinor = 0.0f;
      ComputePrincipalCurvatures(xHessian,xGradient,xCurvMajor,xCurvMinor);

      // Assign curvatures to data arrays 
      xCrvArray[MAJOR]->InsertNextValue(xCurvMajor);
      xCrvArray[MINOR]->InsertNextValue(xCurvMinor);

      xCrvArray[GAUSS]->InsertNextValue(xCurvMajor * xCurvMinor);
      xCrvArray[MEAN]->InsertNextValue(0.5 * (xCurvMinor + xCurvMajor));
      }

    // Store the arrays with the data
    for(i=0; i < 4; i++) 
      mesh->GetPointData()->AddArray(xCrvArray[i]);

    cout << "array " << mesh->GetPointData()->GetArray(3)->GetNumberOfTuples() << endl;

    // Save the array
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInput(mesh);
    writer->SetFileName(fnMeshOutput);
    writer->SetFileTypeToBinary();
    writer->Update();
    }
}
