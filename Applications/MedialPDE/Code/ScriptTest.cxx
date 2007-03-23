#include "ScriptInterface.h"
#include "BasisFunctions2D.h"
#include "MedialAtom.h"
#include "CartesianMedialModel.h"
#include "OptimizationTerms.h"
#include "CoefficientMapping.h"
#include "MedialAtomGrid.h"
#include "PrincipalComponents.h"
#include "TestSolver.h"
#include "vnl/vnl_erf.h"
#include "vnl/vnl_random.h"

#include "vtkOBJReader.h"
#include "vtkBYUWriter.h"
#include "vtkPolyData.h"

#include <string>
#include <iostream>

using namespace std;
using namespace medialpde;

string dirWork = "/home/pauly/data2005/Stanley/data/";
// string dirWork = "/mnt/data2/PUBLIC/Data/Input/StanleySchizophrenia/";

/**
 * This is a test class that is a sphere that pretends to be a floating
 * point image. The idea is to use this ellipsoid for testing image-based
 * derivative code, as the derivatives computed here are exact rather than
 * numerical.
 */
class TestFloatImage : public FloatImage
{
public:
  TestFloatImage(const SMLVec3d &C, double R, double sigma)
    {
    this->C = C;
    this->R = R;
    this->sigma = sigma;
    } 

  // Interpolate the image at a given position
  float Interpolate(const SMLVec3d &x)
    { 
    double u = ( R - (x - C).two_norm() ) / sigma;
    return vnl_erf(u);
    }
    
  // Interpolate the image gradient
  void InterpolateImageGradient(const SMLVec3d &x, SMLVec3f &g)
    {
    SMLVec3d g1;
    InterpolateImageGradient(x, g1);
    g[0] = (float) g1[0];
    g[1] = (float) g1[1];
    g[2] = (float) g1[2];
    }

  // Interpolate the image gradient
  void InterpolateImageGradient(const SMLVec3d &x, SMLVec3d &g)
    {
    double z = (x - C).two_norm();
    double u = ( R - z ) / sigma;
    double a = 2.0 * exp(-u * u) / sqrt(M_PI);
    
    g[0] = a * (- (x[0] - C[0]) / z) / sigma;
    g[1] = a * (- (x[1] - C[1]) / z) / sigma;
    g[2] = a * (- (x[2] - C[2]) / z) / sigma;
    }

  // Integrate positive voxels in the image volume
  double IntegratePositiveVoxels()
    {
    // This is a hack!
    return 4.0 * M_PI * R * R * R / 3.0;
    }

private:
  SMLVec3d C;
  double R, sigma;
};



class TestFunction01 : public EuclideanFunction
{
public:
  double Evaluate(const SMLVec3d &x)
    { return 1; }
};

// Test volume computations
void TestAreaAndVolume(GenericMedialModel *xSolver)
{
  unsigned int nCuts = 6;

  // Compute the area 
  MedialIterationContext *xGrid = xSolver->GetIterationContext();
  double xArea, *xAreaWeights = new double[xGrid->GetNumberOfBoundaryPoints()];
  xArea = ComputeMedialBoundaryAreaWeights(xGrid, xSolver->GetAtomArray(), xAreaWeights);

  cout << "Area : " << xArea << endl;
  TestFunction01 fn1;
  cout << "Verification " << 
    IntegrateFunctionOverBoundary(xGrid, xSolver->GetAtomArray(), 
      xAreaWeights, &fn1) << endl;
  
  // Compute the volume
  double xVol, *xVolWeights = new double[xGrid->GetNumberOfInternalPoints(nCuts)];
  double *xProfileWeights = new double[xGrid->GetNumberOfProfileIntervals(nCuts)];
  
  SMLVec3d *xInternal = new SMLVec3d[xGrid->GetNumberOfInternalPoints(nCuts)];
  ComputeMedialInternalPoints(
    xGrid, xSolver->GetAtomArray(), nCuts, xInternal);
  
  xVol = ComputeMedialInternalVolumeWeights(
    xGrid, xInternal, nCuts, xVolWeights, xProfileWeights);

  cout << "Volume : " << xVol << endl;
  cout << "Verification " << 
    IntegrateFunctionOverInterior(
      xGrid, xInternal, xVolWeights, nCuts, &fn1) << endl;

  // Volume accumulator
  xVol = 0.0;
  
  // Create an internal point iterator
  for(MedialBoundaryPointIterator it(xGrid); !it.IsAtEnd(); ++it) 
    {
    // Evaluate the image at this location
    BoundaryAtom &B = GetBoundaryPoint(it, xSolver->GetAtomArray()); 
    xVol += dot_product(B.X , B.N) * xAreaWeights[it.GetIndex()] / 3.0; 
    }

  // Return match scaled by total weight
  cout << "Verification " << xVol << endl;

  // One more test
  xVol = 0.0;
  for(MedialProfileIntervalIterator itProf(xGrid,nCuts); !itProf.IsAtEnd(); ++itProf)
    xVol += xProfileWeights[itProf.GetIndex()];
  cout << "Profiles add up to " << xVol << endl;

  // Check the Jacobian
  SolutionData S(xGrid, xSolver->GetAtomArray());
  BoundaryJacobianEnergyTerm termJac;
  termJac.ComputeEnergy(&S);
  termJac.PrintReport(cout);
}

void ExportMedialMeshToVTK(
  GenericMedialModel *xModel, ITKImageWrapper<float> *xImage, const char *file);

void ExportBoundaryMeshToVTK(
  GenericMedialModel *xModel, ITKImageWrapper<float> *xImage, const char *file);

/**
 * This test ensures that all the volume elements in a model are computed
 * correctly.
 */
int TestVolumeComputation(const char *file)
{
  // The number of cuts at which to test
  size_t i,j,k,nCuts = 5;
 
  // Success code
  bool flagSuccess = true;

  // Read the model
  MedialPDE mp(file);
  GenericMedialModel *model = mp.GetMedialModel();

  // Create a solution data based on this model
  MedialIterationContext *grid = model->GetIterationContext();
  SolutionData sd(grid, model->GetAtomArray());

  // Generate internal weights in the solution
  sd.UpdateInternalWeights(nCuts);
  sd.UpdateBoundaryWeights();

  // Test 1. Ensure that the volume of every wedge in the model is positive
  double minWeight = sd.xInternalWeights[0], maxWeight = sd.xInternalWeights[0];
  for(i = 1; i < sd.nInternalPoints; i++)
    {
    minWeight = std::min(minWeight, sd.xInternalWeights[i]);
    maxWeight = std::max(maxWeight, sd.xInternalWeights[i]);
    }
  cout << "Max volume element: " << maxWeight << endl;
  cout << "Min volume element: " << minWeight << endl;
  cout << "Avg volume element: " << (sd.xInternalVolume / sd.nInternalPoints) << endl;

  if(minWeight < 0.0)
    flagSuccess = false;

  // Test 2. Ensure that the boundary area elements are also positive
  double minAreaElt = sd.xBoundaryWeights[0], 
         maxAreaElt = sd.xBoundaryWeights[0],
         sumAreaElt = sd.xBoundaryWeights[0];
  for(i = 1; i < model->GetNumberOfBoundaryPoints(); i++)
    {
    minAreaElt = std::min(minAreaElt, sd.xBoundaryWeights[i]);
    maxAreaElt = std::max(minAreaElt, sd.xBoundaryWeights[i]);
    sumAreaElt += sd.xBoundaryWeights[i];
    }
  cout << "Max area element: " << maxAreaElt << endl;
  cout << "Min area element: " << minAreaElt << endl;
  cout << "Avg area element: " << 
    (sd.xBoundaryArea / model->GetNumberOfBoundaryPoints()) << endl;
  cout << "Sum area element: " << sumAreaElt << endl;
  cout << "Surface area: " << sd.xBoundaryArea << endl;

  if(minAreaElt < 0.0 || sumAreaElt != sd.xBoundaryArea)
    flagSuccess = false;

  // Test 3. Compare volume computed by integration to volume estimated by
  // Green's theorem
  double xVolGreen = 0.0;
  for(MedialBoundaryPointIterator bit(grid); !bit.IsAtEnd(); ++bit)
    {
    BoundaryAtom bat = GetBoundaryPoint(bit, model->GetAtomArray());
    xVolGreen += dot_product(bat.X, bat.N) * sd.xBoundaryWeights[bit.GetIndex()];
    }
  xVolGreen /= 3.0;

  cout << "Volume Integral : " << sd.xInternalVolume << endl;
  cout << "Green Thm. Integral : " << xVolGreen << endl;

  // Test 4. Another way to compute Green's integral
  double xVolGreen2 = 0.0;
  for(MedialBoundaryTriangleIterator btt(grid); !btt.IsAtEnd(); ++btt)
    {
    // Get the triangle
    SMLVec3d A = GetBoundaryPoint(btt, model->GetAtomArray(), 0).X;
    SMLVec3d B = GetBoundaryPoint(btt, model->GetAtomArray(), 1).X;
    SMLVec3d C = GetBoundaryPoint(btt, model->GetAtomArray(), 2).X;

    // Get the area-scaled normal vector
    SMLVec3d N = vnl_cross_3d(B-A,C-A) / 2.0;
    SMLVec3d X = (A + B + C) / 3.0;
    xVolGreen2 += dot_product(X, N) / 3.0;
    }
  cout << "Green Thm. Integral 2 : " << xVolGreen2 << endl;

  // Export model for external verification
  ExportBoundaryMeshToVTK(model, NULL, "volume_test.vtk");

  // Return the test score
  return flagSuccess ? 0 : -1;
}

void Test01()
{
  CartesianMPDE *mp = new CartesianMPDE(3, 5, 20, 40);
  //mp->LoadFromParameterFile(fMrep.c_str());
  
  mp->LoadFromDiscreteMRep("/tmp/surf01.txt",-0.3);
  // mp->GenerateSampleModel();
  mp->Solve();
  mp->SaveBYUMesh("temp.byu");

  // Make sure areas and volumes add up
  TestAreaAndVolume(mp->GetMedialModel());
  
  // Load the image and gradients
  // FloatImage img;
  // img.LoadFromFile((dirWork + "avg/average_hippo_blurred_hi.mha").c_str());

  // Match the volume to the image
  // mp->MatchImageByMoments(&img, 5);

  
  
  // RenderMedialPDE(mp);
}

/**
void Test02()
{
  // Decide how many points to interpolate
  unsigned int nSide = 11, nPoints = nSide * nSide;
  double uStep = 1.0 / (nSide - 1);

  // Allocate arrays of points and coordinates
  double xPoints[nPoints], yPoints[nPoints], zPoints[nPoints];
  double uPoints[nPoints], vPoints[nPoints];

  // Create an array of points
  unsigned int i = 0;
  for(unsigned int u = 0; u < nSide; u++)
    for(unsigned int v = 0; v < nSide; v++)
      {
      double uu = u * uStep, vv = v * uStep;

      uPoints[i] = uu;
      vPoints[i] = vv;
      xPoints[i] = 0.5 * uu + 0.25;
      yPoints[i] = vv;
      zPoints[i] = ((uu - 0.5) * (uu - 0.5) + (vv - 0.5) * (vv - 0.5)) * 0.25;
      ++i;
      }

  // Peform the fit
  FourierSurfaceOld s1(5, 5, 3);
  s1.FitData(0, nPoints, uPoints, 1, vPoints, 1, xPoints, 1);
  s1.FitData(1, nPoints, uPoints, 1, vPoints, 1, yPoints, 1);
  s1.FitData(2, nPoints, uPoints, 1, vPoints, 1, zPoints, 1);

  FourierSurface s2(5, 5);
  s2.FitToData(nPoints, 0, uPoints, vPoints, xPoints);
  s2.FitToData(nPoints, 1, uPoints, vPoints, yPoints);
  s2.FitToData(nPoints, 2, uPoints, vPoints, zPoints);

  SMLVec3d T1, T2;
  for(double u = 0.0; u <= 1.0; u+=0.234)
    for(double v = 0.0; v <= 1.0; v+=0.1234)
      {
      T1[0] = s1.Evaluate(u, v, 1, 0, 0);
      T1[1] = s1.Evaluate(u, v, 1, 0, 1);
      T1[2] = s1.Evaluate(u, v, 1, 0, 2);
      s2.EvaluateDerivative(u, v, 1, 0, 0, 3, T2.data_block());
      cout << T1 - T2 << endl;
      }
  
}
*/

void Test03()
{
  // Create a medial PDE object
  CartesianMPDE *mp = new CartesianMPDE(5, 5, 50, 100);
  mp->LoadFromParameterFile((dirWork + "avg/average_mrepL_01.mpde").c_str());
  mp->Solve();

  // Convert the binary image into a floating point and compute match
  FloatImage *img = new FloatImage();
  img->LoadFromFile((dirWork + "avg/average_hippo_blurred.mha").c_str());
  img->LoadGradientFromFile(0, (dirWork + "avg/average_hippo_blurred.dx.mha").c_str());
  img->LoadGradientFromFile(1, (dirWork + "avg/average_hippo_blurred.dy.mha").c_str());
  img->LoadGradientFromFile(2, (dirWork + "avg/average_hippo_blurred.dz.mha").c_str());
  img->SetOutsideValue(-1.0);

  // Compute the gradient of the match
  cout << "Image Match = " << mp->ComputeImageMatch(img) << endl;

  // TODO: Load optimization parameters from a file
  mp->RunOptimization(img, 30, "file.txt");
  mp->SaveToParameterFile((dirWork + "avg/average_mrepL_02.mpde").c_str());
  mp->SaveBYUMesh((dirWork + "avg/average_mrepL_02.byu").c_str());
}

int TestWedgeVolume()
{
  // Define an object with a volume that we know
  double vol1 = WedgeVolume(
    SMLVec3d(0, 0, 0), SMLVec3d(1, 0, 0), SMLVec3d(1, 1, 0),
    SMLVec3d(0, 0, 1), SMLVec3d(1, 0, 1), SMLVec3d(1, 1, 1));

  cout << "Volume 1 = " << vol1 << " and should be 0.5" << endl;
  
  // Compare to the right amount
  return (vol1 == 0.5) ? 0 : -1;
}

void TestCellVolume()
{
  // Compute the volume of a unit cube
  cout << "Unit Cube Volume: " << CellVolume(
    SMLVec3d(0,0,0),SMLVec3d(0,0,1),SMLVec3d(0,1,0),SMLVec3d(0,1,1),
    SMLVec3d(1,0,0),SMLVec3d(1,0,1),SMLVec3d(1,1,0),SMLVec3d(1,1,1)) << endl;

  // Compute the volume of an elongated cube
  cout << "123 Cube Volume : " << CellVolume(
    SMLVec3d(0,0,0),SMLVec3d(0,0,3),SMLVec3d(0,2,0),SMLVec3d(0,2,3),
    SMLVec3d(1,0,0),SMLVec3d(1,0,3),SMLVec3d(1,2,0),SMLVec3d(1,2,3)) << endl;
}


/** Test differential geometry relationships */
int TestDifferentialGeometry(const char *fnMPDE)
{
  // Load a medial PDE for the test
  CartesianMPDE mp(2, 4, 25, 49);
  mp.LoadFromParameterFile(fnMPDE);

  // Pick a point to evaluate at
  double u = 0.6, v = 0.2, eps = 0.0001;

  // Compute the jet at this point
  SMLVec3d X00, X10, X01, X20, X11, X02;
  mp.GetMedialSurface()->EvaluateDerivative(u, v, 0, 0, 0, 3, X00.data_block());
  mp.GetMedialSurface()->EvaluateDerivative(u, v, 1, 0, 0, 3, X10.data_block());
  mp.GetMedialSurface()->EvaluateDerivative(u, v, 0, 1, 0, 3, X01.data_block());
  mp.GetMedialSurface()->EvaluateDerivative(u, v, 2, 0, 0, 3, X20.data_block());
  mp.GetMedialSurface()->EvaluateDerivative(u, v, 1, 1, 0, 3, X11.data_block());
  mp.GetMedialSurface()->EvaluateDerivative(u, v, 0, 2, 0, 3, X02.data_block());

  // Compute the differential geometry
  GeometryDescriptor gd;
  gd.SetTwoJet(
    X00.data_block(), X10.data_block(), X01.data_block(),
    X20.data_block(), X11.data_block(), X02.data_block());

  // Check the symbols
  double E = dot_product(X10, X10);
  double F = dot_product(X10, X01);
  double G = dot_product(X01, X01);
  double Eu = 2.0 * dot_product(X10, X20);
  double Ev = 2.0 * dot_product(X10, X11);
  double Gu = 2.0 * dot_product(X01, X11);
  double Gv = 2.0 * dot_product(X01, X02);
  double Fu = dot_product(X20, X01) + dot_product(X11, X10);
  double Fv = dot_product(X11, X01) + dot_product(X02, X10);

  double g = E * G - F * F;

  double G111 = (G * Eu - 2 * F * Fu + F * Ev) / (2 * g);
  double G121 = (G * Ev - F * Gu) / (2 * g);
  double G221 = (2 * G * Fv - G * Gu - F * Gv) / (2 * g);
  double G112 = (2 * E * Fu - E * Ev - F * Eu) / (2 * g);
  double G122 = (E * Gu - F * Ev) / (2 * g);
  double G222 = (E * Gv - 2 * F * Fv + F * Gu) / (2 * g);

  cout << "E = " << E << " vs " << gd.xCovariantTensor[0][0] << endl;
  cout << "F = " << F << " vs " << gd.xCovariantTensor[0][1] << endl;
  cout << "G = " << G << " vs " << gd.xCovariantTensor[1][1] << endl;
  cout << "g = " << g << " vs " << gd.g << endl;

  cout << "G 11 1 = " << G111 << " vs " << gd.xChristoffelSecond[0][0][0] << endl;
  cout << "G 12 1 = " << G121 << " vs " << gd.xChristoffelSecond[0][1][0] << endl;
  cout << "G 22 1 = " << G221 << " vs " << gd.xChristoffelSecond[1][1][0] << endl;
  cout << "G 11 2 = " << G112 << " vs " << gd.xChristoffelSecond[0][0][1] << endl;
  cout << "G 12 2 = " << G122 << " vs " << gd.xChristoffelSecond[0][1][1] << endl;
  cout << "G 22 2 = " << G222 << " vs " << gd.xChristoffelSecond[1][1][1] << endl;

  gd.PrintSelf(cout);
  
  // Evaluate the nearby points
  /*
  SMLVec3d YPZ, YPP, YZP, YMP, YMZ, YMM, YZM, YPM;
  mp.GetSurface()->EvaluateDerivative(u + eps, v      , 0, 0, 0, 3, YPZ.data_block()); 
  mp.GetSurface()->EvaluateDerivative(u + eps, v + eps, 0, 0, 0, 3, YPP.data_block()); 
  mp.GetSurface()->EvaluateDerivative(u      , v + eps, 0, 0, 0, 3, YZP.data_block()); 
  mp.GetSurface()->EvaluateDerivative(u - eps, v + eps, 0, 0, 0, 3, YMP.data_block()); 
  mp.GetSurface()->EvaluateDerivative(u - eps, v      , 0, 0, 0, 3, YMZ.data_block()); 
  mp.GetSurface()->EvaluateDerivative(u - eps, v - eps, 0, 0, 0, 3, YMM.data_block()); 
  mp.GetSurface()->EvaluateDerivative(u      , v - eps, 0, 0, 0, 3, YZM.data_block()); 
  mp.GetSurface()->EvaluateDerivative(u + eps, v - eps, 0, 0, 0, 3, YPM.data_block()); 

  // Compute the partial derivatives
  SMLVec3d Y00, Y01, Y10, Y20, Y11, Y02;
  Y10 = (YPZ - YMZ) / (2.0 * eps);
  Y01 = (YZP - YZM) / (2.0 * eps);
  Y20 = (YPZ + YMZ - X00) / (eps * eps);
  Y02 = (YZP + YZM - X00) / (eps * eps);
  Y11 = (YPP + YMM - YPM - YMP) / (4 * eps * eps);

  // Compute the partial derivative tensor
  */
  return 0;
}

/**
 * This test routine checks the accuracy of basis function variation
 * computations, under various masks
 */
int TestBasisFunctionVariation(const char *fnMPDE)
{
  int iReturn = 0;

  // Load a medial PDE for the test
  CartesianMPDE mp(fnMPDE);

  // Get the model
  GenericMedialModel *model = mp.GetMedialModel();

  // Set up a selection mask
  vnl_vector<size_t> xMask(model->GetNumberOfCoefficients(), 0);
  vnl_random rnd;
  for(size_t i = 0; i < xMask.size(); i+=4)
    xMask.update(vnl_vector<size_t>(4, rnd.lrand32(0, 1)), i);

  // Set up an array of coefficient masks
  vector<CoefficientMapping *> vm;
  vm.push_back(new IdentityCoefficientMapping(model->GetNumberOfCoefficients()));
  vm.push_back(new AffineTransformCoefficientMapping(model));
  vm.push_back(new SubsetCoefficientMapping(xMask));
  
  char *nm[] = {
    "IdentityCoefficientMapping",
    "AffineTransformCoefficientMapping",
    "SubsetCoefficientMapping" };
  
  // Repeat for each mask
  for(size_t i=0; i<vm.size(); i++)
    {
    // Get the current value of the coefficients and parameters
    vnl_vector<double> C, C0 = model->GetCoefficientArray();
    vnl_vector<double> P(vm[i]->GetNumberOfParameters(), 0.0);

    // Create a random variation vector in P
    double eps = 0.0001;
    vnl_vector<double> dP(vm[i]->GetNumberOfParameters());
    for(size_t j=0; j < dP.size(); j++)
      dP[j] = rand() * 2.0 / RAND_MAX - 1.0;

    // Compute the surface derivative using each variation
    SMLVec4d f1, f2, dfNum, dfAn;

    // Evaluate the surface at finite differences
    model->SetCoefficientArray(vm[i]->Apply(C0, P + eps * dP));
    model->ComputeAtoms();
    mp.GetMedialSurface()->Evaluate(0.5, 0.5, f1.data_block());

    model->SetCoefficientArray(vm[i]->Apply(C0, P - eps * dP));
    model->ComputeAtoms();
    mp.GetMedialSurface()->Evaluate(0.5, 0.5, f2.data_block());

    dfNum = (f1 - f2) * 0.5 / eps;

    // Same evaluation using analytical derivative
    model->SetCoefficientArray(vm[i]->Apply(C0, P));
    model->ComputeAtoms();

    // Get the variation corresponding to dP
    IHyperSurface2D *xVariation = mp.GetMedialSurface()->GetVariationSurface(
      vm[i]->ApplyJacobianInParameters(C0, P, dP).data_block());
    xVariation->Evaluate(0.5, 0.5, dfAn.data_block());

    // Report
    double dMax = (dfAn - dfNum).inf_norm();
    if(dMax > eps)
      iReturn += -1;
    
    cout << "Testing mask " << nm[i] << endl;
    cout << "  right f.d.    : " << f1 << endl;
    cout << "  left f.d.     : " << f2 << endl;
    cout << "  central diff. : " << dfNum << endl;
    cout << "  numeric der.  : " << dfAn << endl;
    cout << "  maximum error : " << dMax << endl;
    }
 
  return iReturn;
}

/**
 * This test routine makes sure that the derivatives computed in the MedialPDE
 * solver and related classes are correct, bt comparing them to central
 * difference derivatives
 *
int TestDerivativesNoImage(const char *fnMPDE, const char *fnPCA)
{
  int iReturn = 0;

  // Create a cartesian medial model for testing
  CartesianMPDE mp(fnMPDE);

  // Extract the medial model itself
  GenericMedialModel *model = mp.GetMedialModel();

  // Test the most convoluted mask that we have
  cout << "**************************************************" << endl;
  cout << "** TESTIING Affine3DAndPCACoefficientMask       **" << endl;
  cout << "**************************************************" << endl;
 
  // We have to set the number of coefficients to match the PCA matrix
  mp.SetNumberOfCoefficients(8, 10);
  
  // Create a PCA based mask based on the test matrix
  double xAPCAPosn[] = 
    { 1.2, 0.1, 0.2,           // Affine Matrix
     -0.2, 0.9, 0.1, 
     -0.1, 0.1, 1.0, 
      5.0, 3.0,-2.0,           // Translation Vector
      0.5, 0.4, 0.3,-0.2, 0.3  // PCA offsets
    };

  // Read principal components from the file
  vnl_matrix<double> pcaMatrix;
  ReadMatrixFile(pcaMatrix, fnPCA);
  PrincipalComponents pca(pcaMatrix);
  cout << "PCA Matrix Size: " << pcaMatrix.rows() << " x " << pcaMatrix.columns() << endl;

  // Create the PCA/Affine optimizer with 5 modes
  PCAPlusAffineCoefficientMapping xPCAMask(model, &pca, 5);
  // PCACoefficientMapping xPCAMask(&pca, 5);

  // Create the starting point vector
  vnl_vector<double> pPCA(xAPCAPosn, 17);

  // Perform the test
  // iReturn += TestGradientComputation(model, &xPCAMask, pPCA, 3);
  iReturn += TestGradientComputation(model, &xPCAMask, 3);

  // We have to set the number of coefficients to match the PCA matrix
  mp.SetNumberOfCoefficients(2, 4);
  
  // Test straight-through gradient computation
  cout << "**************************************************" << endl;
  cout << "** TESTIING IdentityCoefficientMapping          **" << endl;
  cout << "**************************************************" << endl;
 
  IdentityCoefficientMapping xMapping(model);
  iReturn += TestGradientComputation(model, &xMapping, 3);

  // Test affine transform gradient computation
  cout << "**************************************************" << endl;
  cout << "** TESTIING AffineTransformCoefficientMapping   **" << endl;
  cout << "**************************************************" << endl;

  // Don't start at the default position, or we may get a false positive
  double xAffinePosn[] = 
    { 1.2, 0.2, 0.1, 
      0.1, 0.9,-0.2,
      0.3,-0.2, 1.4,
      7.0,-4.5, 2.3 };
  
  AffineTransformCoefficientMapping xAffineMask(model);
  vnl_vector<double> pAffine(xAffinePosn, 12);
  iReturn += TestGradientComputation(mp.GetMedialModel(), &xAffineMask, pAffine, 3);

  // Return the total of the test values
  return iReturn;
}*/

/**
 * This test routine makes sure that the derivatives computed in the MedialPDE
 * solver and related classes are correct, bt comparing them to central
 * difference derivatives
 */
int TestDerivativesNoImage(const char *fnMPDE)
{
  int iReturn = 0;

  // Create a cartesian medial model for testing
  MedialPDE mp(fnMPDE);

  // Extract the medial model itself
  GenericMedialModel *model = mp.GetMedialModel();

  // Test straight-through gradient computation
  cout << "**************************************************" << endl;
  cout << "** TESTIING IdentityCoefficientMapping          **" << endl;
  cout << "**************************************************" << endl;
 
  IdentityCoefficientMapping xMapping(model);
  iReturn += TestGradientComputation(model, &xMapping, 3);

  // Test affine transform gradient computation
  cout << "**************************************************" << endl;
  cout << "** TESTIING AffineTransformCoefficientMapping   **" << endl;
  cout << "**************************************************" << endl;

  // Don't start at the default position, or we may get a false positive
  double xAffinePosn[] = 
    { 1.2, 0.2, 0.1, 
      0.1, 0.9,-0.2,
      0.3,-0.2, 1.4,
      7.0,-4.5, 2.3 };
  
  AffineTransformCoefficientMapping xAffineMask(model);
  vnl_vector<double> pAffine(xAffinePosn, 12);
  iReturn += TestGradientComputation(mp.GetMedialModel(), &xAffineMask, pAffine, 3);

  // Test the most convoluted mask that we have
  cout << "**************************************************" << endl;
  cout << "** TESTIING Affine3DAndPCACoefficientMask       **" << endl;
  cout << "**************************************************" << endl;
 
  // Create a PCA based mask based on the test matrix
  double xAPCAPosn[] = 
    { 0.2, 0.1, 0.2,           // Affine Matrix
     -0.2,-0.1, 0.1, 
     -0.1, 0.1, 0.0, 
      5.0, 3.0,-2.0,           // Translation Vector
      0.5, 0.4, 0.3,-0.2, 0.3  // PCA offsets
    };

  // Create a PCA based on a random matrix
  vnl_matrix<double> pcaMatrix(15, model->GetNumberOfCoefficients());
  vnl_random randy;
  for(size_t q = 0; q < pcaMatrix.rows(); q++) 
    for(size_t k = 0; k < pcaMatrix.columns(); k++)
      pcaMatrix[q][k] = randy.drand32(-0.01, 0.01);
  PrincipalComponents pca(pcaMatrix);

  // Create the PCA/Affine optimizer with 5 modes
  PCAPlusAffineCoefficientMapping xPCAMask(model, &pca, 5);

  // Create the starting point vector
  vnl_vector<double> pPCA(xAPCAPosn, 17);

  // Perform the test
  iReturn += TestGradientComputation(model, &xPCAMask, pPCA, 3);

  // Return the total of the test values
  return iReturn;
}

int TestGradientTiming(const char *fnMPDE)
{
  // Create and read the MPDE
  CartesianMPDE mp(2, 4, 32, 80);
  mp.LoadFromParameterFile(fnMPDE);

  // Test straight-through gradient computation
  IdentityCoefficientMapping xMask(mp.GetMedialModel());
  return TestGradientComputation(mp.GetMedialModel(), &xMask);
}

void MakeFlatTemplate(FourierSurface *xSurface)
{
  size_t n = 10;
  size_t q = (n+1) * (n+1);
  double *u = new double[q], *v = new double[q];
  double *x = new double[q], *y = new double[q];
  double *z = new double[q], *rho = new double[q];

  size_t k = 0;
  for(size_t i = 0; i <= n; i++) for(size_t j = 0; j <= n; j++)
    {
    u[k] = i * 1.0 / n;
    v[k] = j * 1.0 / n;
    x[k] = 6 * u[k];
    y[k] = 12 * v[k];
    z[k] = 0.0;
    rho[k] = -0.45;
    k++;
    }

  xSurface->FitToData(q, 0, u, v, x);
  xSurface->FitToData(q, 1, u, v, y);
  xSurface->FitToData(q, 2, u, v, z);
  xSurface->FitToData(q, 3, u, v, rho);

  delete u; delete v; delete x; delete y; delete z; delete rho;
}

int TestDerivativesWithImage(const char *fnMPDE)
{
  // Return Code
  int iReturn = 0;

  // Load a cartesian medial PDE
  MedialPDE mp(fnMPDE);
  
  // Define a test image (this in not a real thing)
  GenericMedialModel *model = mp.GetMedialModel();
  SMLVec3d C = model->GetCenterOfRotation();
  TestFloatImage img( C, 7.0, 4.0 );

  // Define a test point set and a test radius function for computation
  vnl_vector<double> px(model->GetNumberOfAtoms(), 0.0);
  vnl_vector<double> py(model->GetNumberOfAtoms(), 0.0);
  vnl_vector<double> pz(model->GetNumberOfAtoms(), 0.0);
  vnl_vector<double> rad(model->GetNumberOfAtoms(), 0.0);
  vnl_random randy;
  for(size_t i = 0; i < rad.size(); i++)
    {
    px[i] = model->GetAtomArray()[i].X[0] + randy.drand32(-0.1, 0.1);
    py[i] = model->GetAtomArray()[i].X[1] + randy.drand32(-0.1, 0.1);
    pz[i] = model->GetAtomArray()[i].X[2] + randy.drand32(-0.1, 0.1);
    rad[i] = model->GetAtomArray()[i].R * randy.drand32(0.9, 1.1);
    }

  // Create an array of image match terms
  vector<EnergyTerm *> vt;
  vt.push_back(new MedialBendingEnergyTerm(model));
  vt.push_back(new DistanceToPointSetEnergyTerm(
      model, px.data_block(), py.data_block(), pz.data_block()));
  vt.push_back(new BoundaryImageMatchTerm(&img));
  vt.push_back(new DistanceToRadiusFieldEnergyTerm(model, rad.data_block()));
  vt.push_back(new ProbabilisticEnergyTerm(&img, 4));
  vt.push_back(new MedialAnglesPenaltyTerm(model));
  vt.push_back(new MedialRegularityTerm(model));
  vt.push_back(new BoundaryJacobianEnergyTerm());
  vt.push_back(new BoundaryGradRPenaltyTerm());

  // Create an array of masks
  vector<CoefficientMapping *> vm;
  vm.push_back(new SubsetCoefficientMapping(model->GetRadialCoefficientMask()));
  vm.push_back(new IdentityCoefficientMapping(model));
  vm.push_back(new AffineTransformCoefficientMapping(model));

  // Create labels
  char *nt[] = {
    "MedialBendingEnergyTerm",
    "DistanceToPointSetEnergyTerm",
    "BoundaryImageMatchTerm",
    "DistanceToRadiusFieldEnergyTerm",
    "ProbabilisticEnergyTerm",
    "MedialAnglesPenaltyTerm",
    "MedialRegularityTerm",
    "BoundaryJacobianEnergyTerm", 
    "BoundaryGradRPenaltyTerm" };
  char *nm[] = {
    "RadialOnlyCoefficientMapping",
    "IdentityCoefficientMapping",
    "AffineTransformCoefficientMapping" };

  double stepsize[] = {0.1, 0.1, 0.1, 0.01, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

  // Loop over both options
  size_t i, j;
  for(i = 0; i < vt.size(); i++) for(j = 0; j < vm.size(); j++)
    {
    cout << "-------------------------------------------------------------------" << endl;
    cout << " TESTING " << nt[i] << " ON " << nm[j] << endl;
    cout << "-------------------------------------------------------------------" << endl;
    
    // Set up the test
    MedialOptimizationProblem mop(model, vm[j]);
    mop.AddEnergyTerm(vt[i], 1.0);
    iReturn += TestOptimizerGradientComputation(mop, *vm[j], model, stepsize[i]);
    }

  // Delete both pointers
  for(i = 0; i < vt.size(); i++) 
    delete vt[i];
  
  for(j = 0; j < vm.size(); j++)
    delete vm[j];

  return iReturn;
}

int TestAffineTransform(const char *fnMPDE)
{
  // Load the model from file
  MedialPDE mp(fnMPDE);
  GenericMedialModel *model = mp.GetMedialModel();
  vnl_vector<double> hint = model->GetHintArray();

  // Generate a random affine transform. Because the PDE solution may change
  // under shear and non-uniform scaling, we stick just to the similarity
  // transform (A is an orthogonal matrix)
  AffineTransformDescriptor::Vec b(3, 0.0), c(3, 0.0);
  AffineTransformDescriptor::Mat A(3, 3, 0.0);

  // Generate random vectors V1, V2 for making the matrix A, plus b and c
  SMLVec3d v1, v2;
  vnl_random rnd;
  for(size_t i = 0; i < 3; i++)
    {
    b[i] = rnd.drand32(-1.0, 1.0);
    c[i] = rnd.drand32(-1.0, 1.0);
    v1[i] = rnd.drand32(-1.0, 1.0);
    v2[i] = rnd.drand32(-1.0, 1.0);
    }

  // Create the orthogonal matrix A
  A.set_column(0, v1.normalize());
  A.set_column(1, vnl_cross_3d(A.get_column(0), v2).normalize());
  A.set_column(2, vnl_cross_3d(A.get_column(0), A.get_column(1)).normalize());
  double s = rnd.drand32(0.0, 1.0); A *= s;

  // Print the transform
  cout << "Affine Matrix: " << endl << A << endl;
  cout << "A.t() * A: " << endl << A.transpose() * A << endl;
  cout << "Translation: " << b << endl;
  cout << "Center of Rotation: " << c << endl;

  // Get a list of interior points
  vector<SMLVec3d> X1;

  MedialInternalPointIterator it = model->GetInternalPointIterator(6);
  for( ; !it.IsAtEnd(); ++it)
    {
    SMLVec3d x = GetInternalPoint(it, model->GetAtomArray());
    SMLVec3d Ax = (A * (x - c)) + b + c;
    X1.push_back(Ax); 
    }

  cout << "Obtained initial point data" << endl;

  // Get the coefficients of the model
  vnl_vector<double> xCoeff = model->GetCoefficientArray();

  // Scale the hint by the square of s (hint = phi)
  hint *= s * s;

  // Apply the transform to the model
  const AffineTransformDescriptor *adt = model->GetAffineTransformDescriptor();
  model->SetCoefficientArray(adt->ApplyAffineTransform(xCoeff, A, b, c));
  model->ComputeAtoms(hint.data_block());

  // Max difference accumulator
  double maxdiff = 0.0;

  size_t q = 0;
  for(it = model->GetInternalPointIterator(6); !it.IsAtEnd(); ++it, ++q)
    {
    SMLVec3d x = GetInternalPoint(it, model->GetAtomArray());
    double diff = (x - X1[q]).two_norm();
    maxdiff = std::max(maxdiff, diff);
    }

  // Report
  cout << "Maximum difference in Affine Transform computations: " << maxdiff << endl;
  return (maxdiff > 1.0e-6);
}


int usage()
{
  cout << "testpde: MedialPDE Test Module" << endl;
  cout << "  usage: testpde TEST_ID [parameters] " << endl;
  cout << "  tests: " << endl;
  cout << "    DERIV1 XX.mpde XX.mat      Check analytic derivative PDEs." << endl;
  cout << "    DERIV2 XX.mpde             Check gradient computation in image match terms." << endl;
  cout << "    DERIV3 XX.mpde             Check variations on basis functions." << endl;
  cout << "    DERIV4 XX.mpde             Test diff. geom. operators." << endl;
  cout << "    AFFINE XX.mpde             Test affine transform computation." << endl;
  cout << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  // Different tests that can be executed
  if(argc == 1) return usage();

  // Choose a test depending on the parameters
  if(0 == strcmp(argv[1], "DERIV1") && argc > 2)
    return TestDerivativesNoImage(argv[2]);
  else if(0 == strcmp(argv[1], "DERIV2") && argc > 2)
    return TestDerivativesWithImage(argv[2]);
  else if(0 == strcmp(argv[1], "DERIV3") && argc > 2)
    return TestBasisFunctionVariation(argv[2]);
  else if(0 == strcmp(argv[1], "DERIV4") && argc > 2)
    return TestDifferentialGeometry(argv[2]);
  else if(0 == strcmp(argv[1], "DERIV5") && argc > 2)
    return TestGradientTiming(argv[2]);
  else if(0 == strcmp(argv[1], "WEDGE"))
    return TestWedgeVolume();
  else if(0 == strcmp(argv[1], "VOLUME1") && argc > 2)
    return TestVolumeComputation(argv[2]);
  else if(0 == strcmp(argv[1], "AFFINE") && argc > 2)
    return TestAffineTransform(argv[2]);
  else 
    return usage();
}
