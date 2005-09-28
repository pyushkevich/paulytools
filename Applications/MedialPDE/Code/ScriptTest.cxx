#include "ScriptInterface.h"
#include "BasisFunctions2D.h"
#include "MedialAtom.h"
#include "MedialPDESolver.h"
#include "CartesianMedialAtomGrid.h"
#include "OptimizationTerms.h"
#include "CoefficientMask.h"
#include "TestSolver.h"
#include "vnl/vnl_erf.h"

#include <string>
#include <iostream>

using namespace std;
using namespace medialpde;

string dirWork = "/home/pauly/data2005/Stanley/data/";
//string dirWork = "/mnt/data2/PUBLIC/Data/Input/StanleySchizophrenia/";

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
void TestAreaAndVolume(MedialPDESolver *xSolver)
{
  unsigned int nCuts = 6;

  // Compute the area 
  MedialAtomGrid *xGrid = xSolver->GetAtomGrid();
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
  MedialBoundaryPointIterator *it = xGrid->NewBoundaryPointIterator();
  while(!it->IsAtEnd())
    {
    // Evaluate the image at this location
    BoundaryAtom &B = GetBoundaryPoint(it, xSolver->GetAtomArray()); 
    xVol += dot_product(B.X , B.N) * xAreaWeights[it->GetIndex()] / 3.0; 

    // On to the next point
    ++(*it); 
    }

  // Return match scaled by total weight
  cout << "Verification " << xVol << endl;

  // One more test
  xVol = 0.0;
  MedialProfileIntervalIterator *itProf = xGrid->NewProfileIntervalIterator(nCuts);
  for(; !itProf->IsAtEnd(); ++(*itProf))
    xVol += xProfileWeights[itProf->GetIndex()];
  cout << "Profiles add up to " << xVol << endl;

  // Clean up
  delete it;
  
  // Check the Jacobian
  SolutionData S(xSolver);
  BoundaryJacobianEnergyTerm termJac;
  termJac.ComputeEnergy(&S);
  termJac.PrintReport(cout);
}

void Test01()
{
  MedialPDE *mp = new MedialPDE(3, 5, 20, 40);
  //mp->LoadFromParameterFile(fMrep.c_str());
  
  mp->LoadFromDiscreteMRep("/tmp/surf01.txt",-0.3);
  // mp->GenerateSampleModel();
  mp->Solve();
  mp->SaveBYUMesh("temp.byu");

  // Make sure areas and volumes add up
  TestAreaAndVolume(mp->GetSolver());
  
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
  MedialPDE *mp = new MedialPDE(5, 5, 50, 100);
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
  mp->SetOptimizerToGradientDescent(0.01);
  mp->SetOptimizationToAffine();
  mp->RunOptimization(img, 30);
  mp->SaveToParameterFile((dirWork + "avg/average_mrepL_02.mpde").c_str());
  mp->SaveBYUMesh((dirWork + "avg/average_mrepL_02.byu").c_str());
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

/**
 * This test routine makes sure that the derivatives computed in the MedialPDE
 * solver and related classes are correct, bt comparing them to central
 * difference derivatives
 */
int TestDerivativesNoImage(const char *fnMPDE)
{
  MedialPDE mp(2, 4, 33, 65, 0.5, 0, 0);
  mp.LoadFromParameterFile(fnMPDE);

  // Test straight-through gradient computation
  PassThroughCoefficientMask xMask(mp.GetSurface());
  TestGradientComputation(mp.GetSolver(), &xMask);

  // Test affine transform gradient computation
  AffineTransformCoefficientMask xAffineMask(mp.GetSurface());
  TestGradientComputation(mp.GetSolver(), &xAffineMask);

  return 0;
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
  // Load the Medial PDE
  MedialPDE mp(2, 4, 33, 81, 0.5, 0, 0);
  mp.LoadFromParameterFile(fnMPDE);
  
  // Define a test image (this in not a real thing)
  SMLVec3d C = mp.GetSurface()->GetCenterOfRotation().extract(3);
  TestFloatImage img( C, 7.0, 4.0 );

  // Define another solver to represent the template - this is used in 
  // conjunction with the regularity prior 
  MedialPDE mpTemplate(4, 4, 33, 81);
  MakeFlatTemplate(mpTemplate.GetSurface());
  mpTemplate.Solve();

  // Create an array of image match terms
  vector<EnergyTerm *> vt;
  vt.push_back(new MedialRegularityTerm(
      mpTemplate.GetSolver()->GetAtomArray(), mpTemplate.GetSolver()->GetAtomGrid()));
  vt.push_back(new BoundaryJacobianEnergyTerm());
  vt.push_back(new BoundaryImageMatchTerm(&img));
  vt.push_back(new ProbabilisticEnergyTerm(&img, 4));

  // Create an array of masks
  vector<IMedialCoefficientMask *> vm;
  vm.push_back(new PassThroughCoefficientMask(mp.GetSurface()));
  vm.push_back(new AffineTransformCoefficientMask(mp.GetSurface()));

  // Create labels
  char *nt[] = {
    "MedialRegularityTerm",
    "BoundaryJacobianEnergyTerm", 
    "BoundaryImageMatchTerm",
    "ProbabilisticEnergyTerm" };
  char *nm[] = {
    "PassThroughCoefficientMask",
    "AffineTransformCoefficientMask" };

  // Loop over both options
  size_t i, j;
  for(i = 0; i < vt.size(); i++) for(j = 0; j < vm.size(); j++)
    {
    cout << "-------------------------------------------------------------------" << endl;
    cout << " TESTING " << nt[i] << " ON " << nm[j] << endl;
    cout << "-------------------------------------------------------------------" << endl;
    
    // Set up the test
    MedialOptimizationProblem mop(mp.GetSolver(), vm[j]);
    mop.AddEnergyTerm(vt[i], 1.0);
    TestOptimizerGradientComputation(mop, *vm[j], mp.GetSolver());
    }

  // Delete both pointers
  for(i = 0; i < vt.size(); i++) 
    delete vt[i];
  
  for(j = 0; j < vm.size(); j++)
    delete vm[j];

  return 0;
}


void TestFDMasksInMedialSolver()
{
  // Load a medial PDE
  MedialPDE mp(8, 12, 33, 17);
  
  
  // Define a function

}

int usage()
{
  cout << "testpde: MedialPDE Test Module" << endl;
  cout << "  usage: testpde TEST_ID [parameters] " << endl;
  cout << "  tests: " << endl;
  cout << "    DERIV1 XX.mpde    Check analytic derivative PDEs." << endl;
  cout << "    DERIV2 XX.mpde    Check gradient computation in image match terms." << endl;
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
  else 
    return usage();
}
