#include "ScriptInterface.h"
#include "FourierSurface.h"
#include "BasisFunctions2D.h"
#include "MedialAtom.h"
#include "MedialPDESolver.h"
#include "CartesianMedialAtomGrid.h"
#include "OptimizationTerms.h"

#include <string>
#include <iostream>

using namespace std;
using namespace medialpde;

string dirWork = "/home/pauly/data2005/Stanley/data/";

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
  SMLVec3d *xInternal = new SMLVec3d[xGrid->GetNumberOfInternalPoints(nCuts)];
  ComputeMedialInternalPoints(
    xGrid, xSolver->GetAtomArray(), nCuts, xInternal);
  
  xVol = ComputeMedialInternalVolumeWeights(
    xGrid, xInternal, nCuts, xVolWeights);

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

  // Clean up
  delete it;
  
  // Check the Jacobian
  SolutionData S(xSolver, false);
  BoundaryJacobianEnergyTerm termJac;
  termJac.ComputeEnergy(&S);
  termJac.PrintReport(cout);

  // Return match scaled by total weight
  cout << "Verification " << xVol << endl;
}

void Test01()
{
  // string fMrep = dirWork + "avg/average_mrepL_01.mpde";
  
  
  MedialPDE *mp = new MedialPDE(3, 5, 24);
  //mp->LoadFromParameterFile(fMrep.c_str());
  
  mp->LoadFromDiscreteMRep("/tmp/surf01.txt",-0.5,16);
  mp->Solve();
  mp->SaveBYUMesh("temp.byu");

  // Make sure areas and volumes add up
  TestAreaAndVolume(mp->GetSolver());
  
  RenderMedialPDE(mp);
}

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

void Test03()
{
  // Create a medial PDE object
  MedialPDE *mp = new MedialPDE(5, 5, 10);
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

int main(int argc, char *argv[])
{
  // TestCellVolume();
  // TestCartesianGrid();
  Test01();
  
  return 0;
}
