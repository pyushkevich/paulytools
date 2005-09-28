#include "ScriptInterface.h"
#include "MedialPDESolver.h"
#include "BasisFunctions2D.h"
#include "ITKImageWrapper.h"
#include "MedialPDERenderer.h"
#include "OptimizationTerms.h"
#include "Procrustes.h"
#include "Registry.h"
#include <iostream>
#include <ctime>

using namespace std;

class CodeTimer 
{
public:
  CodeTimer() 
    { tElapsed = 0.0; }

  void Start()
    { tStart = clock(); }
  
  void Stop()
    { tElapsed += (clock() - tStart) * 1.0 / CLOCKS_PER_SEC; }

  void Reset()
    { tElapsed = 0.0; }

  double Read()
    { return tElapsed; }

private:
  clock_t tStart;
  double tElapsed;
};

void TestGradientComputation(
  MedialPDESolver *xSolver, IMedialCoefficientMask *xMask)
{
  // Timers for different routines
  CodeTimer tCentral, tAnalytic;
  
  cout << "*****************************************" << endl;
  cout << "Comparing Numeric Derivatives to Analytic" << endl;
  cout << "*****************************************" << endl;

  // Create the first problem
  MedialPDESolver S1(xSolver->GetGridU(), xSolver->GetGridV());
  S1.SetMedialSurface(xSolver->GetMedialSurface());

  // Create the second problem
  MedialPDESolver S2(xSolver->GetGridU(), xSolver->GetGridV());
  S2.SetMedialSurface(xSolver->GetMedialSurface());

  for(size_t iComp = 0; iComp < xMask->GetNumberOfCoefficients(); iComp++)
    {
    cout << "----------------------------------------" << endl;
    cout << " Coefficient " << iComp << endl;
    cout << "----------------------------------------" << endl;
    
    // Get an adapter corresponding to the i-th component of the surface
    IHyperSurface2D *xEta = xMask->GetComponentSurface(iComp);

    // Create an array of atoms to hold the derivative
    MedialAtom *dAtoms = 
      new MedialAtom[xSolver->GetAtomGrid()->GetNumberOfAtoms()];

    // Solve for Phi at the current location
    xSolver->Solve();

    // Compute the variational derivative
    tAnalytic.Start();
    xSolver->ComputeVariationalDerivative(xEta, dAtoms);
    tAnalytic.Stop();

    // Release the variation
    xMask->ReleaseComponentSurface(xEta);

    // Now, use central differences to compute the same
    double cc = xMask->GetCoefficient(iComp), eps = 0.0001;

    // Solve both problems
    tCentral.Start();
    xMask->SetCoefficient(iComp, cc + eps);
    S1.Solve();

    xMask->SetCoefficient(iComp, cc - eps);
    S2.Solve();
    tCentral.Stop();

    // Restore the coefficient to its correct value
    xMask->SetCoefficient(iComp, cc);

    // Get the atom arrays
    MedialAtom *A1 = S1.GetAtomArray();
    MedialAtom *A2 = S2.GetAtomArray();

    // Compute the values of dPhi/dC between central difference and analytic
    // derivatives
    unsigned int n = xSolver->GetAtomGrid()->GetNumberOfAtoms();
    double xMaxPhiDiff = 0.0;
    double xMaxBndDiff = 0.0;
    double xMaxBndNrmDiff = 0.0;
    double xMaxGRMS = 0.0;

    for(unsigned int i = 0; i < n; i++)
      {
      // Point to the atoms
      MedialAtom &a1 = S1.GetAtomArray()[i];
      MedialAtom &a2 = S2.GetAtomArray()[i];
      MedialAtom &a0 = dAtoms[i];

      // Take the largest difference in Phi
      double dcPhi = 0.5 * (a1.F - a2.F) / eps;
      double daPhi = a0.F;
      if(fabs(dcPhi - daPhi) > xMaxPhiDiff)
        xMaxPhiDiff = fabs(dcPhi - daPhi);

      // Also look at the differences in grad R mag sqr.
      double dcGRMS = 0.5 * (a1.xGradRMagSqr - a2.xGradRMagSqr) / eps;
      double daGRMS = a0.xGradRMagSqr;
      if(fabs(dcGRMS - daGRMS) > xMaxGRMS)
        xMaxGRMS = fabs(dcGRMS - daGRMS);

      // Take the largest difference in boundary derivs
      for(unsigned int k = 0; k < 2; k++)
        {
        SMLVec3d dcBnd = (0.5 / eps) * (a1.xBnd[k].X - a2.xBnd[k].X);
        SMLVec3d daBnd = a0.xBnd[k].X;
        double xBndDiff = (dcBnd - daBnd).two_norm();
        if(xBndDiff > xMaxBndDiff) xMaxBndDiff = xBndDiff;

        SMLVec3d dcBndNrm = (0.5 / eps) * (a1.xBnd[k].N - a2.xBnd[k].N);
        SMLVec3d daBndNrm = a0.xBnd[k].N;
        double xBndNrmDiff = (dcBndNrm - daBndNrm).two_norm();
        if(xBndNrmDiff > xMaxBndNrmDiff) xMaxBndNrmDiff = xBndNrmDiff;
        }
      }

    // Print out the max difference
    cout << "Maximal differences between numeric and analytic derivatives:" << endl;
    cout << "Phi                      : " <<  xMaxPhiDiff << endl;
    cout << "Sqr. Grad. Mag. of Phi   : " <<  xMaxGRMS << endl;
    cout << "Boundary Node Positions  : " <<  xMaxBndDiff << endl;
    cout << "Boundary Node Normals    : " <<  xMaxBndNrmDiff << endl;

    // Now, test the effect on SolutionData objects
    SolutionData sd0(xSolver), sd1(&S1), sd2(&S2);
    PartialDerivativeSolutionData pdsd(&sd0, dAtoms);

    // Compute the boundary weights derivatives
    sd0.UpdateBoundaryWeights();
    sd1.UpdateBoundaryWeights();
    sd2.UpdateBoundaryWeights();
    pdsd.UpdateBoundaryWeights();

    // Compute the area weight difference
    double xMaxBndAreaDiff = 0.0;
    MedialBoundaryPointIterator *itBnd 
      = xSolver->GetAtomGrid()->NewBoundaryPointIterator();
    while(!itBnd->IsAtEnd())
      {
      double aNumeric = 0.5 * ( 
        sd1.xBoundaryWeights[itBnd->GetIndex()] - 
        sd2.xBoundaryWeights[itBnd->GetIndex()]) / eps;
      double aAnalytic = pdsd.xBoundaryWeights[itBnd->GetIndex()];
      double d = fabs(aAnalytic - aNumeric);
      if(d > eps)
        cout << itBnd->GetAtomIndex() << " : " << d << " " << aNumeric << " " << aAnalytic << endl;
      if(d > xMaxBndAreaDiff) xMaxBndAreaDiff = d;
      ++(*itBnd); 
      }
    delete itBnd;

    cout << "Max difference between numeric and analytic derivs (Bnd Area) : " 
      <<  xMaxBndAreaDiff << endl;  

    cout << "Difference in total area derivatives : " 
      << fabs(0.5 * (sd1.xBoundaryArea - sd2.xBoundaryArea) / eps - 
        pdsd.xBoundaryArea) << endl;

    // Compute the difference in internal point volume weights
    size_t nInternal = 4;
    sd0.UpdateInternalWeights(nInternal);
    sd1.UpdateInternalWeights(nInternal);
    sd2.UpdateInternalWeights(nInternal);
    pdsd.UpdateInternalWeights(nInternal);

    // Compute the area weight difference
    double xMaxCellVolumeDiff = 0.0;
    MedialInternalPointIterator *itPoint 
      = xSolver->GetAtomGrid()->NewInternalPointIterator(nInternal);
    while(!itPoint->IsAtEnd())
      {
      double aNumeric = 0.5 * ( 
        sd1.xInternalWeights[itPoint->GetIndex()] - 
        sd2.xInternalWeights[itPoint->GetIndex()]) / eps;
      double aAnalytic = pdsd.xInternalWeights[itPoint->GetIndex()];
      SMLVec3d vNumeric = 0.5 * (
        sd1.xInternalPoints[itPoint->GetIndex()] -
        sd2.xInternalPoints[itPoint->GetIndex()]) / eps;
      SMLVec3d vAnalytic = pdsd.xInternalPoints[itPoint->GetIndex()];
      double d = fabs(aAnalytic - aNumeric);
      // double d = (vAnalytic - vNumeric).two_norm();
      if(d > xMaxCellVolumeDiff) xMaxCellVolumeDiff = d;
      ++(*itPoint); 
      }
    delete itPoint;

    cout << "Max difference between numeric and analytic derivs (Volume Elt) : " 
      <<  xMaxCellVolumeDiff << endl;  

    cout << "Difference in total volumes : " 
      << fabs(0.5 * (sd1.xInternalVolume - sd2.xInternalVolume) / eps - 
        pdsd.xInternalVolume) << endl;
    }

  cout << "========================================" << endl;
  cout << "Central difference computed in " << tCentral.Read() 
    << " sec., Analytic in " << tAnalytic.Read() << " sec." << endl; 
  cout << "========================================" << endl;
}

void TestOptimizerGradientComputation(
  MedialOptimizationProblem &mop, 
  IMedialCoefficientMask &xMask,
  MedialPDESolver *xSolver)
{
  // Set up timers
  CodeTimer tAnalytic, tNumeric;

  // Compute the analytic gradient
  vnl_vector<double> xGradient(xMask.GetNumberOfCoefficients());
  tAnalytic.Start();
  mop.ComputeGradient(xMask.GetCoefficientArray(), xGradient.data_block());
  tAnalytic.Stop();
  
  // Speed up this computation
  xSolver->SetSolutionAsInitialGuess();
  
  // Compute the numeric gradient of the objective function
  double eps = 0.01;
  for(size_t i = 0; i < xMask.GetNumberOfCoefficients(); i++)
    {
    // Get the current value of the coefficient
    double ci = xMask.GetCoefficient(i);

    // Compute the objective function at point 1
    xMask.SetCoefficient(i, ci + eps);
    tNumeric.Start();
    double f1 = mop.Evaluate(xMask.GetCoefficientArray());
    tNumeric.Stop();

    // Compute the objective at point 2
    xMask.SetCoefficient(i, ci - eps);
    tNumeric.Start();
    double f2 = mop.Evaluate(xMask.GetCoefficientArray());
    tNumeric.Stop();

    // Restore the coefficient
    xMask.SetCoefficient(i, ci);
    
    // Print the derivative
    double dNumeric = 0.5 * (f1 - f2) / eps;
    cout << "Derivatives wrt. " << i << ". Numeric = " << dNumeric 
      << ";   Analytic = " << xGradient[i] 
      << ";   Difference = " << fabs(xGradient[i] - dNumeric) << endl;
    }

  // Report difference in time
  cout << "Central difference computed in " << tNumeric.Read() 
    << " sec., Analytic in " << tAnalytic.Read() << " sec." << endl; 
}
