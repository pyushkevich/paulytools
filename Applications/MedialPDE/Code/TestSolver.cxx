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

// An inline function to update max counters
inline void UpdateMax(double xValue, double &xMaxValue)
{
  if(xValue > xMaxValue)
    xMaxValue = xValue;
}

void TestGradientComputationPrintReport(
  double xMaxPhiDiff, double xMaxGRMS, double xMaxBndDiff, 
  double xMaxBndNrmDiff, double xMaxBndAreaDiff, double xMaxCellVolumeDiff,
  double xTotalAreaDiff, double xTotalVolumeDiff)
{
    // Print out the max difference
    cout << "Maximal differences between numeric and analytic derivatives:" << endl;
    cout << "Phi                      : " <<  xMaxPhiDiff << endl;
    cout << "Sqr. Grad. Mag. of Phi   : " <<  xMaxGRMS << endl;
    cout << "Boundary Node Positions  : " <<  xMaxBndDiff << endl;
    cout << "Boundary Node Normals    : " <<  xMaxBndNrmDiff << endl;
    cout << "Boundary Area Elements   : " <<  xMaxBndAreaDiff << endl;
    cout << "Total Boundary Area      : " <<  xMaxBndAreaDiff << endl;
    cout << "Internal Volume Elements : " <<  xTotalAreaDiff << endl;
    cout << "Total Internal Volume    : " <<  xTotalVolumeDiff << endl;
}

int TestGradientComputation(
  MedialPDESolver *xSolver, IMedialCoefficientMask *xMask)
{
  // Timers for different routines
  CodeTimer tCentral, tAnalytic;
  
  // The epsilon for central difference tests
  double eps = 0.0001;
  
  // Max values for different classes of error
  double xGlobalMaxPhiDiff = 0.0, xGlobalMaxGRMS = 0.0;
  double xGlobalMaxBndDiff = 0.0, xGlobalMaxBndNrmDiff = 0.0;
  double xGlobalMaxBndAreaDiff = 0.0, xGlobalMaxCellVolumeDiff = 0.0;
  double xGlobalTotalAreaDiff = 0.0, xGlobalTotalVolumeDiff = 0.0;

  cout << "*****************************************" << endl;
  cout << "Comparing Numeric Derivatives to Analytic" << endl;
  cout << "*****************************************" << endl;

  // Create the first problem
  MedialPDESolver S1(xSolver->GetGridU(), xSolver->GetGridV());
  S1.SetMedialSurface(xSolver->GetMedialSurface());

  // Create the second problem
  MedialPDESolver S2(xSolver->GetGridU(), xSolver->GetGridV());
  S2.SetMedialSurface(xSolver->GetMedialSurface());

  // Create a list of variational surfaces: some components, some random
  // variations
  vector<IHyperSurface2D *> xVariations;
  vector<string> xVariationNames;
  vector< vnl_vector<double> > xVariationCoeffs;

  // Get the current coefficient vector
  vnl_vector<double> x0(
    xMask->GetCoefficientArray(), xMask->GetNumberOfCoefficients());

  // Create the per-component variations
  size_t iComp;
  for( iComp = 0; iComp < xMask->GetNumberOfCoefficients(); iComp++)
    {
    ostringstream oss;
    oss << " Coefficient " << iComp;
    
    vnl_vector<double> xCoeff(xMask->GetNumberOfCoefficients(), 0.0);
    xCoeff[iComp] = 1.0;
    
    xVariations.push_back(xMask->GetComponentSurface(iComp));
    xVariationNames.push_back(oss.str());
    xVariationCoeffs.push_back(xCoeff);
    }

  // Create mixed variations
  for(size_t i = 0; i < 3; i++)
    {    
    vnl_vector<double> dx(xMask->GetNumberOfCoefficients());
    for(size_t j=0; j<dx.size(); j++)
      dx[j] = rand() * 2.0 / RAND_MAX - 1.0;

    xVariations.push_back(xMask->GetVariationSurface(dx.data_block()));
    xVariationNames.push_back(" Random Variation ");
    xVariationCoeffs.push_back(dx);
    }
  
  // Repeat for all variations
  for(size_t iVar = 0; iVar < xVariations.size(); iVar++)
    {
    cout << "----------------------------------------" << endl;
    cout << xVariationNames[iVar] << endl;
    cout << "----------------------------------------" << endl;
    
    // Get an adapter corresponding to the i-th component of the surface
    IHyperSurface2D *xEta = xVariations[iVar];

    // Create an array of atoms to hold the derivative
    MedialAtom *dAtoms = 
      new MedialAtom[xSolver->GetAtomGrid()->GetNumberOfAtoms()];

    // Solve for Phi at the current location
    xSolver->Solve();

    // Compute the variational derivative
    tAnalytic.Start();
    xSolver->ComputeVariationalDerivative(xEta, dAtoms);
    tAnalytic.Stop();

    // Solve both problems
    tCentral.Start();
    xMask->SetCoefficientArray( (x0 + eps * xVariationCoeffs[iVar]).data_block() );
    S1.Solve();

    xMask->SetCoefficientArray( (x0 - eps * xVariationCoeffs[iVar]).data_block() );
    S2.Solve();
    tCentral.Stop();

    // Restore the coefficient to its correct value
    xMask->SetCoefficientArray( x0.data_block() );

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
      UpdateMax(fabs(dcPhi - daPhi), xMaxPhiDiff);

      // Also look at the differences in grad R mag sqr.
      double dcGRMS = 0.5 * (a1.xGradRMagSqr - a2.xGradRMagSqr) / eps;
      double daGRMS = a0.xGradRMagSqr;
      UpdateMax(fabs(dcGRMS - daGRMS), xMaxGRMS);

      // Take the largest difference in boundary derivs
      for(unsigned int k = 0; k < 2; k++)
        {
        SMLVec3d dcBnd = (0.5 / eps) * (a1.xBnd[k].X - a2.xBnd[k].X);
        SMLVec3d daBnd = a0.xBnd[k].X;
        UpdateMax((dcBnd - daBnd).two_norm(), xMaxBndDiff);

        SMLVec3d dcBndNrm = (0.5 / eps) * (a1.xBnd[k].N - a2.xBnd[k].N);
        SMLVec3d daBndNrm = a0.xBnd[k].N;
        UpdateMax((dcBndNrm - daBndNrm).two_norm(), xMaxBndNrmDiff);
        }
      }

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
      UpdateMax(fabs(aAnalytic - aNumeric), xMaxBndAreaDiff);
      ++(*itBnd); 
      }
    delete itBnd;


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
      UpdateMax(fabs(aAnalytic - aNumeric), xMaxCellVolumeDiff);

      ++(*itPoint); 
      }
    delete itPoint;

    // Compare total area and volume derivatives
    double xTotalVolumeDiff = 
      fabs(0.5 * (sd1.xInternalVolume - sd2.xInternalVolume) / eps 
        - pdsd.xInternalVolume);

    double xTotalAreaDiff = 
      fabs(0.5 * (sd1.xBoundaryArea - sd2.xBoundaryArea) / eps 
        - pdsd.xBoundaryArea);

    // Print the report
    TestGradientComputationPrintReport(
      xMaxPhiDiff, xMaxGRMS, xMaxBndDiff, 
      xMaxBndNrmDiff, xMaxBndAreaDiff, xMaxCellVolumeDiff,
      xTotalAreaDiff, xTotalVolumeDiff);

    // Update the global maximum counters
    UpdateMax(xMaxPhiDiff, xGlobalMaxPhiDiff);
    UpdateMax(xMaxGRMS, xGlobalMaxGRMS);
    UpdateMax(xMaxBndDiff, xGlobalMaxBndDiff);
    UpdateMax(xMaxBndNrmDiff, xGlobalMaxBndNrmDiff);
    UpdateMax(xMaxBndAreaDiff, xGlobalMaxBndAreaDiff);
    UpdateMax(xMaxCellVolumeDiff, xGlobalMaxCellVolumeDiff);
    UpdateMax(xTotalAreaDiff, xGlobalTotalAreaDiff);
    UpdateMax(xTotalVolumeDiff, xGlobalTotalVolumeDiff);
    }

  // Release unused variations
  for(iComp = 0; iComp < xMask->GetNumberOfCoefficients(); iComp++)
    xMask->ReleaseComponentSurface(xVariations[iComp]);
  for(iComp = 0; iComp < 3; iComp++)
    xMask->ReleaseVariationSurface(
      xVariations[iComp + xMask->GetNumberOfCoefficients()]);

    
  // Print the final report
  cout << "========================================" << endl;
  cout << "             FINAL REPORT               " << endl;
  cout << "----------------------------------------" << endl;

  TestGradientComputationPrintReport(
    xGlobalMaxPhiDiff, xGlobalMaxGRMS, xGlobalMaxBndDiff, 
    xGlobalMaxBndNrmDiff, xGlobalMaxBndAreaDiff, xGlobalMaxCellVolumeDiff,
    xGlobalTotalAreaDiff, xGlobalTotalVolumeDiff);

  cout << "----------------------------------------" << endl;
  cout << "Numeric Gradient CPU Secs. :" << tCentral.Read() << endl;
  cout << "Analytic Gradient CPU Secs.:" << tAnalytic.Read() << endl; 
  cout << "----------------------------------------" << endl;

  // The return value is based on any of the error terms exceeding epsilon
  int iPass = 
    (xGlobalMaxPhiDiff > eps || xGlobalMaxGRMS > eps ||
     xGlobalMaxBndDiff > eps || xGlobalMaxBndNrmDiff > eps ||
     xGlobalMaxBndAreaDiff > eps || xGlobalMaxCellVolumeDiff > eps ||
     xGlobalTotalAreaDiff > eps || xGlobalTotalVolumeDiff > eps) ? -1 : 0;

  // Describe what happened
  if(iPass == 0)
    cout << "TEST PASSED! (No error exceeds " << eps << ")" << endl;
  else
    cout << "TEST FAILED! (Errors exceed " << eps << ")" << endl;
  cout << "========================================" << endl;

  // Return value
  return iPass;
}

int TestOptimizerGradientComputation(
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

  // Keep track of the maximum error value
  double xMaxError = 0.0;
  
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
    double xDifference = fabs(xGradient[i] - dNumeric);
    cout << "Derivatives wrt. " << i << ". Numeric = " << dNumeric 
      << ";   Analytic = " << xGradient[i] 
      << ";   Difference = " << xDifference << endl;

    // Update the max error tracker
    UpdateMax(xDifference, xMaxError);
    }

  // Report difference in time
  cout << "Maximal Error : " << xMaxError << endl;
  cout << (xMaxError > eps ? "TEST FAILED" : "TEST PASSED") << endl;
  cout << "Central difference computed in " << tNumeric.Read() << " sec." << endl;
  cout << "Analytic gradient computed in  " << tAnalytic.Read() << " sec." << endl; 

  return xMaxError > eps ? -1 : 0;
}
