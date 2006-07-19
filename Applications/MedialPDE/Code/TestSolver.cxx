#include "ScriptInterface.h"
#include "CartesianMedialModel.h"
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
  double xMaxPhiDiff, double xMaxGradRDiff, double xMaxGRMS, double xMaxBndDiff, 
  double xMaxBndNrmDiff, double xMaxBndAreaDiff, double xMaxCellVolumeDiff,
  double xTotalAreaDiff, double xTotalVolumeDiff)
{
    // Print out the max difference
    cout << "Maximal differences between numeric and analytic derivatives:" << endl;
    cout << "Phi                      : " <<  xMaxPhiDiff << endl;
    cout << "Grad R                   : " <<  xMaxGradRDiff << endl;
    cout << "Sqr. Grad. Mag. of Phi   : " <<  xMaxGRMS << endl;
    cout << "Boundary Node Positions  : " <<  xMaxBndDiff << endl;
    cout << "Boundary Node Normals    : " <<  xMaxBndNrmDiff << endl;
    cout << "Boundary Area Elements   : " <<  xMaxBndAreaDiff << endl;
    cout << "Total Boundary Area      : " <<  xMaxBndAreaDiff << endl;
    cout << "Internal Volume Elements : " <<  xTotalAreaDiff << endl;
    cout << "Total Internal Volume    : " <<  xTotalVolumeDiff << endl;
}

int TestGradientComputation(
  GenericMedialModel *xSolver, CoefficientMapping *xMapping, 
  vnl_vector<double> P0, int nRandVar)
{
  size_t iVar;

  // Timers for different routines
  CodeTimer tCentral, tAnalytic, tInitialSolver;
  
  // The epsilon for central difference tests
  double eps = 0.0001;

  // The number of atoms in the solver
  size_t nAtoms = xSolver->GetNumberOfAtoms();

  // The number of parameters in the source coefficient space
  size_t nParams = xMapping->GetNumberOfParameters();
  size_t nCoeffs = xMapping->GetNumberOfCoefficients();
  
  // Max values for different classes of error
  double xGlobalMaxPhiDiff = 0.0, xGlobalMaxGradR = 0.0, xGlobalMaxGRMS = 0.0;
  double xGlobalMaxBndDiff = 0.0, xGlobalMaxBndNrmDiff = 0.0;
  double xGlobalMaxBndAreaDiff = 0.0, xGlobalMaxCellVolumeDiff = 0.0;
  double xGlobalTotalAreaDiff = 0.0, xGlobalTotalVolumeDiff = 0.0;

  cout << "*****************************************" << endl;
  cout << "Comparing Numeric Derivatives to Analytic" << endl;
  cout << "*****************************************" << endl;

  // Create a list of variational surfaces: some components, some random
  // variations
  vector<string> xVariationNames;
  vector< vnl_vector<double> > xVariations;

  // The initial coefficient vector is C0
  vnl_vector<double> C0 = xSolver->GetCoefficientArray();

  // Create the per-component variations
  size_t iComp;
  for( iComp = 0; iComp < nParams; iComp++)
    {
    ostringstream oss;
    oss << " Coefficient " << iComp;
    
    vnl_vector<double> xParamVariation(nParams, 0.0); xParamVariation[iComp] = 1.0;
    
    xVariationNames.push_back(oss.str());
    xVariations.push_back(xParamVariation);
    }

  // Create mixed variations
  for(size_t i = 0; i < nRandVar; i++)
    {    
    vnl_vector<double> dx(nParams);
    for(size_t j = 0; j < nParams; j++)
      dx[j] = rand() * 2.0 / RAND_MAX - 1.0;

    xVariationNames.push_back(" Random Variation ");
    xVariations.push_back(dx);
   }

  // For each variation, allocate a derivative array
  vector<MedialAtom *> dAtomArray;
  
  // Repeat for all variations
  for(iVar = 0; iVar < xVariations.size(); iVar++)
    {
    // Create an array of atoms to hold the derivative data
    MedialAtom *dAtoms = new MedialAtom[nAtoms];
    xSolver->PrepareAtomsForVariationalDerivative(xVariations[iVar], dAtoms);
    dAtomArray.push_back(dAtoms);
    }

  // Solve the medial PDE for all these variations
  tInitialSolver.Start();
  xSolver->ComputeAtoms();
  tInitialSolver.Stop();

  // Make three copies of the atom array
  MedialAtom *A0 = new MedialAtom[nAtoms];
  MedialAtom *A1 = new MedialAtom[nAtoms];
  MedialAtom *A2 = new MedialAtom[nAtoms];

  // Copy the current solution into the first array
  std::copy(xSolver->GetAtomArray(), xSolver->GetAtomArray() + nAtoms, A0);

  // Compute the analytical gradient of the atoms
  tAnalytic.Start();
  xSolver->ComputeAtomGradient(dAtomArray);
  tAnalytic.Stop();

  // Compute all the other derivatives
  for(iVar = 0; iVar < xVariations.size(); iVar++) 
    {
    cout << "----------------------------------------" << endl;
    cout << xVariationNames[iVar] << endl;
    cout << "----------------------------------------" << endl;
    
    // Solve for the forward difference
    xSolver->SetCoefficientArray(xMapping->Apply(C0, P0 + eps * xVariations[iVar]));
    tCentral.Start(); xSolver->ComputeAtoms(A0); tCentral.Stop();
    std::copy(xSolver->GetAtomArray(), xSolver->GetAtomArray() + nAtoms, A1);

    // Solve for the backward difference
    xSolver->SetCoefficientArray(xMapping->Apply(C0, P0 - eps * xVariations[iVar]));
    tCentral.Start(); xSolver->ComputeAtoms(A0); tCentral.Stop();
    std::copy(xSolver->GetAtomArray(), xSolver->GetAtomArray() + nAtoms, A2);

    // Compute the values of dPhi/dC between central difference and analytic
    // derivatives
    double xMaxPhiDiff = 0.0;
    double xMaxBndDiff = 0.0;
    double xMaxBndNrmDiff = 0.0;
    double xMaxGRMS = 0.0;
    double xMaxGradRDiff = 0.0;

    for(unsigned int i = 0; i < nAtoms; i++)
      {
      // Point to the atoms
      MedialAtom &a1 = A1[i];
      MedialAtom &a2 = A2[i];
      MedialAtom &a0 = dAtomArray[iVar][i];

      // Take the largest difference in Phi
      double dcPhi = 0.5 * (a1.F - a2.F) / eps;
      double daPhi = a0.F;
      UpdateMax(fabs(dcPhi - daPhi), xMaxPhiDiff);

      // Also look at the differences in grad R mag sqr.
      double dcGRMS = 0.5 * (a1.xGradRMagSqr - a2.xGradRMagSqr) / eps;
      double daGRMS = a0.xGradRMagSqr;
      UpdateMax(fabs(dcGRMS - daGRMS), xMaxGRMS);

      // Look at the gradR difference
      SMLVec3d dcGradR = 0.5 * (a1.xGradR - a2.xGradR) / eps;
      SMLVec3d daGradR = a0.xGradR;
      UpdateMax((dcGradR - daGradR).two_norm(), xMaxGradRDiff);

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

    cout << endl;

    // Now, test the effect on SolutionData objects
    MedialIterationContext *xGrid = xSolver->GetIterationContext();
    SolutionData sd0(xGrid, A0);
    SolutionData sd1(xGrid, A1);
    SolutionData sd2(xGrid, A2);
    PartialDerivativeSolutionData pdsd(&sd0, dAtomArray[iVar]);

    // Compute the boundary weights derivatives
    sd0.UpdateBoundaryWeights();
    sd1.UpdateBoundaryWeights();
    sd2.UpdateBoundaryWeights();
    pdsd.UpdateBoundaryWeights();

    // Compute the area weight difference
    double xMaxBndAreaDiff = 0.0;
    for(MedialBoundaryPointIterator itBnd(xGrid); !itBnd.IsAtEnd(); ++itBnd)
      {
      double aNumeric = 0.5 * ( 
        sd1.xBoundaryWeights[itBnd.GetIndex()] - 
        sd2.xBoundaryWeights[itBnd.GetIndex()]) / eps;
      double aAnalytic = pdsd.xBoundaryWeights[itBnd.GetIndex()];
      UpdateMax(fabs(aAnalytic - aNumeric), xMaxBndAreaDiff);
      }

    // Compute the difference in internal point volume weights
    size_t nInternal = 4;
    sd0.UpdateInternalWeights(nInternal);
    sd1.UpdateInternalWeights(nInternal);
    sd2.UpdateInternalWeights(nInternal);
    pdsd.UpdateInternalWeights(nInternal);

    // Compute the area weight difference
    double xMaxCellVolumeDiff = 0.0;
    for(MedialInternalPointIterator itPoint(xGrid,nInternal); !itPoint.IsAtEnd(); ++itPoint)
      {
      double aNumeric = 0.5 * ( 
        sd1.xInternalWeights[itPoint.GetIndex()] - 
        sd2.xInternalWeights[itPoint.GetIndex()]) / eps;
      double aAnalytic = pdsd.xInternalWeights[itPoint.GetIndex()];
      UpdateMax(fabs(aAnalytic - aNumeric), xMaxCellVolumeDiff);
      }

    // Compare total area and volume derivatives
    double xTotalVolumeDiff = 
      fabs(0.5 * (sd1.xInternalVolume - sd2.xInternalVolume) / eps 
        - pdsd.xInternalVolume);

    double xTotalAreaDiff = 
      fabs(0.5 * (sd1.xBoundaryArea - sd2.xBoundaryArea) / eps 
        - pdsd.xBoundaryArea);

    // Print the report
    /*
    TestGradientComputationPrintReport(
      xMaxPhiDiff, xMaxGradRDiff, xMaxGRMS, xMaxBndDiff, 
      xMaxBndNrmDiff, xMaxBndAreaDiff, xMaxCellVolumeDiff,
      xTotalAreaDiff, xTotalVolumeDiff); */

    // Update the global maximum counters
    UpdateMax(xMaxPhiDiff, xGlobalMaxPhiDiff);
    UpdateMax(xMaxGRMS, xGlobalMaxGRMS);
    UpdateMax(xMaxGradRDiff, xGlobalMaxGradR);
    UpdateMax(xMaxBndDiff, xGlobalMaxBndDiff);
    UpdateMax(xMaxBndNrmDiff, xGlobalMaxBndNrmDiff);
    UpdateMax(xMaxBndAreaDiff, xGlobalMaxBndAreaDiff);
    UpdateMax(xMaxCellVolumeDiff, xGlobalMaxCellVolumeDiff);
    UpdateMax(xTotalAreaDiff, xGlobalTotalAreaDiff);
    UpdateMax(xTotalVolumeDiff, xGlobalTotalVolumeDiff);

    // Delete the atoms
    delete dAtomArray[iVar];
    }

  // Delete the atom arrays
  delete A0; delete A1; delete A2;

    
  // Print the final report
  cout << "========================================" << endl;
  cout << "             FINAL REPORT               " << endl;
  cout << "----------------------------------------" << endl;

  TestGradientComputationPrintReport(
    xGlobalMaxPhiDiff, xGlobalMaxGradR, xGlobalMaxGRMS, xGlobalMaxBndDiff, 
    xGlobalMaxBndNrmDiff, xGlobalMaxBndAreaDiff, xGlobalMaxCellVolumeDiff,
    xGlobalTotalAreaDiff, xGlobalTotalVolumeDiff);

  cout << "----------------------------------------" << endl;
  cout << "Numeric Gradient CPU Secs. :" << tCentral.Read() << endl;
  cout << "Analytic Gradient CPU Secs.:" << tAnalytic.Read() << endl; 
  cout << "Initial Solution CPU Secs.:" << tInitialSolver.Read() << endl; 
  cout << "----------------------------------------" << endl;

  // The return value is based on any of the error terms exceeding epsilon
  int iPass = 
    (xGlobalMaxPhiDiff > eps || xGlobalMaxGRMS > eps ||
     xGlobalMaxGradR > eps ||
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
  CoefficientMapping &xMapping,
  GenericMedialModel *xSolver)
{
  // Set up timers
  CodeTimer tAnalytic, tNumeric;

  // Get the number of parameters that are being optimized over
  size_t nParams = xMapping.GetNumberOfParameters();

  // Create the initial parameter vector and the vector to hold the gradient
  vnl_vector<double> P(nParams, 0.0), xGradient(nParams, 0.0);

  // We do not want to only perform the test at P = 0 because that may not
  // detect problems for other P values. So instead, we can move along the
  // gradient a certain direction.
  mop.ComputeGradient(P.data_block(), xGradient.data_block());
  P += 0.1 * xGradient;

  // Now we are ready to perform the test.

  // Compute the analytic gradient
  tAnalytic.Start();
  mop.ComputeGradient(P.data_block(), xGradient.data_block());
  tAnalytic.Stop();
  
  // Keep track of the maximum error value
  double xMaxError = 0.0;
  
  // Compute the numeric gradient of the objective function
  double eps = 0.01;
  for(size_t i = 0; i < nParams; i++)
    {
    // Get the current value of the coefficient
    double ci = P[i];

    // Compute the objective function at point 1
    P[i] = ci + eps;
    tNumeric.Start();
    double f1 = mop.Evaluate(P.data_block());
    tNumeric.Stop();

    // Compute the objective at point 2
    P[i] = ci - eps;
    tNumeric.Start();
    double f2 = mop.Evaluate(P.data_block());
    tNumeric.Stop();

    // Restore the coefficient
    P[i] = ci;
    
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
