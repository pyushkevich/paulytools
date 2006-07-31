#include "ScriptInterface.h"
#include "CartesianMedialModel.h"
#include "BasisFunctions2D.h"
#include "ITKImageWrapper.h"
#include "MedialPDERenderer.h"
#include "MedialAtomGrid.h"
#include "OptimizationTerms.h"
#include "Procrustes.h"
#include "Registry.h"
#include <iostream>
#include <iomanip>
#include <ctime>

using namespace std;

// An inline function to update max counters
inline void UpdateMax(double xValue, double &xMaxValue)
{
  if(xValue > xMaxValue)
    xMaxValue = xValue;
}

class DerivativeTestQuantities {
public:
  DerivativeTestQuantities()
    { slenmax = 0; }

  void Update(const std::string &s, double v)
    {
    MapIter it = xMap.find(s);
    if(it == xMap.end())
      {
      xMap.insert(std::make_pair(s, v));
      slenmax = std::max(slenmax, s.size());
      }
    else
      it->second = std::max(it->second, v);
    }

  void Update(const std::string &s, double dv, double v1, double v2, double eps)
    {
    this->Update(s, fabs(dv - (v1 - v2) * 0.5 / eps)); 
    }

  void Update(const std::string &s, SMLVec3d dv, SMLVec3d v1, SMLVec3d v2, double eps)
    {
    this->Update(s, (dv - (v1 - v2) * 0.5 / eps).two_norm()); 
    }

  void Update(const DerivativeTestQuantities &src)
    {
    for(ConstMapIter it = src.xMap.begin(); it != src.xMap.end(); ++it)
      this->Update(it->first, it->second);
    }

  bool Test(double xMaxError)
    {
    for(ConstMapIter it = xMap.begin(); it != xMap.end(); ++it)
      if(it->second > xMaxError) return false;
    return true;
    }

  void PrintReport() 
    {
    for(ConstMapIter it = xMap.begin(); it != xMap.end(); ++it)
      cout << std::setw(slenmax) << it->first << " : " << it->second << endl;
    }

private:
  typedef std::map<std::string, double> MapType;
  typedef MapType::iterator MapIter;
  typedef MapType::const_iterator ConstMapIter;
  MapType xMap;

  size_t slenmax;
};

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
  
  // Holder for all the test quantities. This object will report the maximum
  // value for each quantity over all the tests
  DerivativeTestQuantities dtqGlobal;

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
    // Transform the variation in parameters to a variation in coefficients
    vnl_vector<double> dC = 
      xMapping->ApplyJacobianInParameters(C0, P0, xVariations[iVar]);

    // Create an array of atoms to hold the derivative data
    MedialAtom *dAtoms = new MedialAtom[nAtoms];
    xSolver->PrepareAtomsForVariationalDerivative(dC, dAtoms);
    dAtomArray.push_back(dAtoms);
    }

  // Set the coefficients for the initial solution
  xSolver->SetCoefficientArray(xMapping->Apply(C0, P0));

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
  cout << endl;

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

    // Define an accumulator for all atom-wise errors
    DerivativeTestQuantities dtq;

    // Compute atom-wise errors
    for(unsigned int i = 0; i < nAtoms; i++)
      {
      // Point to the atoms
      MedialAtom &a1 = A1[i];
      MedialAtom &a2 = A2[i];
      MedialAtom &a0 = dAtomArray[iVar][i];

      // Look at such a simple thing as a.X
      dtq.Update("Atom's X", a0.X, a1.X, a2.X, eps);
      dtq.Update("Atom's Xu", a0.Xu, a1.Xu, a2.Xu, eps);
      dtq.Update("Atom's Xv", a0.Xv, a1.Xv, a2.Xv, eps);
      dtq.Update("Atom's Xuu", a0.Xuu, a1.Xuu, a2.Xuu, eps);
      dtq.Update("Atom's Xuv", a0.Xuv, a1.Xuv, a2.Xuv, eps);
      dtq.Update("Atom's Xvv", a0.Xvv, a1.Xvv, a2.Xvv, eps);

      // Look at the difference in the specific terms of metric tensor
      for(size_t u = 0; u < 2; u++) for(size_t v = 0; v < 2; v++)
        {
        dtq.Update("Metric Tensor (Covariant)", 
          a0.G.xCovariantTensor[u][v],
          a1.G.xCovariantTensor[u][v],
          a2.G.xCovariantTensor[u][v], eps);

        dtq.Update("Metric Tensor (Contra)", 
          a0.G.xContravariantTensor[u][v],
          a1.G.xContravariantTensor[u][v],
          a2.G.xContravariantTensor[u][v], eps);
        }

      // Take the largest difference in Phi
      dtq.Update("Phi", a0.F, a1.F, a2.F, eps);

      // Also look at the differences in grad R mag sqr.
      dtq.Update("Sqr. Grad. Mag. of R", 
        a0.xGradRMagSqr, a1.xGradRMagSqr, a2.xGradRMagSqr, eps);

      // Look at the gradR difference
      dtq.Update("Grad R", a0.xGradR, a1.xGradR, a2.xGradR, eps);

      // Take the largest difference in boundary derivs
      // for(unsigned int k = 0; k < 2; k++)
      //   {
      //   dtq.Update("Boundary Node Positions", a0.xBnd[k].X, a1.xBnd[k].X, a2.xBnd[k].X, eps);
      //   dtq.Update("Boundary Node Normals", a0.xBnd[k].N, a1.xBnd[k].N, a2.xBnd[k].N, eps);
      //   }
      }

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
    for(MedialBoundaryPointIterator it(xGrid); !it.IsAtEnd(); ++it)
      {
      dtq.Update("Boundary Area Element",
        pdsd.xBoundaryWeights[it.GetIndex()],
        sd1.xBoundaryWeights[it.GetIndex()],
        sd2.xBoundaryWeights[it.GetIndex()], eps);

      dtq.Update("Boundary Node Positions",
        GetBoundaryPoint(it, dAtomArray[iVar]).X,
        GetBoundaryPoint(it, A1).X,
        GetBoundaryPoint(it, A2).X, eps);

      dtq.Update("Boundary Node Normals",
        GetBoundaryPoint(it, dAtomArray[iVar]).N,
        GetBoundaryPoint(it, A1).N,
        GetBoundaryPoint(it, A2).N, eps);
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
      dtq.Update("Internal Volume Element",
        pdsd.xInternalWeights[itPoint.GetIndex()],
        sd1.xInternalWeights[itPoint.GetIndex()],
        sd2.xInternalWeights[itPoint.GetIndex()], eps);

    // Print the report
    cout << "Maximal differences between numeric and analytic derivatives:" << endl;
    dtq.PrintReport();

    // Update the global maximum counters
    dtqGlobal.Update(dtq);

    // Delete the atoms
    delete dAtomArray[iVar];
    }

  // Delete the atom arrays
  delete A0; delete A1; delete A2;

    
  // Print the final report
  cout << "========================================" << endl;
  cout << "             FINAL REPORT               " << endl;
  cout << "----------------------------------------" << endl;

  cout << "Maximal differences between numeric and analytic derivatives:" << endl;
  dtqGlobal.PrintReport();

  cout << "----------------------------------------" << endl;
  cout << "Numeric Gradient CPU Secs. :" << tCentral.Read() << endl;
  cout << "Analytic Gradient CPU Secs.:" << tAnalytic.Read() << endl; 
  cout << "Initial Solution CPU Secs.:" << tInitialSolver.Read() << endl; 
  cout << "----------------------------------------" << endl;

  // Check if all the errors are reasonable
  bool flagPass = dtqGlobal.Test(eps);

  // Describe what happened
  if(flagPass)
    cout << "TEST PASSED! (No error exceeds " << eps << ")" << endl;
  else
    cout << "TEST FAILED! (Errors exceed " << eps << ")" << endl;
  cout << "========================================" << endl;

  // Return value
  return flagPass ? 0 : -1;
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
  mop.PrintReport(cout);

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
