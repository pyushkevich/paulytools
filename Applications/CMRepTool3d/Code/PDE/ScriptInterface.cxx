#include "ScriptInterface.h"
#include "MedialPDESolver.h"
#include "BasisFunctions2D.h"
#include "ITKImageWrapper.h"
#include "MedialPDERenderer.h"
#include "OptimizationTerms.h"
#include "Procrustes.h"
#include "Registry.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "vnl/vnl_matlab_filewrite.h"
#include "ConjugateGradientMethod.h"
#include "EvolutionaryStrategy.h"

#include <iostream>
#include <fstream>

void ExportMedialMeshToVTK(MedialAtomGrid *xGrid, MedialAtom *xAtoms, const char *file);
void ExportBoundaryMeshToVTK(MedialAtomGrid *xGrid, MedialAtom *xAtoms, const char *file);
//void ExportIntensityFieldToVTK(MedialAtomGrid *xGrid, MedialAtom *xAtoms, 
//  ITKImageWrapper<float> *imgField, const char *file);
using namespace std;


namespace medialpde {

struct DiscreteAtom 
{ 
  double u, v, x, y, z; 
  unsigned int iu, iv; 
};

/***************************************************************************
 * MedialPDE Code
 * ----------------
 *  This is the main application code
 **************************************************************************/
MedialPDE::MedialPDE(unsigned int nBasesU, unsigned int nBasesV, 
  unsigned int xResU, unsigned int xResV)
{
  xSurface = new FourierSurface(nBasesU, nBasesV);
  xSolver = new MedialPDESolver(xResU, xResV);
  xSolver->SetMedialSurface(xSurface);

  eMask = FULL; 
  eOptimizer = GRADIENT;
  eMatch = VOLUME;
  
  xStepSize = 0.1;
  xCoarsenessRho = 1.0;
  xCoarsenessX = 1.0;
  xMeshDumpImprovementPercentage = 0.0;

  flagIntensityPresent = false;
}

MedialPDE::~MedialPDE()
{
  delete xSolver;
  delete xSurface;
}

void MedialPDE::LoadFromDiscreteMRep(const char *file, double xRhoInit)
{
  ifstream fin(file, ios_base::in);
  int iu, iv, iend;
  double x, y, z, d;

  // Vector to store atoms
  typedef vector<DiscreteAtom> AList;
  AList xAtoms;

  // Get the max u and v values
  unsigned int uMax = 0, vMax = 0;

  // Read atoms
  bool done = false;
  while(!done)
    {
    DiscreteAtom atom;

    // Read the atom
    fin >> atom.iu; fin >> atom.iv; fin >> iend; 
    fin >> atom.x; fin >> atom.y; fin >> atom.z; 
    fin >> d; fin >> d; fin >> d; fin >> d; 
    fin >> d; fin >> d; fin >> d;

    if(uMax < atom.iu) uMax = atom.iu;
    if(vMax < atom.iv) vMax = atom.iv;

    if(!fin.good())
      { done = true; }
    else
      {
      xAtoms.push_back(atom);
      }
    }

  // Scale by u and v to unit square
  for(AList::iterator it = xAtoms.begin(); it!=xAtoms.end(); ++it)
    { it->u = it->iu * 1.0 / uMax; it->v = it->iv * 1.0 / vMax; }

  // Create double arrays
  double *xx = new double[xAtoms.size()];
  double *yy = new double[xAtoms.size()];
  double *zz = new double[xAtoms.size()];
  double *uu = new double[xAtoms.size()];
  double *vv = new double[xAtoms.size()];
  double *rr = new double[xAtoms.size()];

  for(unsigned int i = 0; i < xAtoms.size(); i++)
    {
    xx[i] = xAtoms[i].x;
    yy[i] = xAtoms[i].y;
    zz[i] = xAtoms[i].z;
    uu[i] = xAtoms[i].u;
    vv[i] = xAtoms[i].v;
    rr[i] = xRhoInit;
    }

  // Perform the fitting on x, y and z
  xSurface->FitToData(xAtoms.size(), 0, uu, vv, xx);
  xSurface->FitToData(xAtoms.size(), 1, uu, vv, yy);
  xSurface->FitToData(xAtoms.size(), 2, uu, vv, zz);

  // Fit the rho function to a constant
  xSurface->FitToData(xAtoms.size(), 3, uu, vv, rr);

  // Clean up
  delete xx; delete yy; delete zz; delete uu; delete vv; delete rr;
}

void MedialPDE::SaveToParameterFile(const char *file)
{
  // Use the registry to save the data
  Registry R;

  // Store the type of the m-rep specification
  R["Grid.Type"] << "Cartesian";
  R["Grid.Size.U"] << xSolver->GetNumberOfUPoints();
  R["Grid.Size.V"] << xSolver->GetNumberOfVPoints();

  // Store the fourier coefficient information
  R["SurfaceModel"] << "Fourier";
  xSurface->SaveToRegistry(R.Folder("Fourier"));

  // Save the registry
  R.WriteToFile(file);
}

void MedialPDE::LoadFromParameterFile(const char *file)
{
  // Use the registry to save the data
  Registry R(file);

  // Store the type of the m-rep specification
  if(R["Grid.Type"][""] == Registry::StringType("Cartesian"))
    {
    unsigned int m = R["Grid.Size.U"][64];
    unsigned int n = R["Grid.Size.V"][64];

    // Store the fourier coefficient information
    if(R["SurfaceModel"][""] == Registry::StringType("Fourier"))
      {
      // Read the surface from the parameters
      xSurface->ReadFromRegistry(R.Folder("Fourier"));

      // Solve the PDE
      xSolver->Solve();
      }
    else
      { cerr << "Invalid surface model in file" << endl; }
    } 
  else
    { cerr << "Invalid grid type in file" << endl; }
}

void MedialPDE::Solve()
{
  xSolver->Solve();
}

void MedialPDE::SetNumberOfCoefficients(unsigned int m, unsigned int n)
{
  xSurface->SetNumberOfCoefficients(m, n);
  xSolver->SetMedialSurface(xSurface);
}

void MedialPDE::GenerateSampleModel()
{
  // Decide how many points to interpolate
  unsigned int nSide = 21, nPoints = nSide * nSide;
  double uStep = 1.0 / (nSide - 1);

  // Allocate arrays of points and coordinates
  double xPoints[nPoints], yPoints[nPoints], zPoints[nPoints], rhoPoints[nPoints];
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
      rhoPoints[i] = -0.35;
      ++i;
      }

  // Peform the fit
  xSurface->FitToData(nPoints, 0, uPoints, vPoints, xPoints);
  xSurface->FitToData(nPoints, 1, uPoints, vPoints, yPoints);
  xSurface->FitToData(nPoints, 2, uPoints, vPoints, zPoints);
  xSurface->FitToData(nPoints, 3, uPoints, vPoints, rhoPoints);
}

/*
double ComputeImageMatchGradient(
  MedialPDESolver *xSolver, PDEFourierWrapper *xModel, FourierSurface *xSurface,
  FloatImage *image, vnl_vector<double> &xGrad)
{
  unsigned int i, j, k;

  // Get the atom grid
  MedialAtomGrid *xAtomGrid = xSolver->GetAtomGrid();
  MedialAtom *xAtoms = xSolver->GetAtomArray();

  // Create a solution data object representing current solution
  SolutionData *S0 = new SolutionData(xSolver, true);

  // Solve the medial PDE for a set of finite difference offsets
  size_t nCoeff = xSurface->GetNumberOfRawCoefficients();
  double *xCoeff = xSurface->GetRawCoefficientArray();

  // Create an array of solution data objects
  SolutionData **SGrad = new SolutionData *[nCoeff];

  // Use current solution as the guess
  xSolver->SetSolutionAsInitialGuess();

  // Compute the time it takes to compute the derivative m-reps
  double tStart = clock();

  // Repeat for each coefficient
  MedialBoundaryPointIterator *itBnd = xAtomGrid->NewBoundaryPointIterator();
  for(size_t iCoeff = 0; iCoeff < nCoeff; iCoeff++)
    {
    // Increment the coefficient by epsilon
    double xOldValue = xCoeff[iCoeff];
    xCoeff[iCoeff] = xOldValue + xEpsilon;

    // Solve the PDE for the new value
    xSolver->Solve(xModel, 1.0e-13); cout << "." << flush;

    // Compute the solution data
    SGrad[iCoeff] = new SolutionData(xSolver, true);
    }

  // Restore the solution (should be immediate)
  xSolver->Solve(xModel); cout << endl;

  // See how much time elapsed for the medial computation
  double tMedial = clock() - tStart;

  // At this point, we have computed the solution at current X and at a set of 
  // X + eps positions. We can use this to compute the gradient of each of the
  // terms involved in the optimization
  vnl_vector<double> xGradImage(nCoeff), xGradJacobian(nCoeff);
  double xSolution, xImageTerm, xJacobianTerm;
  
  // Compute the gradient of the image match term
  BoundaryImageMatchTerm termImage(image);
  xImageTerm = termImage.ComputeGradient(S0, SGrad, xEpsilon, xGradJacobian);

  // Compute the gradient of the jacobian penalty term
  BoundaryJacobianEnergyTerm termJacobian();
  xJacobianTerm = termJacobian.ComputeGradient(S0, SGrad, xEpsilon, xGradJacobian); 

  // Finish timing
  double tGradient = clock() - tMedial;
  cout << "PDE Time: " << (tMedial / CLOCKS_PER_SEC) <<
   "; Gradient Time: " << (tGradient / CLOCKS_PER_SEC) << endl;

  // Clean up everything
  for(size_t iCoeff = 0; iCoeff < nCoeff; iCoeff++)
    delete SGrad[iCoeff];
  delete S0; delete SGrad;

  // Compute the gradient and the solution
  xGrad = xGradImage + 0.01 * xGradJacobian;
  return xImageTerm + 0.01 * xJacobianTerm;
}
*/

double MedialPDE::ComputeImageMatch(FloatImage *image)
{
  // Create a solution object (reference atoms in the solver)
  SolutionData S(xSolver, false);

  // Create an image match term
  BoundaryImageMatchTerm termImage(image);

  // Compute the energy
  double xMatch = termImage.ComputeEnergy(&S);
  
  // Print a report
  cout << "REPORT: " << endl;
  termImage.PrintReport(cout);

  return xMatch;
}

void GradientDescentOptimization(MedialOptimizationProblem *xProblem, 
  vnl_vector<double> &xSolution, unsigned int nSteps, double xStep)
{
  // Create the initial solution
  size_t nCoeff = xSolution.size();
  vnl_vector<double> xGradient(nCoeff, 0.0);

  // Print report information  
  // ofstream fdump("conjgrad.txt",ios_base::out);
  // fdump << "CONJUGATE GRADIENT OPTIMIZER DUMP" << endl;

  // Iterate, moving in the gradient direction
  for(unsigned int p = 0; p < nSteps; p++)
    {
    // Compute the gradient and the image match
    double xMatch = 
      xProblem->ComputeGradient( xSolution.data_block(), xGradient.data_block());

    // Move in the gradient direction
    xSolution -= xStep * xGradient;

    // Report the current state
    // fdump << "STEP " << p << endl;
    // xProblem->PrintReport(fdump);

    // cout << "Step " << p << ", best value: " << xMatch << endl;
    }
}


void MedialPDE::ExportIterationToVTK(unsigned int iter)
{
  ostringstream oss;
  oss << strDumpPath << ".";
  oss.width(5);
  oss.fill('0');
  oss << iter;
  string fMedial = oss.str() + ".med.vtk";
  string fBoundary = oss.str() + ".bnd.vtk";

  ExportMedialMeshToVTK(
    xSolver->GetAtomGrid(), xSolver->GetAtomArray(), fMedial.c_str());
  ExportBoundaryMeshToVTK(
    xSolver->GetAtomGrid(), xSolver->GetAtomArray(), fBoundary.c_str());
}

void MedialPDE::SaveVTKMesh(const char *fMedial, const char *fBnd)
{
  ExportMedialMeshToVTK(
    xSolver->GetAtomGrid(), xSolver->GetAtomArray(), fMedial);
  ExportBoundaryMeshToVTK(
    xSolver->GetAtomGrid(), xSolver->GetAtomArray(), fBnd);
}

void MedialPDE::ConjugateGradientOptimization(
  MedialOptimizationProblem *xProblem, 
  vnl_vector<double> &xSolution, unsigned int nSteps, double xStep)
{
  size_t nCoeff = xSolution.size();

  // Construct the conjugate gradient optimizer
  ConjugateGradientMethod xMethod(*xProblem, Vector(nCoeff, xSolution.data_block()));
  xMethod.setStepSize(xStep);
  // xMethod.setAdaptableStepSize(true);
  xMethod.setBrentStepTolerance(2.0e-3);

  // Debugging info
  // ofstream fdump("conjgrad.txt",ios_base::out);
  // fdump << "CONJUGATE GRADIENT OPTIMIZER DUMP" << endl;

  // When we last dumped a mesh
  double xLastMeshDumpValue = xMethod.getBestEverValue();
  if(xMeshDumpImprovementPercentage > 0.0)
    ExportIterationToVTK(0);
  
  for(size_t p = 0; p < nSteps; p++)
    {
    if(xMethod.isFinished())
      {
      cout << "CONJGRAD finished on its own!" << endl;
      break;
      }
    xMethod.performIteration();

    // fdump << "STEP " << p << endl;
    // xProblem->PrintReport(fdump);

    cout << "Step " << p << ", best value: " << xMethod.getBestEverValue() << endl;

    // Dump the mesh
    if(xMeshDumpImprovementPercentage > 0.0)
      if(xMethod.getBestEverValue() < 
        (1.0 - xMeshDumpImprovementPercentage) * xLastMeshDumpValue)
        {
        xProblem->evaluate(xMethod.getBestEverX());
        ExportIterationToVTK(p);
        xLastMeshDumpValue = xMethod.getBestEverValue(); 

        xProblem->PrintReport(cout);
        }

    // Save the file 
    //ostringstream oss;
    //oss << "iter_" << p;
    //SaveToParameterFile(oss.str().c_str());
    }

  // fdump.close();

  // Store the best result
  xSolution.copy_in(xMethod.getBestEverX().getDataArray());
}

inline vnl_vector<double> VectorCast(const pauly::Vector &Y)
{
  return vnl_vector<double>(Y.getDataArray(), Y.size());
}

inline pauly::Vector VectorCast(vnl_vector<double> &Y)
{
  return pauly::Vector(Y.size(), Y.data_block());
}

void EvolutionaryOptimization(
  MedialOptimizationProblem *xProblem, 
  vnl_vector<double> &xSolution, 
  unsigned int nSteps)
{
  // Create the initial solution
  size_t nCoeff = xSolution.size();
  vnl_vector<double> xSigma(nCoeff, 0.01);
  
  // Create the initial search space
  GaussianSS xSearch(VectorCast(xSolution), VectorCast(xSigma));
  
  // Construct the evolutionary optimizer
  EvolutionaryStrategy xMethod(*xProblem, xSearch, 2, 4, SELECTION_MuPlusLambda);

  // Set the factor by which the sigmas increase
  xMethod.setSigmaFactor(1.0);

  // Add the inital solution to the mix
  xMethod.setNth(0, xSearch.getMean());

  // Set the delta-sigma vector
  vnl_vector<double> xDeltaSigma(nCoeff, 0.1);
  xMethod.setDeltaSigma(VectorCast(xDeltaSigma));
  
  // Run the optimization
  for(size_t p = 0; p < nSteps; p++)
    {
    xMethod.performIteration();
    cout << "Best Value: " << xMethod.getBestEverValue() << endl;
    }

  // Store the best result
  xSolution.copy_in(xMethod.getBestEverX().getDataArray());
}

void MedialPDE
::RunOptimization(FloatImage *image, unsigned int nSteps)
{
  // Create a coefficient mask
  IMedialCoefficientMask *xMask;
  if(eMask == AFFINE) 
    {
    xMask = new AffineTransformCoefficientMask(xSurface);
    }
  else if(eMask == FULL) 
    {
    xMask = new PassThroughCoefficientMask(xSurface);
    }
  else
    {
    double crs[4] = { 
      xCoarsenessX, xCoarsenessX, xCoarsenessX, xCoarsenessRho };
    xMask = xSurface->NewCoarseToFineCoefficientMask(crs);
    }

  // Create the optimization problem
  MedialOptimizationProblem xProblem(xSolver, xMask);

  // Create an image match term and a jacobian term
  EnergyTerm *xTermImage;
  if(eMatch == VOLUME)
    xTermImage = new VolumeOverlapEnergyTerm(image, 8);
  else
    xTermImage = new BoundaryImageMatchTerm(image);
  
  // Create the other terms
  BoundaryJacobianEnergyTerm xTermJacobian;
  CrestLaplacianEnergyTerm xTermCrest;
  AtomBadnessTerm xTermBadness;
  MedialRegularityTerm xTermRegularize(xSolver->GetAtomArray(), xSolver->GetAtomGrid());

  // Add the terms to the problem
  xProblem.AddEnergyTerm(xTermImage, 1.0);
  xProblem.AddEnergyTerm(&xTermJacobian, 0.005);
  xProblem.AddEnergyTerm(&xTermCrest, 0.005);
  xProblem.AddEnergyTerm(&xTermBadness, 1.0);
  xProblem.AddEnergyTerm(&xTermRegularize, 1.0);

  // Create the initial solution
  size_t nCoeff = xMask->GetNumberOfCoefficients();
  vnl_vector<double> xSolution(xMask->GetCoefficientArray(), nCoeff);
  
  // Initial solution report
  cout << "INITIAL SOLUTION REPORT: " << endl;
  xProblem.Evaluate(xSolution.data_block());
  xProblem.PrintReport(cout);

  // At this point, split depending on the method
  if(eOptimizer == CONJGRAD)
    ConjugateGradientOptimization(&xProblem, xSolution, nSteps, xStepSize);
  else if(eOptimizer == GRADIENT)
    GradientDescentOptimization(&xProblem, xSolution, nSteps, xStepSize);
  else
    EvolutionaryOptimization(&xProblem, xSolution, nSteps);

  // Store the best result
  xMask->SetCoefficientArray(xSolution.data_block());
  xSolver->Solve();

  // Delete the mask
  delete xMask;
  delete xTermImage;
}
class FirstMomentComputer : public EuclideanFunction
{
public:
  FirstMomentComputer(unsigned int i)
    { this->i = i; }

  double Evaluate(const SMLVec3d &X)
    { return X[i]; }

private:
  unsigned int i;
};

class SecondMomentComputer : public EuclideanFunction
{
public:
  SecondMomentComputer(unsigned int i, unsigned int j)
    { this->i = i; this->j = j;}

  double Evaluate(const SMLVec3d &X)
    { return X[i] * X[j]; }

private:
  unsigned int i, j;
};

template<typename TPixel>
class ImageVolumeMatch : public EuclideanFunction
{
public:
  ImageVolumeMatch(ITKImageWrapper<TPixel> *xWrapper)
    { this->xWrapper = xWrapper; }

  double Evaluate(const SMLVec3d &X)
    { return xWrapper->Interpolate(X[0], X[1], X[2], 0.0f); }

private:
  ITKImageWrapper<TPixel> *xWrapper;
};

/** Match an mrep to a floating point (level set) image using moments. */
void MedialPDE::MatchImageByMoments(FloatImage *image, unsigned int nCuts)
{
  // Typedefs
  typedef vnl_matrix_fixed<double,3,3> Mat;
  typedef SMLVec3d Vec;

  // Compute the mean and covariance matrix of the non-zero voxels in the
  // image.
  Mat xCov(0.0f); Vec xMean(0.0f); double xVolume = 0.0;
  size_t i, j, k, f;

  // Iterate over all voxels in the image
  typedef FloatImage::WrapperType::ImageType ImageType;
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;
  ImageType::Pointer xImage = image->xImage->GetInternalImage();
  IteratorType it(xImage, xImage->GetBufferedRegion());
  for( ; !it.IsAtEnd(); ++it)
    {
    // Get the volume of this voxel
    if(it.Get() > 0.0)
      {      
      // Get the spatial position of the point
      itk::ContinuousIndex<double, 3> iRaw(it.GetIndex());
      itk::Point<double, 3> ptPosition;
      xImage->TransformContinuousIndexToPhysicalPoint(iRaw, ptPosition);
      SMLVec3d xPosition(ptPosition.GetDataPointer());

      // Add to the mean and 'covariance'
      xVolume += 1.0;
      xMean += xPosition;
      xCov += outer_product(xPosition, xPosition);
      }
    }

  // Compute the actual volume of the image 
  double xPhysicalVolume = xVolume * 
    xImage->GetSpacing()[0] * xImage->GetSpacing()[1] * xImage->GetSpacing()[2];

  // Scale the mean and covariance by accumulated weights
  xMean = xMean / xVolume;
  xCov = (xCov - xVolume *  outer_product(xMean, xMean)) / xVolume;
 
  // Compute the mean and covariance of the image
  cout << "--- MATCHING BY MOMENTS ---" << endl;
  cout << "Image Volume: " << xPhysicalVolume << endl;
  cout << "Image Mean: " << xMean << endl;
  cout << "Image Covariance: " << endl << xCov << endl;

  // Now compute the same statistics for the m-rep (compute coordinates)
  Mat yCov(0.0); Vec yMean(0.0); double yVolume = 0.0;

  // Create a solution data object to speed up computations
  SolutionData S0(xSolver, false);
  S0.UpdateInternalWeights(nCuts);
  
  for(i = 0; i < 3; i++)
    {
    // Compute the first moment in i-th direction
    FirstMomentComputer fmc(i);
    yMean[i] = IntegrateFunctionOverInterior(
      S0.xAtomGrid, S0.xInternalPoints, S0.xInternalWeights, nCuts, &fmc);

    for(j = 0; j <= i; j++)
      {
      // Compute the second moment in i-th and j-th directions
      SecondMomentComputer smc(i,j);
      yCov[i][j] = yCov[j][i] = IntegrateFunctionOverInterior(
        S0.xAtomGrid, S0.xInternalPoints, S0.xInternalWeights, nCuts, &smc);
      }
    }

  // Compute the actual mean and covariance
  yVolume = S0.xInternalVolume;
  yMean /= yVolume;
  yCov = (yCov - yVolume * outer_product(yMean, yMean)) / yVolume;
  
  cout << "Model Volume: " << yVolume << endl;
  cout << "Model Mean: " << yMean << endl;
  cout << "Model Covariance: " << endl << yCov << endl;

  // Now, compute the alignment that will line up the moments
  // Decompose xCov and yCov into eigenvalues
  vnl_vector<double> Dx(3), Dy(3);
  vnl_matrix<double> Vx(3,3), Vy(3,3);
  vnl_symmetric_eigensystem_compute(xCov, Vx, Dx);
  vnl_symmetric_eigensystem_compute(yCov, Vy, Dy);

  // Compute the scale factor
  double s = sqrt(dot_product(Dx,Dy) / dot_product(Dy,Dy));

  // The best set of coefficients and the associated match value
  size_t iBestSurface; double xBestMatch;

  // Get the numbers of coefficients
  size_t nCoeff = xSurface->GetNumberOfRawCoefficients();
  vnl_vector<double> xRotatedCoeff[8], 
    xInitCoeff(xSurface->GetRawCoefficientArray(), nCoeff);

  // Create a volume match object
  VolumeOverlapEnergyTerm tVolumeMatch(image, nCuts);

  // Compute the best image match over all possible flips of the eigenvalues
  // (there are 8 possible matches, including mirroring)
  for(f = 0; f < 8; f++)
    {
    // Set the flip/scale matrix
    vnl_matrix<double> F(3, 3, 0.0);
    F(0,0) = (f & 1) ? -sqrt(Dx(0) / Dy(0)) : sqrt(Dx(0) / Dy(0));
    F(1,1) = (f & 2) ? -sqrt(Dx(1) / Dy(1)) : sqrt(Dx(1) / Dy(1));
    F(2,2) = (f & 4) ? -sqrt(Dx(2) / Dy(2)) : sqrt(Dx(2) / Dy(2));

    // Compute the rotation+flip matrix
    vnl_matrix<double> R = Vx * F * Vy.transpose();

    // Rotate the surface by matrix R and shift to yMean
    xSurface->SetRawCoefficientArray(xInitCoeff.data_block());
    xSurface->ApplyAffineTransform(R, xMean - yMean, yMean);
    xRotatedCoeff[f] = vnl_vector<double>(xSurface->GetRawCoefficientArray(), nCoeff);

    // Compute the boundary
    xSolver->Solve();

    // Create a solution object and a volume match term
    SolutionData SRot(xSolver, false);
    double xMatch = tVolumeMatch.ComputeEnergy(&SRot);
    
    // Record the best match ever
    if(f == 0 || xMatch < xBestMatch)
      { xBestMatch = xMatch; iBestSurface = f; }

    // Report the match quality
    cout << "ROTATION " << f << " :" << endl;
    tVolumeMatch.PrintReport(cout);
    }

  // Use the best surface as the new surface
  xSurface->SetRawCoefficientArray(xRotatedCoeff[iBestSurface].data_block());
  xSolver->Solve();

  // Test the results
  SolutionData S1(xSolver, false);
  S1.UpdateInternalWeights( nCuts );
  for(i = 0; i < 3; i++)
    {
    // Compute the first moment in i-th direction
    FirstMomentComputer fmc(i);
    yMean[i] = IntegrateFunctionOverInterior(
      S1.xAtomGrid, S1.xInternalPoints, S1.xInternalWeights, nCuts, &fmc);

    for(j = 0; j <= i; j++)
      {
      // Compute the second moment in i-th and j-th directions
      SecondMomentComputer smc(i,j);
      yCov[i][j] = yCov[j][i] = IntegrateFunctionOverInterior(
        S1.xAtomGrid, S1.xInternalPoints, S1.xInternalWeights, nCuts, &smc);
      }
    }

  // Compute the actual mean and covariance
  yVolume = S1.xInternalVolume;
  yMean /= yVolume;
  yCov = (yCov - yVolume * outer_product(yMean, yMean)) / yVolume;
  
  cout << "Model Volume: " << yVolume << endl;
  cout << "Model Mean: " << yMean << endl;
  cout << "Model Covariance: " << endl << yCov << endl;
}

void MedialPDE::SaveBYUMesh(const char *file)
{
  // Export the model as a BYU mesh
  fstream fout(file, ios_base::out);

  // Get the atom grid
  MedialAtomGrid *grid = xSolver->GetAtomGrid();
  MedialAtom *atoms = xSolver->GetAtomArray();

  // Write the number of vertices, edges, parts
  unsigned int nVertices = grid->GetNumberOfBoundaryPoints();
  unsigned int nFaces = grid->GetNumberOfBoundaryQuads();
  unsigned int nEdges = 2 * (nFaces + nVertices - 2);
  fout << 1 << "\t" << nVertices << "\t" << nFaces << "\t" << nEdges << "\t" << endl;

  // Write the number of faces in the first part
  fout << 1 << "\t" << nFaces << endl;

  // Write all the vertices
  MedialBoundaryPointIterator *itVertex = grid->NewBoundaryPointIterator();
  while(!itVertex->IsAtEnd())
    {
    fout << GetBoundaryPoint(itVertex, atoms).X << endl;
    ++(*itVertex);
    }
  delete itVertex;

  // Write all the faces
  MedialBoundaryQuadIterator *itQuad = grid->NewBoundaryQuadIterator();
  while(!itQuad->IsAtEnd())
    {
    int i00 = 1 + itQuad->GetBoundaryIndex(0, 0);
    int i01 = 1 + itQuad->GetBoundaryIndex(0, 1);
    int i11 = 1 + itQuad->GetBoundaryIndex(1, 1);
    int i10 = - (1 + itQuad->GetBoundaryIndex(1, 0));

    fout << i00 << "\t" << i01 << "\t" << i11 << "\t" << i10 << endl;

    ++(*itQuad);
    }
  delete itQuad;

  // Close the file
  fout.close();
}

void MedialPDE::SampleImage(FloatImage *fiInput, FloatImage *fiOut, size_t zSamples)
{
  // Allocate the flattened intensity output image
  typedef FloatImage::WrapperType::ImageType ImageType;
  ImageType::Pointer imgOutSmall = ImageType::New();
  ImageType::RegionType regOutput;

  // Get the image dimensions
  size_t m = xSolver->GetNumberOfUPoints();
  size_t n = xSolver->GetNumberOfVPoints();
    
  // Initialize the output region size
  regOutput.SetSize(0, m);
  regOutput.SetSize(1, n);
  regOutput.SetSize(2, zSamples * 2 + 1);
  
  // Continue creating the image
  imgOutSmall->SetRegions(regOutput);
  imgOutSmall->Allocate();
  imgOutSmall->FillBuffer(0.0f);

  // Get the input image
  ImageType::Pointer imgInput = fiInput->xImage->GetInternalImage();

  // Create the output images for U, V, Tau coordinates
  ImageType::Pointer imgCoord[3];
  for(size_t d = 0; d < 3; d++)
    {
    imgCoord[d] = ImageType::New();
    imgCoord[d]->SetRegions(imgInput->GetBufferedRegion());
    imgCoord[d]->SetSpacing(imgInput->GetSpacing());
    imgCoord[d]->SetOrigin(imgInput->GetOrigin());
    imgCoord[d]->Allocate();
    imgCoord[d]->FillBuffer(-2.0);
    }

  // Create a nearest neighbor interpolator
  //typedef itk::NearestNeighborInterpolateImageFunction<
  //  ImageType, double> InterpolatorType;
  typedef itk::BSplineInterpolateImageFunction<
    ImageType, double> InterpolatorType;
  InterpolatorType::Pointer fInterp = InterpolatorType::New();
  fInterp->SetInputImage(imgInput);

  // Get an internal point iterator
  MedialInternalPointIterator *itPoint = 
    xSolver->GetAtomGrid()->NewInternalPointIterator(zSamples - 1);
  for(; !itPoint->IsAtEnd(); ++(*itPoint))
    {
    // Interpolate the position in 3-space
    SMLVec3d xPoint = GetInternalPoint(itPoint, xSolver->GetAtomArray());

    // Get the corresponding medial atom
    MedialAtom &xAtom = xSolver->GetAtomArray()[itPoint->GetAtomIndex()];

    // Get the i and j coordinates of the atom (using cartesian logic)
    size_t i = (size_t) round(xAtom.u * (m - 1));
    size_t j = (size_t) round(xAtom.v * (n - 1));
    size_t k = itPoint->GetBoundarySide() ? 
      zSamples - itPoint->GetDepth() :
       zSamples + itPoint->GetDepth();

    // Create an index into the output image
    ImageType::IndexType idxTarget;
    idxTarget.SetElement(0, i); 
    idxTarget.SetElement(1, j); 
    idxTarget.SetElement(2, k); 

    // Get the point at which to sample the image
    double tau = itPoint->GetDepth() * 1.0 / itPoint->GetMaxDepth();
    SMLVec3d Z = 
      xAtom.X + tau * xAtom.R * xAtom.xBnd[itPoint->GetBoundarySide()].N;
    
    // Create a point and a continuous index
    itk::Point<double, 3> ptZ(Z.data_block());
    itk::ContinuousIndex<double, 3> idxZ;
    imgInput->TransformPhysicalPointToContinuousIndex(ptZ, idxZ);
    float f = fInterp->EvaluateAtContinuousIndex(idxZ);

    // Sample the input image
    // float f = imgInput->xImage->Interpolate(Z[0], Z[1], Z[2], 0.0f);
    imgOutSmall->SetPixel(idxTarget, f); 

    // This is junk!
    SMLVec3d Z1 = Z + 0.5;
    ImageType::IndexType idxCoord;
    itk::Point<double, 3> ptZ1(Z1.data_block());
    imgCoord[0]->TransformPhysicalPointToIndex(ptZ1, idxCoord);
    imgCoord[0]->SetPixel(idxCoord, xAtom.u);
    imgCoord[1]->SetPixel(idxCoord, xAtom.v);
    imgCoord[2]->SetPixel(idxCoord, itPoint->GetBoundarySide() ? tau : -tau);

    //if(rand() % 80 == 0)
    //  {
    //  cout << "Sample " << Z << " value " << f << endl;
    //  }

    // If the index is at the edge, there are two output pixels
    if(itPoint->IsEdgeAtom())
      {
      idxTarget.SetElement(2, zSamples * 2 - k);
      imgOutSmall->SetPixel(idxTarget, f);
      }
    }

  delete itPoint;

  // Store the output image
  fiOut->xImage->SetInternalImage(imgOutSmall);
  fiOut->xGradient[0]->SetInternalImage(imgCoord[0]);
  fiOut->xGradient[1]->SetInternalImage(imgCoord[1]);
  fiOut->xGradient[2]->SetInternalImage(imgCoord[2]);
}

void MedialPDE::SetIntensityImage(FloatImage *imgSource)
{
  // Copy the internal image
  imgIntensity.xImage->SetInternalImage(imgSource->xImage->GetInternalImage());

  // Set the flag
  flagIntensityPresent = true;
}

void MedialPDE::GetIntensityImage(FloatImage *imgTarget)
{
  imgTarget->xImage->SetInternalImage(imgIntensity.xImage->GetInternalImage());
}

class Timer 
{
public:
  void Start()
    { tStart = clock(); }
  
  void Stop()
    { tElapsed = (clock() - tStart) * 1.0 / CLOCKS_PER_SEC; }

  double Read()
    { return tElapsed; }

private:
  clock_t tStart;
  double tElapsed;
};

void MedialPDE::TestDerivativeComputation(unsigned int iComp)
{
  // Timers for different routines
  Timer tCentral, tAnalytic;
  
  // Get the number of coefficients in use
  unsigned int ncu, ncv, nu, nv;
  xSurface->GetNumberOfCoefficientsUV(ncu, ncv);
  nu = xSolver->GetNumberOfUPoints();
  nv = xSolver->GetNumberOfVPoints();
  
  // Create a dummy fourier surface
  FourierSurface *xEta = new FourierSurface(ncu, ncu);
  xEta->SetRawCoefficient( iComp, 1.0 );

  // Create an array of atoms to hold the derivative
  MedialAtom *dAtoms = 
    new MedialAtom[xSolver->GetAtomGrid()->GetNumberOfAtoms()];
  
  // Solve for Phi 
  xSolver->Solve();

  // Compute the variational derivative
  tAnalytic.Start();
  xSolver->ComputeVariationalDerivative(xEta, (iComp % 4) == 3, dAtoms);
  tAnalytic.Stop();
  
  // Now, use central differences to compute the same
  double cc = xSurface->GetRawCoefficient(iComp), eps = 0.001;

  // Create the first problem
  MedialPDESolver S1(nu, nv);
  S1.SetMedialSurface(xSurface);

  // Create the second problem
  MedialPDESolver S2(nu, nv);
  S2.SetMedialSurface(xSurface);
  
  // Solve both problems
  tCentral.Start();
  xSurface->SetRawCoefficient(iComp, cc + eps);
  S1.Solve();

  xSurface->SetRawCoefficient(iComp, cc - eps);
  S2.Solve();
  tCentral.Stop();
  
  // Get the atom arrays
  MedialAtom *A1 = S1.GetAtomArray();
  MedialAtom *A2 = S2.GetAtomArray();

  // Compute the values of dPhi/dC between central difference and analytic
  // derivatives
  unsigned int n = xSolver->GetAtomGrid()->GetNumberOfAtoms();
  double xMaxPhiDiff = 0.0;
  double xMaxBndDiff = 0.0;

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

    // Take the largest difference in boundary derivs
    for(unsigned int k = 0; k < 2; k++)
      {
      SMLVec3d dcBnd = (0.5 / eps) * (a1.xBnd[k].X - a2.xBnd[k].X);
      SMLVec3d daBnd = a0.xBnd[k].X;
      double xBndDiff = (dcBnd - daBnd).two_norm();
      if(xBndDiff > xMaxBndDiff) xMaxBndDiff = xBndDiff;
      }
    }

  // Print out the max difference
  cout << "Max difference between numeric and analytic derivs (Phi) : " 
    <<  xMaxPhiDiff << endl;

  cout << "Max difference between numeric and analytic derivs (Bnd) : " 
    <<  xMaxBndDiff << endl;

  cout << "Central difference computed in " << tCentral.Read() 
    << " sec., Analytic in " << tAnalytic.Read() << " sec." << endl; 

  /* 
  // Take a look at the atoms themselves
  for(unsigned int j = 0; j < n; j++)
    {
    MedialAtom &a1 = S1.GetAtomArray()[j];
    MedialAtom &a2 = S2.GetAtomArray()[j];
    MedialAtom &a0 = dAtoms[j];
    
    SMLVec3d dcX1 = (0.5 / eps) * (a1.xBnd[0].X - a2.xBnd[0].X);
    SMLVec3d daX1 = a0.xBnd[0].X;

    cout << j << " :: " << dcX1 << "  ---  " << daX1 << endl;
    }
  */

  /* 
  // Compute the right hand side for the boundary equation
  for(unsigned int j = 0; j < n; j++)
    {
    double g1 = dot_product(A1[j].xGradR, A1[j].xGradR);
    double g2 = dot_product(A2[j].xGradR, A2[j].xGradR);
    double d2 = 0.5 * (
      g1 * 4.0 * A1[j].R * A1[j].R - 
      g2 * 4.0 * A2[j].R * A2[j].R ) / eps;
    // cout << "Test 2 : " << j << " -- " << d2 << endl;
    
    // SMLVec3d G1 = 2.0 * A1[j].R * A1[j].xGradR;
    // SMLVec3d G2 = 2.0 * A2[j].R * A2[j].xGradR;
    // cout << "Test 2 : " << j << " -- " << 0.5 * (G1 - G2) / eps << endl;
    
    g1 = A1[j].G.xContravariantTensor[0][1];
    g2 = A2[j].G.xContravariantTensor[0][1];

    
    
    // cout << "Test 2 : " << j << " -- " << 0.5 * (g1 - g2) / eps << endl;
    }
    
  double g1 = S1.ComputeLaplaceBeltrami(0, 4, xSolver->GetPhiField());
  double g2 = S2.ComputeLaplaceBeltrami(0, 4, xSolver->GetPhiField());
  cout << "Test 2 : (0, 4) -- " << 0.5 * (g1 - g2) / eps << endl;
  */
}








class PrincipalComponents
{
public:
  typedef vnl_matrix<double> Mat;
  typedef vnl_vector<double> Vec;

  PrincipalComponents(const Mat &A);

  /** Get the number of principal modes */
  size_t GetNumberOfModes() 
    { return nModes; }

  /** Get the eigenvalues in decreasing order */
  double GetEigenvalue(size_t i) const
    { return lambda[i]; }

  /** Get the eigenvectors in decreasing order of eigenvalue */
  Vec GetEigenvector(size_t i) const
    { return Q.get_row(i); }

  /** Get the shape at a certain position in PCA space. The coordinates in this
   * space are the eigenvectors of PCA and the units are standard deviations
   * (sqrt of the eigenvalues). The input vector can be shorter than the number
   * of modes; the remaining modes will be treated as zero */
  Vec MapToFeatureSpace(const Vec &xComponent) const
    {
    Vec x = mu;
    for(size_t i = 0; i < xComponent.size(); i++)
      x += xComponent[i] * sdev[i] * GetEigenvector(i);
    return x;
    }

  /** Project a point in feature (shape) space into the coefficient space */
  Vec MatToCoefficientSpace(const Vec &xFeature, size_t nCoeff) const
    {
    Vec z(m, 0.0);
    if(nCoeff > nModes) nCoeff = nModes;
    for(size_t i = 0; i < nCoeff; i++)
      z[i] = dot_product(xFeature - mu, GetEigenvector(i)) / sdev[i];
    return z;
    }

  /** Get a feature space vector along one of the principal modes. Useful for 
   * animating the PCA modes */
  Vec MapToFeatureSpace(size_t iComponent, double xComponentValue)
    {
    return mu + xComponentValue * sdev[iComponent] * GetEigenvector(iComponent);
    }

  /**
   * Perform leave-one-out analysis, i.e., how well the PCA can predict the 
   * individual members.
   */
  // static void LeaveOneOutAnalysis(const Mat &A);

private:
  // The mean vector
  Vec mu;

  // The principal component vectors (rows)
  Mat Q;

  // The eigenvalues
  Vec lambda, sdev;

  // The dimensions of the data
  size_t m, n, nModes;
};

/*
PrincipalComponents::LeaveOneOutAnalysis(const Mat &A)
{
  // Repeat for each element of A
  unsigned int i, n = A.rows();
  for(i = 0; i < n; i++)
    {
    // Create a matrix with row i taken out
    Mat B(A.rows()-1, A.columns());
    if(i > 0) B.update(0,0,A.get_n_rows(0, i));
    if(i < n-1) B.update(i,0,A.get_n_rows(i+1,n-i-1));

    // Compute the PCA of B
    PrincipalComponents p(B);

    // Project the left out guy on the principal components
    p.

    }
}
*/


PrincipalComponents::PrincipalComponents(const Mat &A)
{
  size_t i, j;
  m = A.columns();
  n = A.rows();

  // Compute the mean
  mu.set_size(m); mu.fill(0.0);
  for(i = 0; i < n; i++) for(j = 0; j < m; j++)
    mu[j] += A(i,j);
  mu /= n;

  // Compute the zero mean matrix
  Mat Z(n, m);
  for(i = 0; i < n; i++) for(j = 0; j < m; j++)
    Z(i,j) = A(i,j) - mu[j];

  // Split, depending on if m > n or not
  if(m > n)
    {
    // Compute the covariance matrix Z Zt
    Mat Zt = Z.transpose();
    Mat K = ( Z * Zt ) / n;

    // Get the eigensystem of K
    vnl_symmetric_eigensystem<double> eig(K);

    // Take the eigenvalues
    lambda.set_size(n); Q.set_size(n, m); 
    lambda.fill(0.0); Q.fill(0.0);
    
    for(i = 0; i < n; i++) 
      {
      lambda[i] = eig.get_eigenvalue(n - 1 - i);
      if(lambda[i] > 0.0)
        Q.set_row(i, (Zt * eig.get_eigenvector(n - 1 - i)).normalize());
      else break;
      }
    }
  else
    {
    // Compute the real covariance matrix At Z
    Mat K = (Z.transpose() * Z) / n;
    vnl_symmetric_eigensystem<double> eig(K);

    // Copy the eigenvalues and vectors
    lambda.set_size(m); Q.set_size(m,m); 
    lambda.fill(0.0); Q.fill(0.0);
    
    for(i = 0; i < m; i++) 
      {
      lambda[i] = eig.get_eigenvalue(m - 1 - i);
      if(lambda[i] > 0.0)
        Q.set_row(i, eig.get_eigenvector(n - 1 - i));
      else break;
      }
    }

  // Finish up
  nModes = i;
  sdev = lambda.apply(sqrt);
}

/**
 * MEDIAL PCA CODE
 */
void MedialPCA::AddSample(MedialPDE *pde)
{
  // Make sure that the number of coefficients matches
  if(xSurfaces.size())
    assert(xSurfaces.front()->GetNumberOfRawCoefficients()
      == pde->xSurface->GetNumberOfRawCoefficients());

  // Get the coefficients from this medial PDE
  FourierSurface *xNew = new FourierSurface(*(pde->xSurface));
  xSurfaces.push_back(xNew);
  
  // Make a copy of the intensity values associated with this PDE
  if(pde->flagIntensityPresent)
    {
    // Make a copy of the image
    FloatImage *img = new FloatImage();
    img->xImage->SetInternalImage(
      pde->imgIntensity.xImage->GetInternalImage());
    xAppearance.push_back(img);
    }
}

MedialPCA::MedialPCA()
{
  xPCA = NULL;
  xAppearancePCA = NULL;
}

MedialPCA::~MedialPCA()
{
  if(xPCA) delete xPCA;
  if(xAppearancePCA) delete xAppearancePCA;
}

void MedialPCA::ComputePCA()
{
  // The non-shape information must be cleaned from the data using
  // something like the generalized Procrustes method. We begin by
  // interpolating each mpde
  typedef vnl_matrix<double> Mat;
  typedef vnl_vector<double> Vec;

  // These are the parameters to the GPA method
  Mat *A = new Mat[xSurfaces.size()];
  Mat *R = new Mat[xSurfaces.size()];
  Vec *t = new Vec[xSurfaces.size()];
  double *s = new double[xSurfaces.size()];
  size_t i, j;
  
  // Populate the input matrices
  MedialPDESolver xSolver(32, 80);
  for(i = 0; i < xSurfaces.size(); i++)
    {
    // Solve for this surface
    xSolver.SetMedialSurface(xSurfaces[i]);
    xSolver.Solve();

    // Initialize the A matrix
    A[i].set_size(xSolver.GetAtomGrid()->GetNumberOfAtoms(), 3);
    R[i].set_size(3,3);
    t[i].set_size(3);

    // Parse over all the boundary sites
    MedialAtomIterator *it = 
      xSolver.GetAtomGrid()->NewAtomIterator();
    size_t j = 0;
    while(!it->IsAtEnd())
      {
      MedialAtom &atom = xSolver.GetAtomArray()[it->GetIndex()];
      A[i][j][0] = atom.X[0];
      A[i][j][1] = atom.X[1];
      A[i][j][2] = atom.X[2];
      ++(*it); ++j;
      }
    delete it;
    }

  // cout << "Testing Procrustes" << endl;
  // TestProcrustes(A[0]);
  
  // Run the procrustes method
  cout << "Procrustenating..." << endl;
  
  GeneralizedProcrustesAnalysis(xSurfaces.size(), A, R, t, s);

  // A medial PCA to play with
  MedialPDE xJunk(8,12,32,80);
    
  // Rotate each of the models into the common coordinate frame
  for(i = 0; i < xSurfaces.size(); i++)
    {
    xSurfaces[i]->ApplyAffineTransform(s[i] * R[i].transpose(), t[i], Vec(3, 0.0));
    
    MedialPDE xJunk(8,12,32,80);
    xJunk.xSurface = new FourierSurface(*xSurfaces[i]);
    xJunk.xSolver->SetMedialSurface(xJunk.xSurface);
    xJunk.xSurface->SetRawCoefficientArray(xSurfaces[i]->GetRawCoefficientArray());
    xJunk.Solve();
    }

  // Create a principal components object

  // Compute the mean shape and the covariance matrix on the fourier
  // parameters. Since the Fourier basis is orthonormal, doing PCA on the
  // surface and on the Fourier components is identical
  size_t m = xSurfaces[0]->GetNumberOfRawCoefficients();
  size_t n = xSurfaces.size();

  // Populate the data matrix 
  xDataShape.set_size(n, m);
  for(i = 0; i < n; i++) for(j = 0; j < m; j++)
    xDataShape[i][j] = xSurfaces[i]->GetRawCoefficient(j);

  // Compute the principal components of xDataShape
  if(xPCA) delete xPCA;
  xPCA = new PrincipalComponents(xDataShape);
  xPCALocation.set_size(xPCA->GetNumberOfModes());
  xPCALocation.fill(0.0);

  // If there is appearance information, also compute the appearance PCA
  if(xAppearance.size() == n)
    {
    // Create a matrix of appearance values
    typedef FloatImage::WrapperType::ImageType ImageType;
    ImageType::Pointer imgFirst = 
      xAppearance.front()->xImage->GetInternalImage();
    size_t mpix = imgFirst->GetBufferedRegion().GetNumberOfPixels();
    xDataAppearance.set_size(n, mpix);

    // Populate the matrix from the image
    for(i = 0; i < n; i++)
      {
      float *pix = xAppearance[i]->xImage->GetInternalImage()->GetBufferPointer();
      for(j = 0; j < mpix; j++)
        xDataAppearance(i, j) = pix[j];
      }

    // Compute the PCA from the image array
    if(xAppearancePCA) delete xAppearancePCA;
    xAppearancePCA = new PrincipalComponents(xDataAppearance);
    xAppearancePCALocation.set_size(xAppearancePCA->GetNumberOfModes());
    xAppearancePCALocation.fill(0.0);
    }

  // Report the principal components
  cout << "EIGENVALUES: " << endl;
  for(size_t l = 0; l < xPCA->GetNumberOfModes(); l++)
    { cout << xPCA->GetEigenvalue(l) << " "; }
  cout << endl;
  
  /*
  // Generate shapes along the first four eigenvectors
  for(i = 0; i < 4; i++)
    {
    for(int y = 0; y < 60; y++)
      {
      // A shape!
      double c = 0.1 * (y - 30);
      Vec z = xPCA->MapToFeatureSpace(i, c);

      cout << "MODE " << i << ", COEFF " << c << ", LAMBDA " << 
        xPCA->GetEigenvalue(i) << endl;

      MedialPDE xJunk(8,12,32,80);
      xJunk.xSurface = new FourierSurface(*xSurfaces[0]);
      xJunk.xSolver->SetMedialSurface(xJunk.xSurface);
      xJunk.xSurface->SetRawCoefficientArray(z.data_block());
      xJunk.Solve();

      ostringstream sJunk1, sJunk2;
      sJunk1 << "/tmp/pcamed_m" << i << "_" << (y / 10) << (y % 10) << ".vtk";
      sJunk2 << "/tmp/pcabnd_m" << i << "_" << (y / 10) << (y % 10) << ".vtk";    
      xJunk.SaveVTKMesh(sJunk1.str().c_str(), sJunk2.str().c_str());
      }
    }
  */
}

// Move along a given mode a certain number of S.D.
void MedialPCA::SetFSLocationToMean()
{
  xPCALocation.fill(0.0);
  if(xAppearancePCALocation.size())
    xAppearancePCALocation.fill(0.0);
}

// Export shape matrix
void MedialPCA::ExportShapeMatrix(const char *filename)
{
  vnl_matlab_filewrite writer(filename);
  writer.write(xDataShape, "xShape");
}


// Move along a given mode a certain number of S.D.
void MedialPCA::SetFSLocation(unsigned int iMode, double xSigma)
{
  xPCALocation(iMode) = xSigma;
  if(xAppearancePCALocation.size())
    xAppearancePCALocation(iMode) = xSigma;
}

// Generate a sample at the current location
void MedialPCA::GetShapeAtFSLocation(MedialPDE *target)
{
  typedef vnl_matrix<double> Mat;
  typedef vnl_vector<double> Vec;

  // Compute the shape space point
  Vec z = xPCA->MapToFeatureSpace(xPCALocation);

  // Compute the shape
  target->xSurface = new FourierSurface(*xSurfaces[0]);
  target->xSolver->SetMedialSurface(target->xSurface);
  target->xSurface->SetRawCoefficientArray(z.data_block());
  target->Solve();

  // Compute the feature image
  if(xAppearance.size() == xSurfaces.size())
    {
    FloatImage *imgJunk = new FloatImage();
    FloatImage::WrapperType::ImageType *img = imgJunk->xImage->GetInternalImage();
    img->SetRegions(
      xAppearance[0]->xImage->GetInternalImage()->GetBufferedRegion());
    img->Allocate();

    Vec z1 = xAppearancePCA->MapToFeatureSpace(xAppearancePCALocation);
    for(size_t i = 0; i < z1.size(); i++)
      img->GetBufferPointer()[i] = z1[i];

    target->SetIntensityImage(imgJunk);
    delete imgJunk;
    }
}

void RenderMedialPDE(MedialPDE *model)
{
  // Create a renderer
  PDESplineRenderer *rndPDESpline = new PDESplineRenderer(model->xSolver);

  // Initialize the GL environment
  char *argv[] = {"test",NULL};
  GLDisplayDriver::init(1, argv);

  // Add the important renderers
  GLDisplayDriver::addRenderer(new DefaultLightRenderer(), GLDisplayDriver::EYE);
  GLDisplayDriver::addRenderer(new StarfieldRenderer(), GLDisplayDriver::UNIT);

  // Add our spline renderer
  GLDisplayDriver::addRenderer(rndPDESpline);

  // Put the display into Talairach coordinates
  GLDisplayDriver::center = SMLVec3f(102,134,73);
  GLDisplayDriver::scale = 20;

  // Start GLUT
  glutMainLoop();
}


} // namespace
