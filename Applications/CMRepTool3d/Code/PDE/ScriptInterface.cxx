#include "ScriptInterface.h"
#include "MedialPDESolver.h"
#include "BasisFunctions2D.h"
#include "ITKImageWrapper.h"
#include "MedialPDERenderer.h"
#include "OptimizationTerms.h"
#include "Registry.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "ConjugateGradientMethod.h"
#include "EvolutionaryStrategy.h"

#include <iostream>
#include <fstream>

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
MedialPDE::MedialPDE(unsigned int nBasesU, unsigned int nBasesV, unsigned int xResolution)
{
  xSurface = new FourierSurface(nBasesU, nBasesV);
  xSolver = new MedialPDESolver(xResolution * nBasesU + 1, xResolution * nBasesV + 1);
  xSolver->SetMedialSurface(xSurface);
}

MedialPDE::~MedialPDE()
{
  delete xSolver;
  delete xSurface;
}

void MedialPDE::LoadFromDiscreteMRep(const char *file, unsigned int xResolution)
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

  for(unsigned int i = 0; i < xAtoms.size(); i++)
    {
    xx[i] = xAtoms[i].x;
    yy[i] = xAtoms[i].y;
    zz[i] = xAtoms[i].z;
    uu[i] = xAtoms[i].u;
    vv[i] = xAtoms[i].v;
    }

  // Perform the fitting on x, y and z
  xSurface->FitToData(xAtoms.size(), 0, uu, vv, xx);
  xSurface->FitToData(xAtoms.size(), 1, uu, vv, yy);
  xSurface->FitToData(xAtoms.size(), 2, uu, vv, zz);

  // Clean up
  delete xx; delete yy; delete zz; delete uu; delete vv;
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

      // Create a new solver and solve the PDE
      // delete xSolver; xSolver = new MedialPDESolver(m, n);
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

void MedialPDE::GenerateSampleModel()
{
  // Decide how many points to interpolate
  unsigned int nSide = 21, nPoints = nSide * nSide;
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
  xSurface->FitToData(nPoints, 0, uPoints, vPoints, xPoints);
  xSurface->FitToData(nPoints, 1, uPoints, vPoints, yPoints);
  xSurface->FitToData(nPoints, 2, uPoints, vPoints, zPoints);
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
  return termImage.ComputeEnergy(&S);
}
  
void MedialPDE
::GradientDescentOptimization(FloatImage *image, unsigned int nSteps, double xStep)
{
  // Create the optimization problem
  MedialOptimizationProblem xProblem(xSolver, xSurface);

  // Create an image match term and a jacobian term
  BoundaryImageMatchTerm xTermImage(image);
  BoundaryJacobianEnergyTerm xTermJacobian;

  // Add the terms to the problem
  xProblem.AddEnergyTerm(&xTermImage, 1.0);
  xProblem.AddEnergyTerm(&xTermJacobian, 0.005);

  // Create the initial solution
  size_t nCoeff = xSurface->GetNumberOfRawCoefficients();
  vnl_vector<double> xSolution(xSurface->GetRawCoefficientArray(), nCoeff);
  vnl_vector<double> xGradient(nCoeff, 0.0);
  
  // Iterate, moving in the gradient direction
  for(unsigned int p = 0; p < nSteps; p++)
    {
    // Compute the gradient and the image match
    double xMatch = 
      xProblem.ComputeGradient( xSolution.data_block(), xGradient.data_block());

    // Move in the gradient direction
    xSolution -= xStep * xGradient;

    // Report the current state
    cout << "Iteration " << p << "; Match = " << xMatch << endl;
    }

  // Store the last solution
  xSurface->SetRawCoefficientArray(xSolution.data_block());
  xSolver->Solve();
}


void MedialPDE
::ConjugateGradientOptimization(FloatImage *image, unsigned int nSteps)
{
  // Create the optimization problem
  MedialOptimizationProblem xProblem(xSolver, xSurface);

  // Create an image match term and a jacobian term
  BoundaryImageMatchTerm xTermImage(image);
  BoundaryJacobianEnergyTerm xTermJacobian;

  // Add the terms to the problem
  xProblem.AddEnergyTerm(&xTermImage, 1.0);
  xProblem.AddEnergyTerm(&xTermJacobian, 0.005);

  // Create the initial solution
  size_t nCoeff = xSurface->GetNumberOfRawCoefficients();
  vnl_vector<double> xSolution(xSurface->GetRawCoefficientArray(), nCoeff);
  
  // Construct the conjugate gradient optimizer
  ConjugateGradientMethod xMethod(xProblem, Vector(nCoeff, xSolution.data_block()));
  for(size_t p = 0; p < nSteps; p++)
    {
    if(xMethod.isFinished())
      break;
    xMethod.performIteration();
    }

  // Store the best result
  xSurface->SetRawCoefficientArray(xMethod.getBestEverX().getDataArray());
  xSolver->Solve();
}

inline vnl_vector<double> VectorCast(const pauly::Vector &Y)
{
  return vnl_vector<double>(Y.getDataArray(), Y.size());
}

inline pauly::Vector VectorCast(vnl_vector<double> &Y)
{
  return pauly::Vector(Y.size(), Y.data_block());
}

void MedialPDE
::EvolutionaryOptimization(FloatImage *image, unsigned int nSteps)
{
  // Create the optimization problem
  MedialOptimizationProblem xProblem(xSolver, xSurface);

  // Create an image match term and a jacobian term
  BoundaryImageMatchTerm xTermImage(image);
  BoundaryJacobianEnergyTerm xTermJacobian;

  // Add the terms to the problem
  xProblem.AddEnergyTerm(&xTermImage, 1.0);
  xProblem.AddEnergyTerm(&xTermJacobian, 0.005);

  // Create the initial solution
  size_t nCoeff = xSurface->GetNumberOfRawCoefficients();
  vnl_vector<double> xSolution(xSurface->GetRawCoefficientArray(), nCoeff);
  vnl_vector<double> xSigma(nCoeff, 0.01);
  
  // Create the initial search space
  GaussianSS xSearch(VectorCast(xSolution), VectorCast(xSigma));
  
  // Construct the evolutionary optimizer
  EvolutionaryStrategy xMethod(xProblem, xSearch, 2, 4, SELECTION_MuPlusLambda);

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
  xSurface->SetRawCoefficientArray(xMethod.getBestEverX().getDataArray());
  xSolver->Solve();
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

void MedialPDE::MatchImageByMoments(BinaryImage *image)
{
  // Compute the mean and covariance matrix of the non-zero voxels
  typedef vnl_matrix_fixed<double,3,3> Mat;
  typedef SMLVec3d Vec;

  Mat xCov(0.0f);
  Vec xMean(0.0f);
  int n = 0;

  // Iterate over all voxels in the image
  typedef BinaryImage::WrapperType::ImageType ImageType;
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;
  ImageType::Pointer xImage = image->xImage->GetInternalImage();
  IteratorType it(xImage, xImage->GetBufferedRegion());

  while(!it.IsAtEnd())
    {
    if(it.Get() != 0)
      {
      // Get the spatial position of the point
      itk::ContinuousIndex<double, 3> iRaw(it.GetIndex());
      itk::Point<double, 3> ptPosition;
      xImage->TransformContinuousIndexToPhysicalPoint(iRaw, ptPosition);
      SMLVec3d xPosition(ptPosition.GetDataPointer());

      // Add to the mean and 'covariance'
      xMean += xPosition;
      xCov += outer_product(xPosition, xPosition);
      n++;
      }
    ++it;
    }

  // Compute the actual statistics
  double xVolume = n * xImage->GetSpacing()[0] 
    * xImage->GetSpacing()[1] * xImage->GetSpacing()[2];
  double nn = n;
  xMean = xMean / nn;
  xCov = (xCov - nn *  outer_product(xMean, xMean)) / nn;
  
  cout << "Image Volume: " << xVolume << endl;
  cout << "Image Mean: " << xMean << endl;
  cout << "Image Covariance: " << endl << xCov << endl;

  // Now compute the same statistics for the m-rep (compute coordinates)
  Mat yCov; Vec yMean; double yVolume;
  
  unsigned int i, j, k, f;
  for(i = 0; i < 3; i++)
    {
    for(j = 0; j <= i; j++)
      {
      SecondMomentComputer smc(i,j);
      yCov[i][j] = yCov[j][i] = xSolver->IntegrateVolumeMeasure(&smc, 10, yVolume);
      }
    FirstMomentComputer fmc(i);
    yMean[i] = xSolver->IntegrateVolumeMeasure(&fmc, 10, yVolume);
    }

  // Compute the actual mean and covariance
  yMean /= yVolume;
  yCov = (yCov - yVolume * outer_product(yMean, yMean)) / yVolume;
  
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
  unsigned int iBestSurface;
  double xBestMatch;

  // Get the numbers of coefficients
  unsigned int nCoeff = xSurface->GetNumberOfRawCoefficients();
  vnl_vector<double> xRotatedCoeff[8], 
    xInitCoeff(xSurface->GetRawCoefficientArray(), nCoeff);

  // Compute the best image match over all possible flips of the eigenvalues
  // (there are 8 possible matches, including mirroring)
  for(f = 0; f<8; f++)
    {
    // Set the flip/scale matrix
    vnl_matrix<double> F(3, 3, 0.0), S(3, 3, 0.0);
    F(0,0) = (f & 1) ? -s : s;
    F(1,1) = (f & 2) ? -s : s;
    F(2,2) = (f & 4) ? -s : s;

    // Compute the rotation+flip matrix
    vnl_matrix<double> R = Vx * F * Vy.transpose();

    // Create a temporary problem
    xSurface->SetRawCoefficientArray(xInitCoeff.data_block());
    xSurface->ApplyAffineTransform(R, xMean, yMean);
    xRotatedCoeff[f] = vnl_vector<double>(xSurface->GetRawCoefficientArray(), nCoeff);

    // Solve the problem
    xSolver->Solve();

    // Integrate the image match over the interior
    ImageVolumeMatch<unsigned char> ivm(image->xImage);
    double xMatch = xSolver->IntegrateVolumeMeasure(&ivm, 10, xVolume);

    if(f == 0 || xMatch > xBestMatch)
      {
      xBestMatch = xMatch;
      iBestSurface = f;
      }

    cout << "Rotation " << f << " match value " << (xMatch/xVolume) << endl;
    }

  // Use the best surface as the new surface
  xSurface->SetRawCoefficientArray(xRotatedCoeff[iBestSurface].data_block());
  xSolver->Solve();
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
