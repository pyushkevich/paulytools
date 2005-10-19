#ifndef __ScriptInterface_h_
#define __ScriptInterface_h_

#include <string>
#include <vector>
#include <smlmath.h>

class MedialPDESolver;
class FourierSurface;
class MedialOptimizationProblem;
class IMedialCoefficientMask;
class PrincipalComponents;
class PCACoefficientMask;
template <typename TPixel> class ITKImageWrapper;

void TestTriangleAreaPartialDerivative();
void TestTetrahedronVolumePartialDerivative();
void TestFDStuff();

namespace medialpde {

class BinaryImage;
class FloatImage;

/** This is a high level representation of a 3D non-binary image */
class FloatImage
{
public:

  FloatImage();
  ~FloatImage();

  // Load the image from a file
  void LoadFromFile(const char *file);

  // Load the image and the gradient from a path
  void LoadFromPath(const char *filebase, const char *ext);

  // Save the image to a file
  void SaveToFile(const char *file);

  // Load the image and the gradient from a path
  void SaveToPath(const char *filebase, const char *ext);

  // Load the image from a file
  void LoadGradientFromFile(unsigned int iComponent, const char *file);

  // Save the image to a file
  void SaveGradientToFile(unsigned int iComponent, const char *file);

  // Compute the image by processing a binary image
  void SetToGradientMagnitude(BinaryImage *imgSource, double xSigma);

  // Compute the image by blurring a binary image. The surface is represented
  // by the 0.0 level set, 1.0 is inside, -1.0 is outside of the object. The
  // optional error scale parameter applies the erf() function to the image,
  // to bring the values nearer to -1.0 or 1.0. This is useful for volume
  // overlap measures
  void SetToBlurredBinary(BinaryImage *imgSource, double xSigma, double xErrorScale = 0.0);

  // Compute the volume of the interior of the image. The image must represent
  // an follows: internal points near 1.0, external near -1.0, interface 0.0
  double ComputeObjectVolume();

  // Compute the integral of positive voxels in the image. This is the upper
  // bound for the probabilistic image match function that we use. In other
  // words, if the cm-rep model and the object in an image are perfectly
  // matched, the integral of the image inside the model will equal the output
  // of this function
  virtual double IntegratePositiveVoxels();

  // Interpolate the image at a given position
  virtual float Interpolate(const SMLVec3d &x);

  // Interpolate the image gradient
  virtual void InterpolateImageGradient(const SMLVec3d &x, SMLVec3f &g);

  // Interpolate the image gradient
  virtual void InterpolateImageGradient(const SMLVec3d &x, SMLVec3d &g);

  // Check whether the gradient information is available
  bool IsGradientAvailable();

  // Set the default 'outside' value (value of pixels outside of the image)
  void SetOutsideValue(float xValue);

  // Set the origin of the image, since some readers (Analyze) don't set it correctly
  // Do this after loading the image!
  void SetImageOrigin(double ox, double oy, double oz);

private:
  // Internally stored image
  typedef ITKImageWrapper<float> WrapperType;
  WrapperType *xImage;

  // Optional gradient images - used to compute image match with gradient
  WrapperType *xGradient[3];

  // The outside value
  float xOutsideValue;

  // Allow the medial PDE to access this class
  friend class MedialPDE;
  friend class MedialPCA;
  friend class FloatImageEuclideanFunctionAdapter;
};

/** This is a binary image representation. */
class BinaryImage
{
public:
  BinaryImage();
  ~BinaryImage();
  
  // Load the image from a file
  void LoadFromFile(const char *file);

private:
  // Internally stored image
  typedef ITKImageWrapper<unsigned char> WrapperType;
  WrapperType *xImage;

  // MedialPDE has access to private members
  friend class MedialPDE;
  friend class FloatImage;
};

/** This is the highest level class for working with the PDE application */
class MedialPDE
{
public:
  /** Create a new instance of the medial PDE object with a given number of 
   * surface coefficients and a given number of samples */
  MedialPDE(unsigned int nBasesU, unsigned int nBasesV, 
    unsigned int xResU, unsigned int xResV, 
    double xFineScale = 0.0, 
    unsigned int xFineU = 0, unsigned int xFineV = 0);

  ~MedialPDE();
  
  /** Initialize the surface using a discrete m-rep */
  void LoadFromDiscreteMRep(const char *file, double xInitRho);

  /** Save to a parameter file */
  void SaveToParameterFile(const char *file);

  /** Load from a parameter file */
  bool LoadFromParameterFile(const char *file);

  /** Set the number of fourier coefficients */
  void SetNumberOfCoefficients(unsigned int m, unsigned int n);

  /** Some default initialization */
  void GenerateSampleModel();

  /** Compute the match between the model and a floating point image */
  double ComputeImageMatch(FloatImage *image);

  /** Compute the boundary jacobian penalty */
  double ComputeBoundaryJacobianPenalty(bool verbose = false);
  
  /** Fit the model to the binary image by matching moments of inertia */
  void MatchImageByMoments(FloatImage *image, unsigned int nCuts);

  /** Perform gradient descent */
  void RunOptimization(FloatImage *image, unsigned int nSteps);
  
  /** Save the model as a BYU mesh */
  void SaveBYUMesh(const char *file);

  /** Save to a VTK file */
  void SaveVTKMesh(const char *fileMedial, const char *fileBoundary);

  /** Compute the radius function after surface/pho update */
  bool Solve();

  /** Set the optimization mode to affine */
  void SetOptimizationToAffine()
    { eMask = AFFINE; }

  /** 
   * Set optimization to a subset of the coefficients. The first ncu, ncv 
   * coefficients in X and Rho will be used. This makes sense for Fourier
   * style masks
   */
  void SetOptimizationToCoarseToFine(
    int ncuX, int ncvX, int ncuRho, int ncvRho)
    { 
    eMask = COARSE_TO_FINE;
    xCoarseFineParams[0] = ncuX;
    xCoarseFineParams[1] = ncvX;
    xCoarseFineParams[2] = ncuRho;
    xCoarseFineParams[3] = ncvRho;
    }

  /** Set the optimization mode to deformable */
  void SetOptimizationToDeformable()
    { eMask = FULL; }

  /** 
   * Set optimization mode to use the PCA matrix. The matrix file must be
   * specified separately
   */
  void SetOptimizationToPCA(size_t nModes)
    { eMask = PCA; nPCAModes = nModes; }

  /** Set the optimization mode */
  void SetOptimizerToGradientDescent(double xStep)
    { eOptimizer = GRADIENT; xStepSize = xStep; }
  
  void SetOptimizerToConjugateGradientDescent(double xStep)
    { eOptimizer = CONJGRAD; xStepSize = xStep; }

  void SetOptimizerToEvolutionaryMethod(double xStep)
    { eOptimizer = EVOLUTION; }

  /** Set the match type */
  void SetMatchToVolumeOverlap()
    { eMatch = VOLUME; }

  void SetMatchToBoundaryGradient()
    { eMatch = BOUNDARY; }

  /** Set the optimizer dump path. As the optimizer runs, it will periodically
   * dump meshes into path.number.med.vtk and path.number.bnd.vtk */
  void EnableMeshDump(const char *path, double xImprovement = 0.01)
    { strDumpPath = path; xMeshDumpImprovementPercentage = xImprovement; }

  void DisableMeshDump()
    { xMeshDumpImprovementPercentage = 0.0; }

  /** Enable the optimizer dump */
  void EnableOptimizerDump(const char *file)
    { strOptimizerDump = file; }

  void DisableOptimizerDump() 
    { strOptimizerDump = ""; }

  /** Should not be here! */
  MedialPDESolver *GetSolver() { return xSolver; }

  // Get a pointer to the surface
  FourierSurface *GetSurface() { return xSurface; }

  /** This method takes an image and an m-rep and samples the image using the 
   * m-reps coordinate system, with linear interpolation. The result is stored in 
   * a separate image. */
  void SampleImage(FloatImage *imgInput, FloatImage *imgOutput, size_t zSamples);

  /** This method is the opposite of SampleImage. Given an image in the medial coordinate
   * system (it must match the sampling rate of the cm-rep), this with map this image back
   * into the ambient space */
  void SampleReferenceFrameImage(FloatImage *imgInput, FloatImage *imgOutput, size_t zSamples);

  /** Assign an intensity image to the cm-rep */
  void SetIntensityImage(FloatImage *imgIntensity);

  /** Get an intensity image from the cm-rep */
  void GetIntensityImage(FloatImage *imgIntensity);

  /** Set the PCA data for PCA-based optimization */
  void SetPCAMatrix(size_t ncu, size_t ncv, const char *fnMatrix);

  /** Create a PCACoefficientMask based on this medial PDE and the previously
   * specified coefficient matrix */
  IMedialCoefficientMask *CreatePCACoefficientMask(size_t nModes);

  /** Release resources from a coefficient mask */
  void ReleasePCACoefficientMask(IMedialCoefficientMask *xMask);

private:
  // Optimization modes and optimizers
  enum OptimizerType { CONJGRAD, GRADIENT, EVOLUTION };
  enum MaskType { AFFINE, COARSE_TO_FINE, FULL, PCA }; 
  enum MatchType { VOLUME, BOUNDARY };
  
  // The solver
  MedialPDESolver *xSolver;

  // The surface
  FourierSurface *xSurface;

  // Parameters for coarse-to-fine optimization
  int xCoarseFineParams[4];

  // Properties associated with different modes
  double xStepSize;

  // A file where the mesh info is dumped
  std::string strDumpPath, strOptimizerDump;
  double xMeshDumpImprovementPercentage;

  // Current modes
  OptimizerType eOptimizer;
  MaskType eMask;
  MatchType eMatch;

  // An image containing grayscale intensities for this m-rep
  FloatImage imgIntensity;

  // Whether the intensities have been loaded
  bool flagIntensityPresent;

  // Friend functions
  void ExportIterationToVTK(unsigned int iIter);
  void ConjugateGradientOptimization(MedialOptimizationProblem *xProblem, 
    vnl_vector<double> &xSolution, unsigned int nSteps, double xStep);
  friend void RenderMedialPDE(MedialPDE *);
  friend class MedialPCA;

  // PCA matrix (used to create PCA coefficient masks)
  vnl_matrix<double> mPCAMatrix;

  // Number of coefficients used in the PCA matrix
  size_t ncuPCA, ncvPCA, nPCAModes;
};

class MedialPCA
{
public:
  MedialPCA();
  ~MedialPCA();
  
  // Add a sample to the PCA
  void AddSample(MedialPDE *pde);
  
  // Compute the PCA
  void ComputePCA();
  
  // Move current sample location to the mean
  void SetFSLocationToMean();
  
  // Move along a given mode a certain number of S.D.
  void SetFSLocation(unsigned int iMode, double xSigma);
  
  // Generate a sample at the current location
  void GetShapeAtFSLocation(MedialPDE *target);

  // Perform leave-one-out analysis using PCA
  // void LeaveOneOutAnalysis();

  // Export the PCA data for use in Mathematica, etc
  void ExportShapeMatrix(const char *filename);

private:
  std::vector< FloatImage* > xAppearance;
  std::vector< FourierSurface* > xSurfaces;
  vnl_vector<double> xPCALocation;  
  vnl_vector<double> xAppearancePCALocation;  

  vnl_matrix<double> xDataShape, xDataAppearance;

  PrincipalComponents *xPCA, *xAppearancePCA;
};

/** Function to visualize the medial PDE result */
void RenderMedialPDE(MedialPDE *model);

} // Namespace medial PDE!

#endif
