#ifndef __ScriptInterface_h_
#define __ScriptInterface_h_

class MedialPDESolver;
class FourierSurface;
class MedialOptimizationProblem;
class IMedialCoefficientMask;
template <typename TPixel> class ITKImageWrapper;

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
  // by the 0.0 level set, 1.0 is inside, -1.0 is outside of the object
  void SetToBlurredBinary(BinaryImage *imgSource, double xSigma);

  // Compute the volume of the interior of the image. The image must represent
  // an follows: internal points near 1.0, external near -1.0, interface 0.0
  double ComputeObjectVolume();

  // Check whether the gradient information is available
  bool IsGradientAvailable();

  // Set the default 'outside' value (value of pixels outside of the image)
  void SetOutsideValue(float xValue);

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
  MedialPDE(unsigned int nBasesU, unsigned int nBasesV, unsigned int xResolution);
  ~MedialPDE();
  
  /** Initialize the surface using a discrete m-rep */
  void LoadFromDiscreteMRep(
    const char *file, double xInitRho, unsigned int xResolution);

  /** Save to a parameter file */
  void SaveToParameterFile(const char *file);

  /** Load from a parameter file */
  void LoadFromParameterFile(const char *file);

  /** Some default initialization */
  void GenerateSampleModel();

  /** Compute the match between the model and a floating point image */
  double ComputeImageMatch(FloatImage *image);
  
  /** Fit the model to the binary image by matching moments of inertia */
  void MatchImageByMoments(FloatImage *image, unsigned int nCuts);

  /** Perform gradient descent */
  void RunOptimization(FloatImage *image, unsigned int nSteps);
  
  /** Save the model as a BYU mesh */
  void SaveBYUMesh(const char *file);

  /** Compute the radius function after surface/pho update */
  void Solve();

  /** Set the optimization mode to affine */
  void SetOptimizationToAffine()
    { eMask = AFFINE; }

  /** Set the optimization mode to deformable */
  void SetOptimizationToDeformable(double crsX, double crsRho)
    { eMask = COARSE_TO_FINE; xCoarsenessX = crsX; xCoarsenessRho = crsRho; }

  void SetOptimizationToDeformable()
    { eMask = FULL; }

  /** Set the optimization mode */
  void SetOptimizerToGradientDescent(double xStep)
    { eOptimizer = GRADIENT; }
  
  void SetOptimizerToConjugateGradientDescent(double xStep)
    { eOptimizer = CONJGRAD; }

  void SetOptimizerToEvolutionaryMethod(double xStep)
    { eOptimizer = EVOLUTION; }

  /** Should not be here! */
  MedialPDESolver *GetSolver() { return xSolver; }

private:
  // Optimization modes and optimizers
  enum OptimizerType { CONJGRAD, GRADIENT, EVOLUTION };
  enum MaskType { AFFINE, COARSE_TO_FINE, FULL }; 
  
  // The solver
  MedialPDESolver *xSolver;

  // The surface
  FourierSurface *xSurface;

  // Properties associated with different modes
  double xCoarsenessX, xCoarsenessRho, xStepSize;

  // Current modes
  OptimizerType eOptimizer;
  MaskType eMask;

  // Friend functions
  friend void RenderMedialPDE(MedialPDE *);
};

/** Function to visualize the medial PDE result */
void RenderMedialPDE(MedialPDE *model);

} // Namespace medial PDE!

#endif