#ifndef __ScriptInterface_h_
#define __ScriptInterface_h_

class MedialPDESolver;
class FourierSurface;
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

  // Save the image to a file
  void SaveToFile(const char *file);

  // Load the image from a file
  void LoadGradientFromFile(unsigned int iComponent, const char *file);

  // Save the image to a file
  void SaveGradientToFile(unsigned int iComponent, const char *file);

  // Compute the image by processing a binary image
  void SetToGradientMagnitude(BinaryImage *imgSource, double xSigma);

  // Compute the image by blurring a binary image. The surface is represented
  // by the 0.0 level set, 1.0 is inside, -1.0 is outside of the object
  void SetToBlurredBinary(BinaryImage *imgSource, double xSigma);

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
  void LoadFromDiscreteMRep(const char *file, unsigned int xResolution);

  /** Save to a parameter file */
  void SaveToParameterFile(const char *file);

  /** Load from a parameter file */
  void LoadFromParameterFile(const char *file);

  /** Some default initialization */
  void GenerateSampleModel();

  /** Compute the match between the model and a floating point image */
  double ComputeImageMatch(FloatImage *image);
  
  /** Perform gradient descent */
  void GradientDescentOptimization(FloatImage *image, unsigned int nSteps, double xStep);
  void ConjugateGradientOptimization(FloatImage *image, unsigned int nSteps);
  void EvolutionaryOptimization(FloatImage *image, unsigned int nSteps);
  
  /** Fit the model to the binary image by matching moments of inertia */
  void MatchImageByMoments(BinaryImage *image);

  /** Save the model as a BYU mesh */
  void SaveBYUMesh(const char *file);

  /** Compute the radius function after surface/pho update */
  void Solve();

private:
  MedialPDESolver *xSolver;
  FourierSurface *xSurface;

  // Friend functions
  friend void RenderMedialPDE(MedialPDE *);
};

/** Function to visualize the medial PDE result */
void RenderMedialPDE(MedialPDE *model);

} // Namespace medial PDE!

#endif
