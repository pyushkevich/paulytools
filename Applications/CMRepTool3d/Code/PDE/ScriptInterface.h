#ifndef __ScriptInterface_h_
#define __ScriptInterface_h_

class MedialPDESolver;
class FourierSurface;
class PDEFourierWrapper;
class ImageCubeITKFloat;

/** This is a high level representation of a 3D non-binary image */
class Image3f
{
public:
  Image3f();
  ~Image3f();

  // Load the image from a file
  void LoadFromFile(const char *file);

private:
  // Internally stored image
  ImageCubeITKFloat *xImageCube;

  // Allow the medial PDE to access this class
  friend class MedialPDE;
};

/** This is the highest level class for working with the PDE application */
class MedialPDE
{
public:
  MedialPDE(unsigned int nBasesU, unsigned int nBasesV, unsigned int xResolution);
  ~MedialPDE();
  
  void LoadFromDiscreteMRep(const char *file, unsigned int xResolution);
  void GenerateSampleModel();
  float ComputeImageMatch(Image3f *image);
  void Solve();

private:
  MedialPDESolver *xSolver;
  PDEFourierWrapper *xModel;
  FourierSurface *xSurface;
};

#endif
