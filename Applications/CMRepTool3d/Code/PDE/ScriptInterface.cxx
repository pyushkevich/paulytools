#include "ScriptInterface.h"
#include "MedialPDESolver.h"
#include "FourierSurface.h"
#include "ImageCubeITK.h"
#include "itkImageFileReader.h"
#include "itkLinearInterpolateImageFunction.h"

#include <iostream>
#include <fstream>

using namespace std;

// Just extend the image cube to avoid templates
class ImageCubeITKFloat : public ImageCubeITK<float>
{
};

struct DiscreteAtom 
{ 
  double u, v, x, y, z; 
  unsigned int iu, iv; 
};

/** Write a match class for matching our solutions to images */
class MRepImageMatcher : virtual public BoundaryMeasure
{
public:
  // Compute the match between an m-rep and image
  virtual float ComputeMatch(IMedialSurfacePatch *xPatch)
    {
    float xMatch, xArea;
    xMatch = xPatch->IntegrateBoundaryMeasure(this, xArea);
    return xMatch / xArea;
    }
};

class MRepEdgeMatcher : public MRepImageMatcher
{
public:
  // Initialize with an image
  void SetEdgeImage(AbstractImage3D *xImage)
    { this->xImage = xImage; }

  // Compute the match
  float computeBoundaryMeasure(const MedialPoint &mp, int d)
    {
    const float *x = mp.bp[d].X.data_block();
    return xImage->interpolateVoxel(x[0], x[1], x[2]);
    }
  
  float computeCrestBoundaryMeasure(const MedialPoint &mp)
    { return computeBoundaryMeasure(mp, 0); }
  
private:
  AbstractImage3D *xImage;
};

class PDEFourierWrapper : virtual public IMedialPDEProblem
{
public:
  PDEFourierWrapper(FourierSurface *xSurface)
    {
    this->xSurface = xSurface;
    }

  void ComputeJet2(MedialPoint *mp)
    {
    float u = mp->u, v = mp->v;

    for(unsigned int d = 0; d < 3; d++)
      {
      mp->F[d]   = xSurface->Evaluate(u, v, 0, 0, d);
      mp->Fu[d]  = xSurface->Evaluate(u, v, 1, 0, d);
      mp->Fv[d]  = xSurface->Evaluate(u, v, 0, 1, d);
      mp->Fuu[d] = xSurface->Evaluate(u, v, 2, 0, d);
      mp->Fuv[d] = xSurface->Evaluate(u, v, 1, 1, d);
      mp->Fvv[d] = xSurface->Evaluate(u, v, 0, 2, d);
      }

    // Compute the normal
    SMLVec3f Nraw = cross_3d(mp->Xu(), mp->Xv());
    SMLVec3f N = Nraw / Nraw.magnitude();
    memcpy(mp->F3raw.data_block(), Nraw.data_block(), 3 * sizeof(float));
    memcpy(mp->F3.data_block(), N.data_block(), 3 * sizeof(float));
    }

  double ComputeLaplacian(double u, double v)
    {
    double ramp = (0.5 - abs(0.5 - u)) * (0.5 - abs(0.5 - v));
    double test = rand() * 0.5 / RAND_MAX;
    return -1.0; 
    }
  
private:
  FourierSurface *xSurface;
};

MedialPDE::MedialPDE(unsigned int nBasesU, unsigned int nBasesV, unsigned int xResolution)
{
  xSurface = new FourierSurface(nBasesU, nBasesV, 3);
  xModel = new PDEFourierWrapper(xSurface);
  xSolver = new MedialPDESolver(xResolution * nBasesU + 1, xResolution * nBasesV + 1);
}

MedialPDE::~MedialPDE()
{
  delete xSolver;
  delete xSurface;
  delete xModel;
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
  float *xx = new float[xAtoms.size()];
  float *yy = new float[xAtoms.size()];
  float *zz = new float[xAtoms.size()];
  float *uu = new float[xAtoms.size()];
  float *vv = new float[xAtoms.size()];

  for(unsigned int i = 0; i < xAtoms.size(); i++)
    {
    xx[i] = xAtoms[i].x;
    yy[i] = xAtoms[i].y;
    zz[i] = xAtoms[i].z;
    uu[i] = xAtoms[i].u;
    vv[i] = xAtoms[i].v;
    }

  // Perform the fitting on x, y and z
  unsigned int stride = sizeof(DiscreteAtom);
  xSurface->FitData(0, xAtoms.size(), uu, 1, vv, 1, xx, 1, 5, 5);
  xSurface->FitData(1, xAtoms.size(), uu, 1, vv, 1, yy, 1, 5, 5);
  xSurface->FitData(2, xAtoms.size(), uu, 1, vv, 1, zz, 1, 5, 5);

  // Clean up
  delete xx; delete yy; delete zz; delete uu; delete vv;
}

void MedialPDE::Solve()
{
  xSolver->Solve(xModel);
}

void MedialPDE::GenerateSampleModel()
{
  // Decide how many points to interpolate
  unsigned int nSide = 11, nPoints = nSide * nSide;
  double uStep = 1.0 / (nSide - 1);

  // Allocate arrays of points and coordinates
  float xPoints[nPoints], yPoints[nPoints], zPoints[nPoints];
  float uPoints[nPoints], vPoints[nPoints];

  // Create an array of points
  unsigned int i = 0;
  for(unsigned int u = 0; u < nSide; u++)
    for(unsigned int v = 0; v < nSide; v++)
      {
      float uu = u * uStep, vv = v * uStep;

      uPoints[i] = uu;
      vPoints[i] = vv;
      xPoints[i] = 0.5 * uu + 0.25;
      yPoints[i] = vv;
      zPoints[i] = ((uu - 0.5) * (uu - 0.5) + (vv - 0.5) * (vv - 0.5)) * 0.25;
      ++i;
      }

  // Peform the fit
  xSurface->FitData(0, nPoints, uPoints, 1, vPoints, 1, xPoints, 1, 4, 4);
  xSurface->FitData(1, nPoints, uPoints, 1, vPoints, 1, yPoints, 1, 4, 4);
  xSurface->FitData(2, nPoints, uPoints, 1, vPoints, 1, zPoints, 1, 4, 4);
}

float MedialPDE::ComputeImageMatch(Image3f *image)
{
  MRepEdgeMatcher matcher;
  matcher.SetEdgeImage(image->xImageCube);
  return matcher.ComputeMatch(xSolver);
}


Image3f::Image3f()
{
  xImageCube = new ImageCubeITKFloat;
}
  
Image3f::~Image3f()
{
  delete xImageCube;
}

void Image3f::LoadFromFile(const char *file)
{
  // Load an image from ITK
  typedef itk::Image<float, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(file);
  fltReader->Update();
  ImageType::Pointer imgInput = fltReader->GetOutput();

  // Create an image cube from this image
  xImageCube->SetImage(imgInput, 1.0);

  // Output
  cout << "Loaded image of dimensions " << imgInput->GetBufferedRegion() << endl;
}
