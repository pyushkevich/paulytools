#include "MedialPDESolver.h"
#include "FourierSurface.h"
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "mspline.h"
#include "glengine.h"
#include "vnl/vnl_cross.h"

// #include "imaging.txx"
#include "ImageCubeITK.h"

using namespace std;

// Define an arbitrary surface to test the differential geometry
// computation, and later to test the medial PDEs
class MyProblem1 : virtual public IMedialPDEProblem {
public:
  void ComputeJet2(MedialPoint *mp)
    {
    float u = mp->u;
    float v = mp->v;
    mp->F[0] = 0.3*u;
    mp->F[1] = v;
    mp->F[2] = 0.03*(-3 - 3*u + u*u + 4*v + 2*u*v - v*v);
    mp->Fu[0] = 0.3;
    mp->Fu[1] = 0;
    mp->Fu[2] = 0.03*(-3 + 2*u + 2*v);
    mp->Fv[0] = 0;
    mp->Fv[1] = 1;
    mp->Fv[2] = 0.03*(4 + 2*u - 2*v);
    mp->Fuu[0] = 0;
    mp->Fuu[1] = 0;
    mp->Fuu[2] = 0.06;
    mp->Fuv[0] = 0;
    mp->Fuv[1] = 0;
    mp->Fuv[2] = 0.06;
    mp->Fvv[0] = 0;
    mp->Fvv[1] = 0;
    mp->Fvv[2] = -0.06;

    // Compute the normal vector
    mp->F3raw[0] = mp->Fu[1] * mp->Fv[2] - mp->Fu[2] * mp->Fv[1];
    mp->F3raw[1] = mp->Fu[2] * mp->Fv[0] - mp->Fu[0] * mp->Fv[2];
    mp->F3raw[2] = mp->Fu[0] * mp->Fv[1] - mp->Fu[1] * mp->Fv[0];

    double zNorm = sqrt(1.0 / 
      (mp->F3raw[0]*mp->F3raw[0] + 
       mp->F3raw[1]*mp->F3raw[1] + 
       mp->F3raw[2]*mp->F3raw[2]));

    mp->F3[0] = zNorm * mp->F3raw[0];
    mp->F3[1] = zNorm * mp->F3raw[1];
    mp->F3[2] = zNorm * mp->F3raw[2];
    }

  double ComputeLaplacian(double u, double v)
    { return -1.0; }
};

// Define an arbitrary surface to test the differential geometry
// computation, and later to test the medial PDEs
class MyProblem2 : virtual public IMedialPDEProblem {
public:
  void ComputeJet2(MedialPoint *mp)
    {
    float u = mp->u;
    float v = mp->v;
    mp->F[0] = u; mp->F[1] = v; mp->F[2] = 0.0;
    mp->Fu[0] = 1; mp->Fu[1] = 0; mp->Fu[2] = 0;
    mp->Fv[0] = 0; mp->Fv[1] = 1; mp->Fv[2] = 0;
    mp->Fuu[0] = 0; mp->Fuu[1] = 0; mp->Fuu[2] = 0;
    mp->Fuv[0] = 0; mp->Fuv[1] = 0; mp->Fuv[2] = 0;
    mp->Fvv[0] = 0; mp->Fvv[1] = 0; mp->Fvv[2] = 0;
    mp->F3[0] = mp->F3raw[0] = 0.0;
    mp->F3[1] = mp->F3raw[1] = 0.0;
    mp->F3[2] = mp->F3raw[2] = 1.0;
    }
  double ComputeLaplacian(double u, double v)
    { return -1.0; }
};

class PDESplineWrapper : virtual public IMedialPDEProblem
{
public:
  /** A wrapper around a medial spline */
  PDESplineWrapper(DynamicBSpline2D *spline, unsigned int level)
  { 
    this->spline = spline;
    this->level = level;
    sgd = SplineGridDefinition::newSGDForCache(spline, level); 
  }
    
  /* Function to get the 2-jet at u, v pair */
  void ComputeJet2(MedialPoint *mp)
  {
    // Find the u, v index for the input coordinates
    unsigned int iu = sgd->getIndexForUV(0, mp->u);
    unsigned int iv = sgd->getIndexForUV(1, mp->v);

    cout << mp->u << " : " << iu << " ;" << mp->v << " : " << iv << endl;

    // Interpolate the grid point for this medial point
    spline->interpolateGridPoint(*sgd, iu, iv, 0, 0, 0, 2, mp->F.data_block());
    spline->interpolateGridPoint(*sgd, iu, iv, 1, 0, 0, 2, mp->Fu.data_block());
    spline->interpolateGridPoint(*sgd, iu, iv, 0, 1, 0, 2, mp->Fv.data_block());
    spline->interpolateGridPoint(*sgd, iu, iv, 2, 0, 0, 2, mp->Fuu.data_block());
    spline->interpolateGridPoint(*sgd, iu, iv, 1, 1, 0, 2, mp->Fuv.data_block());
    spline->interpolateGridPoint(*sgd, iu, iv, 0, 2, 0, 2, mp->Fvv.data_block());

    // Compute the normal
    SMLVec3f Nraw = cross_3d(mp->Xu(), mp->Xv());
    SMLVec3f N = Nraw / Nraw.magnitude();
    memcpy(mp->F3raw.data_block(), Nraw.data_block(), 3 * sizeof(float));
    memcpy(mp->F3.data_block(), N.data_block(), 3 * sizeof(float));
  }
  
  double ComputeLaplacian(double u, double v)
    { return -0.7; }  
    
private:
  // Reference to the medial spline
  DynamicBSpline2D *spline;
  SplineGridDefinition *sgd;
  unsigned int level;
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

void glDrawQuadStripElements(unsigned short width,unsigned short height) 
{
  // Allocate the index array
  int size = width*2;
  unsigned short *index = new unsigned short[size];

  unsigned short iStart = 0,iStart2 = width;

  for (int j=0;j<height-1;j++)
    {
    int iIndex = 0;
    for (int i=0;i<width;i++)
      {
      index[iIndex++] = iStart++;
      index[iIndex++] = iStart2++;
      }
    glDrawElements(GL_QUAD_STRIP,size,GL_UNSIGNED_SHORT,index);
    }

  delete index;
}

void glDrawWireframeElements(unsigned short width, unsigned short height)
{
  // Horizontal elements as lines
  unsigned int size = height + width;
  unsigned short *index = new unsigned short[size];

  // Draw the vertical lines
  unsigned short iStart = 0;
  for(int i = 0; i < width; i++)
    {
    int iIndex = 0;
    for(int j = 0; j < height; j++)
      { index[iIndex++] = iStart++; }
    glDrawElements(GL_LINE_STRIP, height, GL_UNSIGNED_SHORT, index);
    }

  // Draw the horizontal lines
  for(int j = 0; j < height; j++)
    {
    iStart = j;
    int iIndex = 0;
    for(int i = 0; i < width; i++)
      { index[iIndex++] = iStart; iStart+=height; }
    glDrawElements(GL_LINE_STRIP, width, GL_UNSIGNED_SHORT, index);
    }
}

class PDESplineRenderer : public GLDisplayListRenderer
{
public:
  PDESplineRenderer(MedialPDESolver *solver)
    {
    this->solver = solver;
    matMedial = new GLMaterial(GL_FRONT_AND_BACK, 
      GLColor(0.1), GLColor(0.4, 0.4, 0.8), GLColor(0.15), 64); 
      
    matBoundary = new GLMaterial(GL_FRONT_AND_BACK, 
      GLColor(0.1), GLColor(0.4));
    }

  virtual void build()
    {
    // Set the display attributes
    glPushAttrib(GL_LIGHTING_BIT);
    glEnable(GL_LIGHTING);

    // Center the object
    glPushMatrix();
    glTranslated(-0.5, -0.5, -0.0);

    // Enable vector arrays
    glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT); 
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    // Get the number of points
    unsigned int m = solver->GetNumberOfUPoints();
    unsigned int n = solver->GetNumberOfVPoints();

    // Initialize the medial colors
    matMedial->apply();

    // Supply the vertex list (medial surface)
    MedialPoint *mp = solver->GetAtom(0, 0);
    glVertexPointer(3, GL_FLOAT, sizeof(MedialPoint), mp->F.data_block());
    glNormalPointer(GL_FLOAT, sizeof(MedialPoint), mp->F3.data_block());

    // Build the quad array
    // glDrawQuadStripElements(m, n);
    glDrawQuadStripElements(m, n);
    
    // Display the boundary
    matBoundary->apply();

    // Supply the first boundary array
    glVertexPointer(3, GL_FLOAT, sizeof(MedialPoint), mp->bp[0].X.data_block());
    glNormalPointer(GL_FLOAT, sizeof(MedialPoint), mp->bp[0].N.data_block());

    // Draw the wireframe
    glDrawWireframeElements(m, n);

    // Supply the second boundary array
    glVertexPointer(3, GL_FLOAT, sizeof(MedialPoint), mp->bp[1].X.data_block());
    glNormalPointer(GL_FLOAT, sizeof(MedialPoint), mp->bp[1].N.data_block());

    // Draw the wireframe
    glDrawWireframeElements(m, n);

    // Restore the client state
    glPopClientAttrib();

    // Restore the GL state
    glPopMatrix();
    glPopAttrib();
    }

  

private:
  MedialPDESolver *solver;
  GLMaterial *matMedial, *matBoundary;
};

void testFourierFit(FourierSurface *xFourier)
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
  xFourier->FitData(0, nPoints, uPoints, 1, vPoints, 1, xPoints, 1, 4, 4);
  xFourier->FitData(1, nPoints, uPoints, 1, vPoints, 1, yPoints, 1, 4, 4);
  xFourier->FitData(2, nPoints, uPoints, 1, vPoints, 1, zPoints, 1, 4, 4);
}

struct DiscreteAtom 
{ 
  double u, v, x, y, z; 
  unsigned int iu, iv; 
};

void FitDiscreteMRep(const char *file, FourierSurface *surface)
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
  surface->FitData(0, xAtoms.size(), uu, 1, vv, 1, xx, 1, 5, 5);
  surface->FitData(1, xAtoms.size(), uu, 1, vv, 1, yy, 1, 5, 5);
  surface->FitData(2, xAtoms.size(), uu, 1, vv, 1, zz, 1, 5, 5);

  // Clean up
  delete xx; delete yy; delete zz; delete uu; delete vv;
}

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
  void SetEdgeImage(GrayImage *xImage);

  // Compute the match
  float computeBoundaryMeasure(const MedialPoint &mp, int d);
  float computeCrestBoundaryMeasure(const MedialPoint &mp)
    { return computeBoundaryMeasure(mp, 0); }
  
private:
  GrayImage *xImage;
};

#include "itkImageFileReader.h"
#include "itkLinearInterpolateImageFunction.h"

void TestImageIO(const char *fname)
{
  cout << "Reading image " << fname << endl;

  // Load an image from ITK
  typedef itk::Image<short, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(fname);
  fltReader->Update();
  ImageType::Pointer imgInput = fltReader->GetOutput();

  cout << "Creating image cube" << fname << endl;

  // Create an image cube from this image
  ImageCubeITK<short> cubeInput;
  cubeInput.SetImage(imgInput, 0);

  // Get the size object
  ImageType::SizeType szImage = imgInput->GetBufferedRegion().GetSize();

  cout << "Generating random indices" << fname << endl;

  //  Generate an array of random indices into the image
  unsigned int nTests = 1000000, iTest, q;
  vector<SMLVec3f> xIndex(nTests);
  for(iTest = 0; iTest < nTests; iTest++)
    {
    xIndex[iTest][0] = rand() * (szImage[0] * 1.0f) / RAND_MAX;
    xIndex[iTest][1] = rand() * (szImage[1] * 1.0f) / RAND_MAX;
    xIndex[iTest][2] = rand() * (szImage[2] * 1.0f) / RAND_MAX;
    }

  
  // Interpolate the image at these locations using the cube
  clock_t tCube, tITK;
  float xTestCube = 0, xTestITK = 0;

  cout << "Running ITK interpolation" << fname << endl;
  
  // Interpolate using the ITK code
  typedef itk::LinearInterpolateImageFunction<ImageType, float> FuncType;
  FuncType::Pointer funInterp = FuncType::New();
  funInterp->SetInputImage(imgInput);

  tITK = clock();
  for(q = 0; q < 100; q++)
    for(iTest = 0; iTest < nTests; iTest++)
      {
      itk::ContinuousIndex<float, 3> idx(xIndex[iTest].data_block());
      xTestITK += funInterp->EvaluateAtContinuousIndex(idx);
      }
  tITK = (clock() - tITK) * 1.0 / CLOCKS_PER_SEC;

  cout << "Running image cube interpolation" << fname << endl;
  tCube = clock();
  for(q = 0; q < 100; q++)
    for(iTest = 0; iTest < nTests; iTest++)
      {
      float *p = xIndex[iTest].data_block();
      xTestCube += cubeInput.interpolateVoxel(p[0], p[1], p[2]);
      }
  tCube = (clock() - tCube) * 1.0 / CLOCKS_PER_SEC;
  
  cout << "ITK returned " << xTestITK << " in " << (1000 * tITK) << " ms." << endl;
  cout << "CUB returned " << xTestCube << " in " << (1000 * tCube) << " ms." << endl;
}

int main(int argc, char *argv[])
{
  // Command line parameter variables
  unsigned int xResolution = 8;
  char *fnImage = NULL, *fnModel = NULL;

  // Read command line parameters 
  for(unsigned int iArg = 1; iArg < argc; iArg++)
    {
    if(strcmp(argv[iArg],"-r") == 0 && iArg < argc - 1)
      xResolution = atoi(argv[++iArg]);
    if(strcmp(argv[iArg],"-i") == 0 && iArg < argc - 1)
      fnImage = argv[++iArg];
    if(strcmp(argv[iArg],"-m") == 0 && iArg < argc - 1)
      fnModel = argv[++iArg];
    }

  // Test the imaging
  if(fnImage)
    TestImageIO(fnImage);
  
  /* Step 1. Test the Differential Geometry code */
  double u = 0.3, v = 0.4;
  
  // Build and display a spline
  unsigned int xPointsU = 6, xPointsV = 6;
  DynamicBSpline2D *spline = new DynamicBSpline2D(xPointsU, xPointsV, 3, 0);
  for(unsigned int i = 0; i <= xPointsU; i++)
    for(unsigned int j = 0; j <= xPointsV; j++)
      spline->setControl(i, j, 2, 0.0);
    
  PDESplineWrapper psw(spline, 3);
    
  // Compute the Jet
  MyProblem1 p1;
  MedialPoint mp; mp.u = u; mp.v = v;
  p1.ComputeJet2(&mp);

  // Another problem and jet
  MyProblem2 p2;
  p2.ComputeJet2(&mp);

  // Compute the differential geometry
  GeometryDescriptor<float> gd;
  gd.SetJet(&mp);
  gd.PrintSelf(cout);

  // Try a fourier problem
  FourierSurface xFourier(5,5,3);
  srand(clock());
  for(unsigned int k = 0; k < 4; k++)
    for(unsigned int l = 0; l < 4; l++)
      {
      double v = (rand() * 1.0 / RAND_MAX) - 0.5;
      v *= pow(0.5, (double)(k+l));
      // v *= exp((double)(-k - l)) / exp(2.0);
      xFourier.SetCoefficient(k, l, 2, v);
      cout << "u: " << k << ";  v: " << l << "  v: " << v << endl;
      }

  // Test data fitting
  if(fnModel)
    FitDiscreteMRep(fnModel, &xFourier);
  else
    testFourierFit(&xFourier);

  // Create a wrapper around the Fourier problem
  PDEFourierWrapper pfw(&xFourier);
  
  /* Step 2. Evaluate the equation at each site */
  // Initialize the solver
  unsigned int m = xResolution * xPointsU;
  unsigned int n = xResolution * xPointsV;
  MedialPDESolver mps(m+1, n+1);

  time_t t = clock();
  mps.Solve(&pfw);
  cout << "*** Elapsed Time: " << ((clock() - t) * 1.0 / CLOCKS_PER_SEC) << " ms." << " ***" << endl;
  
  // Save the atoms
  /** 
  for(unsigned int i = 0; i < m; i++)
    for(unsigned int j = 0; j < n; j++)
      {
      MedialPoint *mp = mps.GetAtom(i, j);
      cout << "Medial Point [" << i << "][" << j << "] : ";
      cout << " R = " << mp->R() << "; ";
      cout << " B = { " << mp->bp[0].X << ", " << mp->bp[1].X << "} " << endl;
      }
      */

  // Create a renderer
  PDESplineRenderer *rndPDESpline = new PDESplineRenderer(&mps);

  // Initialize the GL environment
  GLDisplayDriver::init(argc, argv);

  // Add the important renderers
  GLDisplayDriver::addRenderer(new DefaultLightRenderer());
  GLDisplayDriver::addRenderer(new StarfieldRenderer());

  // Add our spline renderer
  GLDisplayDriver::addRenderer(rndPDESpline);

  // Start GLUT
  glutMainLoop();
}




