//#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
//#include <crtdbg.h>

// My includes
#include "glui.h"
#include <sstream>
#include <fstream>
#include "imaging.h"
#include "imatch.h"
#include "optimization.h"
#include <registry.h>
#include "colors.h"
#include "splgeo.h"
#include <cstdarg>
#include "smlmath.h"
#include <vnl/vnl_rotation_matrix.h>
#include <vnl/vnl_matrix_fixed.txx>

#ifdef _MSC_VER
  #include "windows.h"
#endif

#ifndef M_PI
const double M_PI = acos(-1.0);
#endif

// Chdir support
#ifdef WIN32
  #include <direct.h>
  #define CHDIR(a) _chdir(a)
  #define GETCWD(a,b) _getcwd(a,b)
#else
  #include <unistd.h>
  #define CHDIR(a) chdir(a)
  #define GETCWD(a,b) getcwd(a,b)
#endif

#define RESLN 5

// Instantiate VNL matrices
template class vnl_matrix_fixed<float,2,2>;
template class vnl_matrix_fixed<float,4,4>;


/******************************************************************************
    Global Variables
 ******************************************************************************/
// The spline that can be edited
MSpline *spline;
SplineDataCache *splineData;

// Spline renderer
BSplineRenderer *rndSpline;

// Mode manager
ModeManager *modeManager = new ModeManager();

// Some fonts
GLFont *fntCourier18 = new GLFont("Unispace",18,false,false);
GLFont *fntCourier12 = new GLFont("Courier New",12,false,false);

// Command Shell
CommandShell *commandShell = NULL;

// Status window;
GLFadingWindow rndFadeStatus(fntCourier18,GLColor(0.8,0.8,0.5));

// A distance transform image
ImagingSystem imaging;

// The current image object
AbstractImage3D *image3D;

// The type of image that is currently in use

// Settings
Registry settingsAll,*settings,*settingsUI;

// Global Event Handler
GlobalEventHandler globalEventHandler;

// Matlab handler
// Matlab *matlab = NULL;

// Script processor
ScriptProcessor scriptProcessor;

// Geodesic thingy
GeodesicRenderer geodesicRenderer;

// Whether we are in interactive mode or not
bool flagInteractive = true;

/******************************************************************************
    Small Functions
 ******************************************************************************/
void glNormal(const SMLVec3f &v) {
  glNormal3f(v[0],v[1],v[2]);
}

void glVertex(const SMLVec3f &v) {
  glVertex3f(v[0],v[1],v[2]);
}

void processCommand(const char *str,ostream &sout);

/**
 * This is a utility method for calling glDrawElements for strides of data
 */
void glDrawElements(GLenum mode,unsigned short startIdx,short count,short* pattern,int patternSize) {
  int size = count*patternSize;
  unsigned short *index = new unsigned short[size];

  int k =0;
  for (int i=0;i<size;i++)
    {
    index[i] = (i==0) ? startIdx : index[i-1] + pattern[(i-1) % patternSize];
    }

  glDrawElements(mode,size,GL_UNSIGNED_SHORT,index);

  delete index;
}

/**
 * This is a utility method to draw a quad patch using glDrawElements
 */
void glDrawQuadStripElements(unsigned short width,unsigned short height) {
  // Allocate the index array
  int size = width*2;
  unsigned short *index = new unsigned short[size];

  unsigned short iStart = 0,iStart2 = width;

  for (int j=0;j<height;j++)
    {
    int iIndex = 0;
    for (int i=0;j<width;i++)
      {
      index[iIndex++] = iStart++;
      index[iIndex++] = iStart2++;
      }
    glDrawElements(GL_QUAD_STRIP,size,GL_UNSIGNED_SHORT,index);
    }

  delete index;
}

/*********************************************************
 Light Definition
  *******************************************************/
extern GLfloat EyeAz,EyeDist,EyeEl;

void LightRenderer::onDraw() {

  // Transform the view into eye coordinates
  glPushMatrix();

  glRotatef(-EyeAz, 0, 1, 0);
  glRotatef(-EyeEl, 1, 0, 0);
  glTranslatef(0, 0, EyeDist);

  // Run GL initialization code - establish a scene
  GLfloat light_position[] = { 0.0f, 0.0f, 1.0f, 0.0f};
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);

  glPopMatrix();
}

SplineIndexArraySet::SplineIndexArraySet(int level,int maxLevel,int m,int n) {
  // Resize into the number of patches
  idxPatchQuad.resize(m,n);

  // Patch size
  int dArray = (1 << level) + 1;
  int dPatch = (1 << maxLevel) + 1;
  int step = 1 << (maxLevel-level);

  for (int jPatch=0;jPatch<n;jPatch++)
    {
    for (int iPatch=0;iPatch<m;iPatch++)
      {
      // Resize the array
      idxPatchQuad(iPatch,jPatch).resize(dArray*2,dArray-1);

      int idxLow = 0;
      int idxHigh = dPatch * step;

      // Fill out the array values
      for (int j=0;j<dArray-1;j++)
        {
        for (int i=0;i<dArray*2;i+=2)
          {
          idxPatchQuad(iPatch,jPatch)(i,j) = idxLow;
          idxPatchQuad(iPatch,jPatch)(i+1,j) = idxHigh;

          idxHigh += step;
          idxLow += step;
          }

        idxHigh += dPatch*(step-1);
        idxLow += dPatch*(step-1);
        }
      }
    }
}

SplineIndexArraySet::~SplineIndexArraySet() {
}

/******************************************************************************
    BSplineRenderer
 ******************************************************************************/

BSplineRenderer::BSplineRenderer(DynamicBSpline2D *spline,SplineDataCache *dataCache, int level)
{
  // Set the spline references
  this->spline = spline;
  this->dataCache = dataCache;
  this->splineGrid = dataCache->getGrid();

  // Set the max and min levels
  maxLevel = dataCache->getMaxLevel();
  minLevel = level;

  // Save the steps
  dPatchMax = (1 << maxLevel) + 1;
  dPatchMin = (1 << minLevel) + 1;

  // Set patch sizes
  pw = spline->dim(0)-2;
  ph = spline->dim(1)-2;

  // Get grid sizes
  gw = splineGrid->getSize(0);
  gh = splineGrid->getSize(1);

  // Create the quad array sets for each level
  arraySet = new SplineIndexArraySet*[maxLevel+1];
  for (int l=0;l<=maxLevel;l++)
    {
    arraySet[l] = new SplineIndexArraySet(l,maxLevel,pw,ph);
    }

  // Initialize the 'general' renderers
  rndControl.parent = this;
  rndControlIndex.parent = this;
  rndWire.parent = this;

  // Initialize the 'patch' renderers
  rndPatchFlat.resize(pw,ph);
  rndPatchSeeThrough.resize(pw,ph);   
  rndPatchUVIndex.resize(pw,ph);  
  rndPatchTrimCurve.resize(pw,ph);
  rndPatchBoundary[0].resize(pw,ph);  
  rndPatchBoundary[1].resize(pw,ph);  

  UV.resize(pw,ph);
  cmBoundary[0].resize(pw,ph);
  cmBoundary[1].resize(pw,ph);

  for (int i=0;i<pw;i++)
    {
    for (int j=0;j<ph;j++)
      {
      rndPatchFlat(i,j).init(this,i,j,SplineDataCache::MEDIAL);
      rndPatchSeeThrough(i,j).init(this,i,j,SplineDataCache::MEDIAL);
      rndPatchUVIndex(i,j).init(this,i,j,SplineDataCache::MEDIAL);
      rndPatchTrimCurve(i,j).init(this,i,j,SplineDataCache::BOUNDARY);
      rndPatchBoundary[0](i,j).init(this,i,j,SplineDataCache::BOUNDARY);
      rndPatchBoundary[0](i,j).side = 0;
      rndPatchBoundary[1](i,j).init(this,i,j,SplineDataCache::BOUNDARY);
      rndPatchBoundary[1](i,j).side = 1;


      // Initialize the UV index array            
      UV(i,j).resize(dPatchMax,dPatchMax);

      float scale = 1.0f / (dPatchMax-1);

      for (int iGrid=0;iGrid<dPatchMax;iGrid++)
        {
        for (int jGrid=0;jGrid<dPatchMax;jGrid++)
          {
          UV(i,j)(iGrid,jGrid)[0] = iGrid * scale;
          UV(i,j)(iGrid,jGrid)[1] = jGrid * scale;
          UV(i,j)(iGrid,jGrid)[2] = 1.0;
          }
        }


      }
    }

  // Create the materials
  makeMaterials();

  // Set the surface mode
  surfMode = SURFACE;
  colorMapMode = CMAP_NONE;
}

BSplineRenderer::~BSplineRenderer() {
  for (int l=0;l<=maxLevel;l++)
    {
    delete arraySet[l];
    }
  delete[] arraySet;
}

void BSplineRenderer::makeMaterials() {
  // Initialize the materials
  for (int i=0;i<3;i++)
    {
    matFlatSurface[i].ambient = GLColor(0.1,0.1,0.1);
    matFlatSurface[i].diffuse = GLColor(0.4,0.6+0.1*i,0.2);
    matFlatSurface[i].specular = GLColor(0.15,0.15,0.15);
    matFlatSurface[i].shininess = 64;

    matSeeThruSurface[i].ambient = GLColor(0.1,0.1,0.1,0.7);
    matSeeThruSurface[i].diffuse = GLColor(0.4+0.1*i,0.4+0.1*i,0.4+0.1*i,0.7);
    matSeeThruSurface[i].specular = GLColor(0.2,0.2,0.2,0.15);
    matSeeThruSurface[i].shininess = 32;

    matBoundary[i].ambient = GLColor(0.1,0.1,0.1,0.3);
    matBoundary[i].diffuse = GLColor(0.6,0.6+0.1*i,0.2,0.3);
    matBoundary[i].specular = GLColor(0.6,0.6,0.6,0.15);
    matBoundary[i].shininess = 128;
    }
}

void HSVtoRGB( float *r, float *g, float *b, float h, float s, float v )
{
  int i;
  float f, p, q, t;
  if ( s == 0 )
    {
    // achromatic (grey)
    *r = *g = *b = v;
    return;
    }
  h /= 60;            // sector 0 to 5
  i = (int) h;
  f = h - i;          // factorial part of h
  p = v * ( 1 - s );
  q = v * ( 1 - s * f );
  t = v * ( 1 - s * ( 1 - f ) );
  switch ( i )
    {
    case 0:
      *r = v;
      *g = t;
      *b = p;
      break;
    case 1:
      *r = q;
      *g = v;
      *b = p;
      break;
    case 2:
      *r = p;
      *g = v;
      *b = t;
      break;
    case 3:
      *r = p;
      *g = q;
      *b = v;
      break;
    case 4:
      *r = t;
      *g = p;
      *b = v;
      break;
    default:        // case 5:
      *r = v;
      *g = p;
      *b = q;
      break;
    }
}

void BSplineRenderer::makeBoundaryColorPatch(PatchDataCache *pdc, int d, int iPatch, int jPatch) {
  // If there is no map, then return
  if (colorMapMode == CMAP_NONE)
    return;

  // Get the number of pixels to paint
  int n = pdc->MP.width() * pdc->MP.height() + pdc->getCrestSize();

  // Refresh the curvatures of the data cache
  pdc->refreshCurvatures(pdc->getMaxLevel());

  // Resize the map to fit the data
  cmBoundary[d](iPatch,jPatch).resize(pdc->MP.width(),pdc->MP.height(),pdc->MP.auxSize());

  // Get the hue base for colorization
  int hueBase = settings->getIntValue("colormap.hueBase",0);
  int clrPatch = (iPatch%2) + (jPatch%2);

  // Make a Gaussian color patch
  if (colorMapMode == CMAP_GAUSS)
    {
    for (int i=0;i<n;i++)
      {
      SMLVec3f rgb;
/*          
            if(pdc->in(i) & 0x01) {

                

                float k1 = pdc->MP(i).bp[d].kappa1;
                float k2 = pdc->MP(i).bp[d].kappa2;
                float r = pdc->MP(i).R();

                if(k2 >= 0) {
                    rgb.x = 0.2;
                    rgb.y = 0.2;
                    rgb.z = 0.8;
                }
                else if(k2 + 1.0/r < -0.00001f) {
                    rgb.x = 0.8;
                    rgb.y = 0.8;
                    rgb.z = 0.8;
                }
                else if(k1 >= 0) {
                    rgb.x = 0.8;
                    rgb.y = 0.8;
                    rgb.z = 0.2;
                }
                else {
                    rgb.x = 0.8;
                    rgb.y = 0.2;
                    rgb.z = 0.2;
                }
            }

            else {
                rgb.x = 0.2;
                rgb.y = 0.8;
                rgb.z = 0.8;
            }

            cmBoundary[d](iPatch,jPatch)(i) = rgb;
*/          

      /*
      cmBoundary[d](iPatch,jPatch)(i).z = pdc->sinTheta2(i);
      cmBoundary[d](iPatch,jPatch)(i).y = pdc->sinTheta2(i);
      cmBoundary[d](iPatch,jPatch)(i).x = 1;
      */

      if (pdc->in(i) & 0x01)
        {
        BoundaryPoint &bp = pdc->MP(i).bp[d];
        float gc = bp.kappa1 * bp.kappa2;

        float hue = 60 + 120*hueBase + 60 * atan(gc * 0.1f);
        float val = 0.5f;
        float sat = 0.6f + 0.1f * clrPatch; // 0.8 + 0.1 * ((iPatch%2)+(jPatch%2));
        HSVtoRGB(&rgb[0],&rgb[1],&rgb[2],hue,sat,val);
        cmBoundary[d](iPatch,jPatch)(i) = rgb;
        }

      }
    }
  else if (colorMapMode == CMAP_KOEN)
    {

    for (int i=0;i<n;i++)
      {
      // Koenderink color mapping
      SMLVec3f rgb;

      if (pdc->in(i) & 0x01)
        {
        float k1 = pdc->MP(i).bp[d].kappa1;
        float k2 = pdc->MP(i).bp[d].kappa2;
        if (k1*k1+k2*k2 > 100000)
          {
          k1 = k1+0.1;
          }
        float C = log(k1*k1+k2*k2) / M_PI;
        float S = - 2 * atan((k1+k2)/(k1-k2)) / M_PI;

        // int idx = ((iPatch%2)+(jPatch%2));
        // float hue = 60 + 240 + 60 * S;
        float hue = 60 + 120*hueBase + 60 * S;
        float val = 0.5 + atan(C) / M_PI;
        float sat = 1.0f; // 0.8 + 0.1 * ((iPatch%2)+(jPatch%2));
        HSVtoRGB(&rgb[0],&rgb[1],&rgb[2],hue,sat,val);
        }
      cmBoundary[d](iPatch,jPatch)(i) = rgb;
      }
    }

  else if (colorMapMode == CMAP_FLETCH)
    {

    for (int i=0;i<n;i++)
      {
      // Koenderink color mapping
      SMLVec3f rgb;

      if (pdc->in(i) & 0x01)
        {
        MedialPoint &mp = pdc->MP(i);
        BoundaryPoint &bp = pdc->MP(i).bp[d];

        float am = dot_product(vnl_cross_3d(mp.Xu(),mp.Xv()),mp.N());
        float ab = dot_product(vnl_cross_3d(bp.Xu,bp.Xv),bp.N);
        float r = d ?  ab / am : -ab /am;

        // Map the ratio to a value between -1 and 1
        if (r >= 0)
          {
          float r1 = 2 * atan(log(fabs(r))/log(2.0)) / M_PI;                
          float hue = fmod((60 + 120*hueBase+60*r1),360);
          float sat = 1.0f;
          float val = 0.75f;
          HSVtoRGB(&rgb[0],&rgb[1],&rgb[2],hue,sat,val);
          }
        else
          {
          float hue = fmod(240.0 + 120.0*hueBase,360.0);
          float sat = 1.0f;
          float val = 1.0f;
          HSVtoRGB(&rgb[0],&rgb[1],&rgb[2],hue,sat,val);
          }
        }
      cmBoundary[d](iPatch,jPatch)(i) = rgb;
      }
    }

  else if (colorMapMode == CMAP_IMAGE)
    {

    // Get the image
    AbstractImage3D *img = imaging.getImage();
    if (img == NULL)
      return;

    // Different effects by image type
    if (imaging.getType() == ImagingSystem::DISTANCE)
      {

      SMLVec3f rgb;

      // Go through all the trimmed points
      int mode = settings->getIntValue("colormap.distanceToRGB.mode",0);
      
      if(mode == 0)
        {
        for (int i=0;i<pdc->idxTrimmed.size();i++)
          {

          // Get to the boundary point in question                
          int idx = pdc->idxTrimmed[i];
          BoundaryPoint &bp = pdc->MP(idx).bp[d];

          // Interpolate the image match at the point
          float dist = img->interpolateVoxel(bp.X[0],bp.X[1],bp.X[2]);               
          cmBoundary[d](iPatch,jPatch)(idx) = convertDistanceToRGB(dist);
          }
        }
      else
        {
        for (int i=0;i<pdc->idxTrimmed.size();i++)
          {

          // Get to the boundary point in question                
          int idx = pdc->idxTrimmed[i];
          BoundaryPoint &bp = pdc->MP(idx).bp[d];

          // Interpolate the image match at the point
          float dist = img->interpolateVoxel(bp.X[0],bp.X[1],bp.X[2]);               
          SMLVec3f &clr = cmBoundary[d](iPatch,jPatch)(idx);
          
          convertTScoreToRGB(dist, clr[0], clr[1], clr[2]);
          }
        }
      }
    else if (imaging.getType() == ImagingSystem::GRAYSCALE || imaging.getType() == ImagingSystem::BINARY)
      {

      // Get the scaling factor for the display
      int scale = settings->getIntValue("colormap.gradientMatchScale",0xffff);
      float scaleFactorVal = 1.0f / scale;

      SMLVec3f rgb;

      // Go through all the trimmed points
      for (int i=0;i<pdc->idxTrimmed.size();i++)
        {

        // Get to the boundary point in question                
        int idx = pdc->idxTrimmed[i];
        BoundaryPoint &bp = pdc->MP(idx).bp[d];

        // Interpolate the image match at the point
        SMLVec3f G;
        img->interpolateVoxelGradient(bp.X[0],bp.X[1],bp.X[2],G.data_block());
        float g = dot_product(G,bp.N) * scaleFactorVal;

        float hue = fmod(120 * (hueBase + g), 360);
        float sat = 1.0f;
        float val = 0.5 + 0.5 * g;

        HSVtoRGB(&rgb[0],&rgb[1],&rgb[2],hue,sat,val);
        cmBoundary[d](iPatch,jPatch)(idx) = rgb;
        }
      }
    }
}

void BSplineRenderer::onDraw() {
  glPushMatrix();
  glTranslated(-0.5,-0.5,-0.5);

  // Make sure that all data in the cache is up to date (because we render all the data on each pass)
  dataCache->refreshMedialAllPatches(minLevel);

  if (GLDisplayDriver::displayMode == GLDisplayDriver::NORM)
    {
    // Draw the solid surface
    if (surfMode & SURFACE)
      {
      for (int i=0;i<pw;i++)
        {
        for (int j=0;j<ph;j++)
          {
          rndPatchFlat(i,j).onDraw();                 
          }
        }
      }

    // Draw the wire frame
    if (surfMode & WIRE)
      {
      rndWire.onDraw();
      }

    // Draw control points
    if (surfMode & CONTROL)
      {
      rndControl.onDraw();
      }
    // Draw the boundary surface
    if (surfMode & IBSURFACE_LEFT)
      {

      // Make sure that all data in the cache is up to date (because we render all the data on each pass)
      dataCache->refreshBoundaryAllPatches(minLevel);

      for (int i=0;i<pw;i++)
        {
        for (int j=0;j<ph;j++)
          {
          rndPatchBoundary[0](i,j).onDraw();
          }
        }
      }

    // Draw the boundary surface
    if (surfMode & IBSURFACE_RIGHT)
      {

      // Make sure that all data in the cache is up to date (because we render all the data on each pass)
      dataCache->refreshBoundaryAllPatches(minLevel);

      for (int i=0;i<pw;i++)
        {
        for (int j=0;j<ph;j++)
          {
          rndPatchBoundary[1](i,j).onDraw();
          }
        }
      }

    // if(!((surfMode & IBSURFACE_LEFT) && (surfMode & IBSURFACE_RIGHT))) {
    if (surfMode)
      {
      // Make sure that all data in the cache is up to date (because we render all the data on each pass)
      dataCache->refreshBoundaryAllPatches(minLevel);

      for (int i=0;i<pw;i++)
        {
        for (int j=0;j<ph;j++)
          {
          rndPatchTrimCurve(i,j).onDraw();
          }
        }
      }
    }

  else if (GLDisplayDriver::displayMode == GLDisplayDriver::TRAN)
    {
    // Draw the see through surface
    if (surfMode & SEETHRU)
      {
      for (int i=0;i<pw;i++)
        {
        for (int j=0;j<ph;j++)
          {
          rndPatchSeeThrough(i,j).onDraw();
          // rndPatchUVIndex(i,j).onDraw();
          }
        }
      }

    }

  // cout << "d";


  glPopMatrix();
}

void BSplineRenderer::setSurfaceMode(int mode) {
  surfMode = mode;
}

void BSplineRenderer::setColorMapMode(ColorModes mode) {
  colorMapMode = mode;

  switch (mode)
    {
    case CMAP_NONE: 
      echoMessage("Color coding off"); 
      break;
    case CMAP_GAUSS: 
      echoMessage("Gaussian curvature color coding"); 
      break;
    case CMAP_KOEN:
      echoMessage("Koenderink's curvature color coding"); 
      break;
    case CMAP_FLETCH:
      echoMessage("Tom Fletcher's legality color coding"); 
      break;
    case CMAP_IMAGE:
      echoMessage("Image value color coding"); 
      break;
    };

  for (int i=0;i<rndPatchBoundary->height();i++)
    {
    for (int j=0;j<rndPatchBoundary->width();j++)
      {
      rndPatchBoundary[0](j,i).invalidate();
      rndPatchBoundary[1](j,i).invalidate();
      }
    }
}

int round(float f) {
  float ff = floor(f);
  return(f - ff < 0.5) ? (int) ff : (int)(ff+1);
}

void BSplineRenderer::findSampleUV(int x,int y,int &u,int &v)
{
  // Static variables to remember the last patch that the mouse was on
  static int iPatch = -1, jPatch = -1;
  float uvb[] = {0,0,0};
  int hit = 0;

  // Push the attributes
  glPushAttrib(GL_COLOR_BUFFER_BIT | GL_PIXEL_MODE_BIT);

  glDrawBuffer(GL_BACK);
  glReadBuffer(GL_BACK);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // GLDisplayDriver::worldTransform();
  glPushMatrix();
  glTranslated(-0.5,-0.5,-0.5);

  // Try the last patch we were on
  if (iPatch >= 0)
    {
    rndPatchUVIndex(iPatch,jPatch).onDraw();
    glReadPixels(x,y,1,1,GL_RGB,GL_FLOAT,uvb);
    hit++;
    }
  for (int i=0;i<pw;i++)
    {
    for (int j=0;j<ph;j++)
      {
      if (uvb[2]==1)
        break;

      if (i!=iPatch || j!=jPatch)
        {
        rndPatchUVIndex(i,j).onDraw();          
        glReadPixels(x,y,1,1,GL_RGB,GL_FLOAT,uvb);
        iPatch = i;
        jPatch = j;
        hit++;
        }
      }
    }

  glPopMatrix();

  // Check if it is blue
  if (uvb[2]==1)
    {
    u = round(uvb[0] * (splineGrid->patchStart(0,iPatch+1) - splineGrid->patchStart(0,iPatch)) + splineGrid->patchStart(0,iPatch));
    v = round(uvb[1] * (splineGrid->patchStart(1,jPatch+1) - splineGrid->patchStart(1,jPatch)) + splineGrid->patchStart(1,jPatch));
    }
  else
    {
    iPatch = jPatch = u = v = -1;
    }

  // cout << "f";

  glPopAttrib();
}


void BSplineRenderer::findControlPoint(int x,int y,int &ci,int &cj)
{
  glDrawBuffer(GL_BACK);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  GLDisplayDriver::worldTransform();

  glPushMatrix();
  glTranslated(-0.5,-0.5,-0.5);
  rndControlIndex.onDraw();
  glPopMatrix();

  glReadBuffer(GL_BACK);

  // Read the pixel
  GLubyte uvb[3];
  glReadPixels(x,y,1,1,GL_RGB,GL_UNSIGNED_BYTE,uvb);

  // Check if it is blue
  if (uvb[2]==128)
    {
    ci = uvb[0];
    cj = uvb[1];
    }
  else
    {
    ci = cj = -1;
    }

  glPopAttrib();

}



void BSplineRenderer::buildControl(int dl) {
  int i;

  glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);
  glEnableClientState(GL_VERTEX_ARRAY);

  // Get the control point array
  spline->getControlArray(C,0);

  // Specify the vertex array
  glVertexPointer(3,GL_FLOAT,0,C.getData());

  glNewList(dl,GL_COMPILE);

  glPushAttrib(GL_LINE_BIT | GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT);

  glDisable(GL_LIGHTING);
  //glEnable(GL_LINE_STIPPLE);
  //glLineStipple(3,0xC3C3);

  glColor3dv(clrControlPoly);

  // Draw lines connecting the control points
  for (i=0;i<C.width();i++)
    {
    short pattern = C.width();
    glDrawElements(GL_LINE_STRIP,i,C.height(),&pattern,1);
    }
  for (i=0;i<C.height();i++)
    {
    short pattern = 1;
    glDrawElements(GL_LINE_STRIP,i*C.width(),C.width(),&pattern,1);
    }

  /*
  // glBegin(GL_POINTS);
  for(int i=0;i<=B->dim(0);i++) {
      for(int j=0;j<=B->dim(1);j++) {
          // glVertex3f(B->getControl(i,j,0),B->getControl(i,j,1),B->getControl(i,j,2));
          glPushMatrix();
          glTranslated(B->getControl(i,j,0),B->getControl(i,j,1),B->getControl(i,j,2));
          glutSolidSphere(0.01,10,4);
          glPopMatrix();
      }
  }
  // glEnd();
  */

  glPopAttrib();

  glEndList();

  glPopClientAttrib();
}


void BSplineRenderer::buildControlIndex(int dl) {

  glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);

  // Get the control point array
  spline->getControlArray(C,0);

  // Specify the vertex array
  glVertexPointer(3,GL_FLOAT,0,C.getData());

  // Create a color index array
  GLubyte *clr = new GLubyte[3*C.width()*C.height()];
  int idx = 0;
  for (int j=0;j<C.height();j++)
    {
    for (int i=0;i<C.width();i++)
      {
      clr[idx++] = i;
      clr[idx++] = j;
      clr[idx++] = 128;
      }
    }

  glColorPointer(3,GL_UNSIGNED_BYTE,0,clr);

  glNewList(dl,GL_COMPILE);

  glPushAttrib(GL_POINT_BIT | GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT);

  glDisable(GL_LIGHTING);
  glPointSize(20);

  // Draw lines connecting the control points
  glDrawArrays(GL_POINTS,0,C.height()*C.width());

  glPopAttrib();

  glEndList();

  glPopClientAttrib();

  delete clr;
}

// This method draws the surface in index mode.  The color at each pixel is determined by
// u, v values at that pixel
void BSplineRenderer::buildUVIndexPatch(int dl,int iPatch,int jPatch) {
  // Refresh the medial data
  dataCache->refreshMedialPatch(iPatch,jPatch,minLevel);
  PatchDataCache *pdc = dataCache->patch(iPatch,jPatch);

  glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT); 
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);

  // Specify the vertex array
  glVertexPointer(3,GL_FLOAT,sizeof(MedialPoint),
                  pdc->MP.getData()->F.data_block());

  // Specify the color array - this array should not change with the life of the spline
  glColorPointer(3,GL_FLOAT,0,UV(iPatch,jPatch).getData());

  glNewList(dl,GL_COMPILE);

  glPushAttrib(GL_LIGHTING_BIT);
  glDisable(GL_LIGHTING);
  glShadeModel(GL_FLAT);

  // Now there is an array of vertices and an array of normals.  We can render it using quad strips
  for (int row=0;row<(1 << minLevel);row++)
    {
    unsigned int *idxArray = arraySet[minLevel]->getPatchQuadStripArray(iPatch,jPatch,row);
    glDrawElements(GL_QUAD_STRIP,(2 << minLevel)+2,GL_UNSIGNED_INT,idxArray);
    }

  glPopAttrib();

  glEndList();

  glPopClientAttrib();
}

void BSplineRenderer::buildWire(int dl) {/*
    // Refresh the medial surface
    dataCache->refreshMedialAllPatches(minLevel);

    glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    // Specify the vertex arrays
    glVertexPointer(3,GL_FLOAT,0,dataCache->X.getData());
    
    // Specify the vertex arrays
    glNormalPointer(GL_FLOAT,0,dataCache->N.getData());

    glNewList(dl,GL_COMPILE);

    glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT | GL_COLOR_BUFFER_BIT);

    //glEnable(GL_LINE_SMOOTH);
    //glEnable(GL_BLEND);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //glLineWidth(0.5);

    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);

    glColor3d(0.2,0.8,0.3);

    short pattern = 1;
    for(int u=0;u<gh;u++) {
        if(u%16)
            continue;
        glDrawElements(GL_LINE_STRIP,u*gw,gw,&pattern,1);
    }

    pattern = gw;
    for(u=0;u<gw;u++) {
        if(u%16)
            continue;
        glDrawElements(GL_LINE_STRIP,u,gh,&pattern,1);
    }

    glPopAttrib();

    glEndList();

    glPopClientAttrib();*/
}

/**
 * Make a display list for a flat patch at a given resolution
 */
void BSplineRenderer::buildSurfacePatch(int dl,int iPatch,int jPatch,GLMaterial materials[]) {

  // Refresh the medial data
  dataCache->refreshBoundaryPatch(iPatch,jPatch,minLevel);
  PatchDataCache *pdc = dataCache->patch(iPatch,jPatch);

  // Specify array data
  glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);

  // Specify the vertex array
  glVertexPointer(3,GL_FLOAT,sizeof(MedialPoint),pdc->MP(0).F.data_block());

  // Specify the color array
  glNormalPointer(GL_FLOAT,sizeof(MedialPoint),pdc->MP(0).F3.data_block());

  glNewList(dl,GL_COMPILE);

  // Specify the material and lighting properties
  glPushAttrib(GL_LIGHTING_BIT | GL_POLYGON_BIT);
  glEnable(GL_LIGHTING);
  glDisable(GL_COLOR_MATERIAL);
  glFrontFace(GL_CW);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);

  // glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

  // Apply the appropriate material for this patch
  materials[(iPatch%2) + (jPatch%2)].apply();

  // Draw the quad strips
  for (int row=0;row<(1 << minLevel);row++)
    {
    unsigned int *ia = arraySet[minLevel]->getPatchQuadStripArray(iPatch,jPatch,row);

    int iStart = 0;
    int w = (2 << minLevel)+2;
    //glDrawElements(GL_QUAD_STRIP,w,GL_UNSIGNED_INT,ia);

    while (iStart < w-3)
      {
      while (iStart < w-3 && !((pdc->in(ia[iStart]) & 0x01) && (pdc->in(ia[iStart+1]) & 0x01)))
        iStart+=2;

      int iEnd = iStart+2;
      while (iEnd < w && ((pdc->in(ia[iEnd]) & 0x01) && (pdc->in(ia[iEnd+1]) & 0x01)))
        iEnd+=2;

      if (iEnd-iStart >= 4)
        {
        glDrawElements(GL_QUAD_STRIP,iEnd-iStart,GL_UNSIGNED_INT,ia+iStart);
        }

      // Update the start
      iStart = iEnd+2;
      }
    }

  if (pdc->idxTrimStripIdx > 0)
    {
    glDrawElements(GL_TRIANGLES,pdc->idxTrimStripIdx,GL_UNSIGNED_INT,pdc->idxTrimStrip);
    }

  glPopAttrib();

  glEndList();

  glPopClientAttrib();    
}

/**
 * Make a display list for a flat patch at a given resolution
 */
void BSplineRenderer::buildBoundaryPatch(int dl,int iPatch,int jPatch,int side,GLMaterial materials[]) {
  // Refresh the medial data
  dataCache->refreshBoundaryPatch(iPatch,jPatch,minLevel);
  PatchDataCache *pdc = dataCache->patch(iPatch,jPatch);

  // Specify array data
  glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);

  // Specify the vertex array
  glVertexPointer(3,GL_FLOAT,sizeof(MedialPoint),pdc->MP(0).bp[side].X.data_block());

  // Specify the color array
  glNormalPointer(GL_FLOAT,sizeof(MedialPoint),pdc->MP(0).bp[side].N.data_block());

  // Specify color maps
  if (colorMapMode)
    {
    glEnableClientState(GL_COLOR_ARRAY);
    makeBoundaryColorPatch(pdc,side,iPatch,jPatch);
    glColorPointer(3,GL_FLOAT,0,cmBoundary[side](iPatch,jPatch).getData());
    }

  glNewList(dl,GL_COMPILE);

  // Specify the material and lighting properties
  glPushAttrib(GL_LIGHTING_BIT | GL_POLYGON_BIT | GL_COLOR_BUFFER_BIT);
  glEnable(GL_LIGHTING);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);

  // glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  glDisable(GL_CULL_FACE);
  // Apply the appropriate material for this patch
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,1);

  // The winding of the polygons depends on the side
  glFrontFace(side ? GL_CW : GL_CCW);

  // glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

  glDisable(GL_BLEND);

  if (colorMapMode)
    {
    // materials[(iPatch%2) + (jPatch%2)].apply();
    glEnable(GL_COLOR_MATERIAL);
    // glDisable(GL_LIGHTING);
    }
  else
    {
    materials[(iPatch%2) + (jPatch%2)].apply();
    }

  // Draw the quad strips
  for (int row=0;row<(1 << minLevel);row++)
    {
    unsigned int *ia = arraySet[minLevel]->getPatchQuadStripArray(iPatch,jPatch,row);


    int iStart = 0;
    int w = (2 << minLevel)+2;



    while (iStart < w-3)
      {
      while (iStart < w-3 && !((pdc->in(ia[iStart]) & 0x01) && (pdc->in(ia[iStart+1]) & 0x01)))
        iStart+=2;

      int iEnd = iStart+2;
      while (iEnd < w && ((pdc->in(ia[iEnd]) & 0x01) && (pdc->in(ia[iEnd+1]) & 0x01)))
        iEnd+=2;

      if (iEnd-iStart >= 4)
        {
        glDrawElements(GL_QUAD_STRIP,iEnd-iStart,GL_UNSIGNED_INT,ia+iStart);
        }

      // Update the start
      iStart = iEnd+2;
      }
    }       

  // Draw the triangles of the trimming strip
  // glDisable(GL_LIGHTING);
  // glColor3d(1,0,0);    
  // glDrawElements(GL_POINTS,pdc->idxTrimStripIdx,GL_UNSIGNED_INT,pdc->idxTrimStrip);
  // glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  // glColor3d(1,1,0);
  // glDrawElements(GL_TRIANGLES,pdc->idxTrimStripIdx,GL_UNSIGNED_INT,pdc->idxTrimStrip);
  if (pdc->idxTrimStripIdx > 0)
    {
    glDrawElements(GL_TRIANGLES,pdc->idxTrimStripIdx,GL_UNSIGNED_INT,pdc->idxTrimStrip);
    }
/*
    if(colorMapMode) {
        glDisable(GL_LIGHTING);
        glColor3d(1,1,1);
        glBegin(GL_LINES);
        
        for(int p=0;p<pdc->MP.size();p++) {
            if(pdc->in(p) & 0x01) {
                BoundaryPoint &bp = pdc->MP(p).bp[side];
                glVertex(bp.X);
                glVertex(bp.X + bp.N * 0.2 * pdc->MP(p).R());
            }
        }
        glEnd();
    }
*/



  glPopAttrib();

  glEndList();

  glPopClientAttrib();    
}


/**
 * Make a display list for a flat patch at a given resolution
 */
void BSplineRenderer::buildPatchTrimCurve(int dl,int iPatch,int jPatch) {
  // Refresh the medial data
  dataCache->refreshBoundaryPatch(iPatch,jPatch,minLevel);
  PatchDataCache *pdc = dataCache->patch(iPatch,jPatch);

  // Specify array data
  glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);
  glEnableClientState(GL_VERTEX_ARRAY);
  // glEnableClientState(GL_NORMAL_ARRAY);

  // Specify the vertex array
  glVertexPointer(3,GL_FLOAT,sizeof(MedialPoint),pdc->MP(0).F.data_block());

  glNewList(dl,GL_COMPILE);

  // Specify the material and lighting properties
  glPushAttrib(GL_LIGHTING_BIT);

  glDisable(GL_LIGHTING);
  glColor3dv(clrTrimCurve);

  int t0 = 0;
  for (int cv=0;cv<pdc->idxTrimCurvesIdx;cv++)
    {
    glDrawArrays(GL_LINE_STRIP,pdc->MP.auxOffset(t0),pdc->idxTrimCurves[cv]-t0);
    t0 = pdc->idxTrimCurves[cv];
    }

  /*
  if(pdc->getCrestSize() > 0) {
      glDrawArrays(GL_LINE_STRIP,pdc->MP.auxOffset(0),pdc->getCrestSize());
  }
  */
  glPopAttrib();

  glEndList();

  glPopClientAttrib();    
}



class UVImageRnd : public GLDisplayListRenderer
  {
public:
  void build();
  void onDraw();

  void reset() {
    ts = -2;
  }
/*
    UVImage() {
        reset();
    }*/

private:
  long ts;
  };

void UVImageRnd::onDraw() {
  if (ts < spline->getControlTimeStamp())
    {
    invalidate();
    }
  GLDisplayListRenderer::onDraw();
}

void UVImageRnd::build() {
  int i,j;

  glPushMatrix();
  glScaled(320,320,1);

  glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT | GL_COLOR_BUFFER_BIT);
  glEnable(GL_BLEND);
  glEnable(GL_LINE_SMOOTH);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);    
  glLineWidth(3);

  glDisable(GL_LIGHTING);

  if(imaging.getType() == ImagingSystem::DISTANCE)
    {
    for (i=0;i<splineData->patch.width();i++)
      {
      for (j=0;j<splineData->patch.height();j++)
        {
        PatchDataCache *pdc = splineData->patch(i,j);

        glBegin(GL_POINTS);
        for(int k = 0; k < pdc->MP.size(); k++)
          {
          MedialPoint &mp = pdc->MP(k);
          if(mp.sinTheta2 > 0)
            {
            float tavg = 0.0f, t;
            for(t = -1; t < 0.0; t+=0.1)
              {
              SMLVec3f xt = mp.X() - (t * mp.R())* mp.bp[0].N;
              float dist = 
                imaging.getImage()->interpolateVoxel(xt[0],xt[1],xt[2]);
              tavg += dist;
              }
            for(t = 0; t <= 1.0; t+=0.1)
              {
              SMLVec3f xt = mp.X() + (t * mp.R()) * mp.bp[1].N;
              float dist = 
                imaging.getImage()->interpolateVoxel(xt[0],xt[1],xt[2]);
              tavg += dist;
              }

            tavg /= 20;
            
            float r, g, b;
            convertTScoreToRGB(tavg, r, g, b);
            glColor3d(r, g, b);
            glVertex2d(mp.u,mp.v);
            }
          }
        glEnd();
        }
      }
    }

  glColor3dv(clrTrimCurve);

  for (i=0;i<splineData->patch.width();i++)
    {
    for (j=0;j<splineData->patch.height();j++)
      {
      PatchDataCache *pdc = splineData->patch(i,j);

      int t0 = 0;
      for (int cv=0;cv<pdc->idxTrimCurvesIdx;cv++)
        {
        int t;
        glBegin(GL_LINE_STRIP);
        for (t=t0;t<pdc->idxTrimCurves[cv];t++)
          {
          glVertex2d(pdc->MP.aux(t).u,pdc->MP.aux(t).v);
          }
        t0 = t;
        glEnd();
        }
      }
    }
  

  ts = spline->getControlTimeStamp();

  glPopAttrib();

  glPopMatrix();  
}

UVImageRnd uvImage;





GLOutputWindow::GLOutputWindow(GLFont *font,GLColor color) {
  this->font = font;
  this->color = color;
  anchorTop = false;
  anchorX = anchorY = 5;
}

void GLOutputWindow::clear() {
  lines.clear();
}

void GLOutputWindow::addLine(const char *format,...) {
  va_list args;
  char buffer[400];

  va_start(args, format);
  vsprintf(buffer, format, args);
  va_end(args);

  char *t = strtok(buffer,"\n\r");
  while (t)
    {
    lines.push_back(string(t));
    t = strtok(NULL,"\n\r");
    }


}

void GLOutputWindow::setAnchor(bool top,int x,int y) {
  this->anchorTop = top;
  this->anchorX = x;
  this->anchorY = y;
}

void GLOutputWindow::onDraw() {
  // const int fsHeader = 18;
  const int fsBody = font->getHeight();
  const int fsHeader = fsBody + 4;
  const int fsSpace = font->getHeight() / 6;

  glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
  glDisable(GL_LIGHTING);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glColor4fv(color.fv);

  // First of all, if there are no lines, do not draw anything
  if (lines.size() > 0)
    {
    // Count the number of lines
    int ls = (int)lines.size();
    int h = fsHeader + fsBody*(ls-1) + fsSpace * (ls+1);

    int ancY = (anchorY < 0 || (anchorTop && anchorY==0)) ? GLDisplayDriver::height + anchorY : anchorY; 
    int y = anchorTop ? ancY : ancY + h;
    int x = anchorX;

    list<string>::iterator it;
    for (it=lines.begin();it!=lines.end();it++)
      {
      char *text = (char *) it->c_str();
      if (it==lines.begin())
        {
        glRasterPos2i(x,y);
        font->print(text,it->length());
        //stroke_output(x,y,fsHeader,text);
        y -= fsHeader+fsSpace;
        }
      else
        {
        //stroke_output(x,y,fsBody,text);
        glRasterPos2i(x,y);
        font->print(text,it->length());
        y -= fsBody+fsSpace;
        }
      }
    }

  glPopAttrib();
}

void GLFadingWindow::onDraw() {
  if (isReset)
    {
    tLastChange = clock();
    tView = tLastChange + 2000 + lines.size() * 1000;
    tFade = tView + 3000;
    isReset = false;
    }
  else
    {
    glPushMatrix();
    if (clock() < tView)
      {
      color.fv[3] = 1.0f;
      GLOutputWindow::onDraw();
      }
    else if (clock() < tFade)
      {
      float alpha = (tFade-clock())*1.0f / (tFade-tView);
      color.fv[3] = alpha;
      GLOutputWindow::onDraw();
      }
    glPopMatrix();
    }
}

void GLFadingWindow::clear() {
  GLOutputWindow::clear();
  isReset = true;
}

GLFadingWindow::GLFadingWindow(GLFont *font,GLColor color) : GLOutputWindow(font,color) {
  isReset = true;
}

list<ModeIcon *> ModeIcon::allIcons;

ModeIcon::ModeIcon(string img1,string img2,string modeName,int icIndex)
: GLCommandIcon(img1.c_str(),img2.c_str(),-144,144*(icIndex+1),128,128), mode(modeName) 
{
  allIcons.push_back(this);
}

ModeIcon::~ModeIcon() {
  allIcons.erase(find(allIcons.begin(),allIcons.end(),this));
}

void ModeIcon::onStateChange() {
  if (clickState)
    {
    // Disable all other icons
    list<ModeIcon *>::iterator it;
    for (it=allIcons.begin();it!=allIcons.end();it++)
      {
      ModeIcon *icon = *it;
      if (icon!=this)
        icon->setState(false);
      }

    // Set the current mode
    modeManager->setMode(mode);
    }
  else
    {
    // Set the current mode
    modeManager->setMode("null");
    }
}

ModeManager::ModeManager() {
  activeMode = NULL;
}

void ModeManager::addMode(string name,ModeHandler *handler) {
  modes[name] = handler;
}

string ModeManager::getActiveMode() {
  return activeModeName;
}

void ModeManager::setMode(string name) {
  // Uninstall the current mode
  if (activeMode)
    {
    activeMode->stop();
    GLDisplayDriver::removeListener(activeMode,GLDisplayDriver::MOTION | GLDisplayDriver::BUTTON | GLDisplayDriver::KEYS | GLDisplayDriver::SPECIAL);
    }

  activeMode = modes[name];
  activeModeName = name;

  if (activeMode)
    {
    try
      {
      activeMode->start();
      GLDisplayDriver::addListener(activeMode,GLDisplayDriver::MOTION | GLDisplayDriver::BUTTON | GLDisplayDriver::KEYS | GLDisplayDriver::SPECIAL);
      }
    catch (string exc)
      {
      activeMode = NULL;
      activeModeName = "null";
      rndFadeStatus.addLine(exc.c_str()); 
      }
    }

  rndFadeStatus.addLine(string("Entering mode: " + name).c_str());
}

// Draw the border and the header of the window, pass the rest on to the clients
void MoveableWindow::onDraw() {

  glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);

  glPushMatrix();
  glTranslated(position[0],position[1],position[2]);

  glDisable(GL_LIGHTING);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // Draw the background of the window
  glColor4d(0.1,0.3,0.4,0.2);
  glBegin(GL_QUADS);
  glVertex3d(0,0,0);
  glVertex3d(size[0],0,0);
  glVertex3d(size[0],size[1],0);
  glVertex3d(0,size[1],0);
  glVertex3d(0,0,0);
  glEnd();

  // Border color depends on the focus of the window
  if (hasFocus())
    glColor3d(0.3,0.8,0.9);
  else
    glColor3d(0.2,0.2,0.2);

  glEnable(GL_LINE_SMOOTH);

  // Draw the border
  glBegin(GL_LINE_STRIP);
  glVertex3d(0,0,0);
  glVertex3d(size[0],0,0);
  glVertex3d(size[0],size[1],0);
  glVertex3d(0,size[1],0);
  glVertex3d(0,0,0);
  glEnd();

  // Draw the underline
  glBegin(GL_LINES);
  glVertex3d(0,szHat+1,0);
  glVertex3d(size[0],szHat+1,0);
  glEnd();

  // Draw the title
  glRasterPos2d(position[0] + 5,position[1] + szHat - 2);   
  fntCourier12->print(title.c_str(),title.length());

  // Translate more, to leave a corner
  glTranslated(3,szHat+3,0);
  drawClient((int)(size[0] - 6),(int)(size[1] - (szHat + 6)));

  // Pop matrix
  glPopMatrix();
  glPopAttrib();
}

bool MoveableWindow::handleKeys(unsigned char key,int x,int y) {
  return false;
}

bool MoveableWindow::handleButton(int button,int state,int x,int y) {
  return true;
}

bool MoveableWindow::handleMotion(int x,int y) {
  return true;
}

MoveableWindow::MoveableWindow(int x,int y,int w,int h) {
  position[0] = x;
  position[1] = y;
  position[2] = 0;
  size[0] = w;
  size[1] = h;
  size[2] = 0;
}

string getTextFromClipboard() {
  string clip = "";
#ifdef _MSC_VER_
  if (OpenClipboard(NULL))
    {
    HGLOBAL hglb = GetClipboardData(CF_TEXT);
    char *text = (char *)hglb;  
    clip += text;
    CloseClipboard();
    }
#endif
  return clip;
}


void CommandShell::print(const char *text,int length) {
  for (int i=0;i<length;i++)
    {
    if (lines.back().length() >= 80 || text[i] == '\n' || text[i] == '\r')
      {
      lines.push_back("");
      }
    else
      {
      char c = text[i];
      lines.back().append(&c,1);
      }
    }
  needRedraw = true;
}

void CommandShell::onDraw() {

  static clock_t tBlink = clock();

  glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);

  glPushMatrix();
  glTranslated(pos[0],pos[1],0);

  glDisable(GL_LIGHTING);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glColor4d(0.1,0.3,0.4,0.2);

  glBegin(GL_QUADS);
  glVertex3d(0,0,0);
  glVertex3d(size[0],0,0);
  glVertex3d(size[0],size[1],0);
  glVertex3d(0,size[1],0);
  glVertex3d(0,0,0);
  glEnd();

  if (hasFocus)
    glColor3d(0.3,0.8,0.9);
  else
    glColor3d(0.2,0.2,0.2);

  glEnable(GL_LINE_SMOOTH);

  glBegin(GL_LINE_STRIP);
  glVertex3d(0,0,0);
  glVertex3d(size[0],0,0);
  glVertex3d(size[0],size[1],0);
  glVertex3d(0,size[1],0);
  glVertex3d(0,0,0);
  glEnd();

  string lastLine = ">> " + prompt;

  if (needRedraw)
    {
    outwin.clear();

    int ls = (int)lines.size();
    int iStart = ls > rows ? ls - rows : 0;     
    for (int i=iStart;i<ls;i++)
      {
      outwin.addLine(lines[i].c_str());
      }
    outwin.addLine(lastLine.c_str());

    needRedraw = false;
    }

  outwin.onDraw();

  // Draw cursor (maybe)
  if (clock()-tBlink > 1000)
    {
    tBlink = clock();
    }
  else if (clock()-tBlink > 250)
    {
    glBegin(GL_LINES);
    int cp = outwin.font->getStringWidth(lastLine.substr(0,cursor+3).c_str());
    glVertex2d(cp+7,32);
    glVertex2d(cp+7,16);
    glEnd();
    }

  glPopMatrix();

  glPopAttrib();

  glutPostRedisplay();
}

CommandShell::CommandShell() : outwin(fntCourier12,GLColor(0.3,0.8,0.9)) {
  hasFocus = true;
  needRedraw = true;

  // Load the position and the size of the command shell from settings registry
  pos[0] = settingsUI->getIntValue("commandShell.position.x",0);
  pos[1] = settingsUI->getIntValue("commandShell.position.y",600);
  size[0] = settingsUI->getIntValue("commandShell.size.x",800);
  size[1] = settingsUI->getIntValue("commandShell.size.y",360);

  // Compute the number of rows
  rows = (int)(size[1] / 14);
  cols = 100;

  cursor = 0;
  idxHistory = 0;
  lines.push_back("Command Shell v 0.0.1");
  lines.push_back("");

  // Try to load the history list
  int rhSize = settingsUI->getIntValue("commandShell.historySize",0);
  int rhEnd = settingsUI->getIntValue("commandShell.historyEnd",-1);
  for (int i=1;i<=rhSize;i++)
    {
    int idx = (rhEnd+i) % rhSize;
    history.push_back(settingsUI->getStringValue("commandShell.history[%d]","null",idx));
    idxHistory++;
    }
}

bool CommandShell::handleButton(int button,int state,int x,int y) {
  y = GLDisplayDriver::height - y;

  // Gain mouse focus if the mouse is clicked inside
  if (button == GLUT_LEFT_BUTTON)
    {
    // Once the mouse is up, no longer care about motion events
    if (state == GLUT_UP)
      GLDisplayDriver::removeListener(this,GLDisplayDriver::MOTION);

    if (pos[0] <= x && x < pos[0] + size[0] && pos[1] <= y && y < pos[1] + size[1])
      {
      if (state == GLUT_UP)
        {
        hasFocus = true;
        return true;
        }
      else
        {
        // We want to get motion events
        GLDisplayDriver::addListener(this,GLDisplayDriver::MOTION);
        return true;
        }
      }
    else
      {
      hasFocus = false;
      }
    }
  return false;
}

bool CommandShell::handleMotion(int x,int y) {
  if (GLDisplayDriver::dragButton == GLUT_LEFT_BUTTON)
    {
    // Get the deltas for the dragging
    int dx = x - GLDisplayDriver::dragX;
    int dy = GLDisplayDriver::dragY - y;

    // Adjust the screen position
    pos[0] += dx;
    pos[1] += dy;

    // Save position to settings
    settingsUI->setDoubleValue("commandShell.position.x",pos[0]);
    settingsUI->setDoubleValue("commandShell.position.y",pos[1]);

    // Need to redraw
    glutPostRedisplay();
    return true;
    }

  return false;
}

bool CommandShell::handleSpecial(int key,int x,int y) {
  if (hasFocus)
    {
    if (key == GLUT_KEY_UP)
      {
      if (idxHistory > 0)
        {
        idxHistory--;
        prompt = history[idxHistory];
        cursor = prompt.length();               
        }
      }
    else if (key == GLUT_KEY_DOWN)
      {
      if (idxHistory+1 < history.size())
        {
        idxHistory++;
        prompt = history[idxHistory];
        cursor = prompt.length();
        }
      }
    else if (key == GLUT_KEY_LEFT)
      {
      if (cursor > 0)
        cursor--;
      }
    else if (key == GLUT_KEY_RIGHT)
      {
      if (cursor < prompt.length())
        cursor++;
      }
    else if (key == GLUT_KEY_HOME)
      {
      cursor = 0;
      }
    else if (key == GLUT_KEY_END)
      {
      cursor = prompt.length();
      }
    else
      {
      return false;
      }

    // Update the screen
    needRedraw = true;
    glutPostRedisplay();

    return true;
    }
  return false;
}

bool CommandShell::handleKeys(unsigned char key,int x,int y) {
  if (hasFocus)
    {
    // Check modifiers - modified keys are passed on down
    if (glutGetModifiers() & GLUT_ACTIVE_CTRL || glutGetModifiers() & GLUT_ACTIVE_ALT)
      {
      if (key == 22)
        {
        string paste = getTextFromClipboard();
        prompt.insert(cursor,paste);
        cursor+=paste.length();
        }
      else
        {
        return false;
        }
      }

    // If a regular key is pressed 
    else if (key >= 32 && key <= 128)
      {
      prompt.insert(prompt.begin()+cursor,key);
      cout << prompt << endl;
      // prompt += key;
      cursor++;
      }

    else if (key=='\r' || key=='\n')
      {
      // print("\n");
      lines.push_back(prompt.c_str());
      lines.push_back("");

      // Append to history
      if (history.size() == 0 || prompt != history.back())
        {
        history.push_back(prompt.c_str());
        idxHistory = history.size();

        // Save the registry's history
        int rhSize = settingsUI->getIntValue("commandShell.historySize",0);
        int rhEnd = settingsUI->getIntValue("commandShell.historyEnd",-1);
        if (rhEnd+1 < REGHISTSIZE)
          {
          settingsUI->setIntValue("commandShell.historyEnd",++rhEnd);
          if (rhEnd+1 > rhSize)
            settingsUI->setIntValue("commandShell.historySize",rhEnd+1);
          settingsUI->setStringValue("commandShell.history[%d]",prompt.c_str(),rhEnd);
          }
        else
          {
          settingsUI->setIntValue("commandShell.historyEnd",0);
          settingsUI->setStringValue("commandShell.history[0]",prompt.c_str());
          }
        }

      // Process the command from command line
      ostringstream oss;
      processCommand(prompt.c_str(),oss);
      string result = oss.str();

      if (result.length() > 0)
        {
        // print("\n");
        print(result.c_str(),result.length());
        }

      prompt = "";
      cursor=0;           
      }

    else if (key=='\b')
      {
      if (cursor > 0)
        {
        prompt = prompt.substr(0,cursor-1) + prompt.substr(cursor);
        cursor--;
        }
      }

    else
      {
      // Not a key we can handle!
      return false;
      }


    needRedraw = true;
    return true;

    glutPostRedisplay();        
    }
  return false;
}

GlobalKeyBinding::GlobalKeyBinding(int key,int state,string command) {
  this->command = command;
  this->key = key;
  this->state = state;

  // Generate Description
  if (state & GLUT_ACTIVE_CTRL)
    desc += "Ctrl-";
  if (state & GLUT_ACTIVE_ALT)
    desc += "Alt-";
  if (state & GLUT_ACTIVE_SHIFT)
    desc += "Shift-";
  if (key > 32 && key < 128)
    desc += key;
  else if (key >= GLUT_KEY_F1 && key < GLUT_KEY_F10)
    {
    desc += "F";
    desc += (1 + key - GLUT_KEY_F1);
    }
  else if (key >= GLUT_KEY_F10 && key <= GLUT_KEY_F12)
    {
    desc += "F1";
    desc += (key - GLUT_KEY_F10);
    }
  else if (key == GLUT_KEY_DOWN)
    desc += "Down Arrow";
  else if (key == GLUT_KEY_UP)
    desc += "Up Arrow";
  else if (key == GLUT_KEY_LEFT)
    desc += "Left Arrow";
  else if (key == GLUT_KEY_RIGHT)
    desc += "Right Arrow";
  else if (key == GLUT_KEY_HOME)
    desc += "Home";
  else if (key == GLUT_KEY_END)
    desc += "End";
  else if (key == GLUT_KEY_PAGE_DOWN)
    desc += "PgDn";
  else if (key == GLUT_KEY_PAGE_UP)
    desc += "PgUp";
  else if (key == GLUT_KEY_INSERT)
    desc += "Insert";
}

void executeCommand(string command) {
  ostringstream oss;
  processCommand(command.c_str(),oss);

  if(flagInteractive)
    {
    // Display in the command shell
    if (commandShell->isVisible())
      {
      commandShell->print(command.c_str(),command.length());
      commandShell->print("\n",1);
      commandShell->print(oss.str().c_str(),oss.str().length());
      }
    else
      {
      // Show command to user
      rndFadeStatus.clear();
      rndFadeStatus.addLine(command.c_str());
      rndFadeStatus.addLine(oss.str().c_str());
      }
    }
  else
    {
    cout << command << endl;
    cout << oss.str() << endl;
    }
}

void echoMessage(string message) {
  rndFadeStatus.clear();
  rndFadeStatus.addLine(message.c_str());
}

void GlobalEventHandler::addKeyBinding(int key,int state,string command) {
  bind.push_back(GlobalKeyBinding(key,state,command));
}

// Handle global key binding
bool GlobalEventHandler::handleKeys(unsigned char key,int x,int y) {
  int state = glutGetModifiers();

  // Find matching key
  list<GlobalKeyBinding>::iterator it;
  for (it = bind.begin();it!=bind.end();it++)
    {
    if (it->key == key && it->state == state)
      {
      // Execute command
      executeCommand(it->command.c_str());
      return true;
      }
    }

  return false;
}

bool GlobalEventHandler::handleSpecial(int key,int x,int y) {
  if (key==GLUT_KEY_F1)
    {
    rndFadeStatus.clear();
    rndFadeStatus.addLine("Main Keys:");
    rndFadeStatus.addLine("   F2              Bring up the command console");
    rndFadeStatus.addLine("   Esc             Exit the program");
    rndFadeStatus.addLine("");
    rndFadeStatus.addLine("Key-Command Bindings:");
    // Find matching key
    list<GlobalKeyBinding>::iterator it;
    for (it = bind.begin();it!=bind.end();it++)
      {
      rndFadeStatus.addLine("   %-16s%s",it->desc.c_str(),it->command.c_str());
      }
    }
  if (key==GLUT_KEY_F2)
    {
    if (commandShell->isVisible())
      {
      commandShell->setVisible(false);
      GLDisplayDriver::removeListener(commandShell,GLDisplayDriver::KEYS | GLDisplayDriver::BUTTON | GLDisplayDriver::SPECIAL | GLDisplayDriver::MOTION); 
      }
    else
      {
      commandShell->setVisible(true);
      GLDisplayDriver::addListener(commandShell,GLDisplayDriver::KEYS | GLDisplayDriver::BUTTON | GLDisplayDriver::SPECIAL); 
      }
    glutPostRedisplay();
    return true;
    }
  else
    {
    int state = glutGetModifiers();

    // Find matching key
    list<GlobalKeyBinding>::iterator it;
    for (it = bind.begin();it!=bind.end();it++)
      {
      if (it->key == key && it->state == state)
        {
        if (it->key == key && it->state == state)
          {
          executeCommand(it->command.c_str());
          return true;
          }
        }
      }
    }
  return false;
}

/**
 * This subroutine initializes the modes and mode icons
 */
void initModes() {
  // Create modes
  modeManager->addMode("surfselect",new SurfaceSelectionMH());
  modeManager->addMode("controledt",new ControlPointEditMH());
  modeManager->addMode("optimize",new OptimizationMH());
  modeManager->addMode("rigid",new RigidMatchMode());
  modeManager->addMode("null",NULL);

  // Create mode icons
  GLDisplayDriver::addCommandIcon(new ModeIcon("robot.tga","robot_act.tga","surfselect",0));
  GLDisplayDriver::addCommandIcon(new ModeIcon("woild.tga","woild_act.tga","controledt",1));
}

extern void testMultiply();

class SliceRenderer : public GLDisplayListRenderer
  {
public:
  int sliceNo;
  int dim;

  SliceRenderer();
  void build();
  void onDraw();
  };

SliceRenderer::SliceRenderer() {
  setVisible(false);
  sliceNo = 0;
  dim = 0;
}

void SliceRenderer::onDraw() {
  if (imaging.getImage() == NULL)
    return;

  static clock_t t = clock();
  if (clock() - t > 250)
    {
    sliceNo = (sliceNo+1) % imaging.getImage()->size(dim);
    invalidate();
    t = clock_t();
    }

  GLDisplayListRenderer::onDraw();
}

void SliceRenderer::build() {
  AbstractImage3D *image = imaging.getImage();

  if (!image)
    return;

  int res;
  float trimS,trimT;

  unsigned char *slice = image->getSlice(dim,sliceNo,res,trimS,trimT);

  glRasterPos2d(0,0);
  glDrawPixels(image->size((dim+1)%3),image->size((dim+2)%3),GL_LUMINANCE,GL_UNSIGNED_BYTE,slice);
}

// SliceRenderer sliceRenderer;



/********************************
 Image Box Renderer
 *******************************/
class ImageBoxRenderer : public GLDisplayListRenderer
  {
private:
  unsigned int texName;
  float mpRed[256],mpBlue[256],mpGreen[256],mpAlpha[256];

  void buildColorMap();
public:
  int sliceNo;
  int dim;

  ImageBoxRenderer();
  void onDraw();
  void build();

  void nextSlice() {
    AbstractImage3D *image = imaging.getImage();
    if (image)
      {
      sliceNo = (sliceNo+1) % image->size(dim);
      invalidate();
      }
  }

  void lastSlice() {
    AbstractImage3D *image = imaging.getImage();
    if (image)
      {
      sliceNo = sliceNo<=0 ? image->size(dim)-1 : sliceNo-1;
      invalidate();
      }
  }

  // Whether the image is rolling or not
  bool sliceMoving;

  };

void ImageBoxRenderer::buildColorMap() {
  // 0 and 1 are special
  mpBlue[0] = mpGreen[0] = mpRed[0] = 1.0f;
  mpBlue[1] = mpGreen[1] = mpRed[1] = 1.0f;

  for (int i=2;i<256;i++)
    {
    mpRed[i] = (255.0f - i) / 253.0f;
    mpGreen[i] = (255 - i) / 253.0f;
    mpBlue[i] = 0.25f;
    }

  for (int j=0;j<256;j++)
    mpAlpha[j] = 1.0f;
}

void ImageBoxRenderer::onDraw() {
  AbstractImage3D *image = imaging.getImage();    
  if (!image) return;

  static clock_t t = clock();
  if (clock() - t > 250 && sliceMoving)
    {
    sliceNo = (sliceNo+1) % image->size(dim);
    invalidate();
    t = clock_t();
    }

  GLDisplayListRenderer::onDraw();    
}

ImageBoxRenderer::ImageBoxRenderer() {
  setVisible(false);
  sliceNo = 0;
  dim = 0;    
  texName = 0;
  buildColorMap();
  sliceMoving = true;
}

void ImageBoxRenderer::build() {
  AbstractImage3D *image = imaging.getImage();    
  if (!image) return;

  glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_PIXEL_MODE_BIT);
  glDisable(GL_LIGHTING);

  glPushMatrix();
  glLoadIdentity();

  // Draw the surrounding cube
  SMLMatrix4f IST = image->IS.transpose();

  glTranslated(-0.5,-0.5,-0.5);
  glMultMatrixf(IST.data_block());

  int d1 = (dim+1)%3;
  int d2 = (dim+2)%3;

  // Box coordinate computation
  int x[3][2][2];
  for (int q=0;q<2;q++)
    {
    for (int p=0;p<2;p++)
      {
      x[dim][q][p] = sliceNo;
      x[d1][q][p] = q * image->size(d1); 
      x[d2][q][p] = p * image->size(d2); 
      }
    }

  // Get the slice data
  int sliceRes;
  float trimS,trimT;
  unsigned char *slice = image->getSlice(dim,sliceNo,sliceRes,trimS,trimT);

  //for(int i=0;i<sliceRes*sliceRes;i++)
  //  slice[i] *= 256;

  // Allow textures
  glEnable(GL_TEXTURE_2D);
  // glEnable(GL_LIGHTING);

  // Create the texture
  if (texName == 0)
    glGenTextures(1, &texName);
  glBindTexture(GL_TEXTURE_2D,texName);
  glPixelStorei(GL_UNPACK_ALIGNMENT,1);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER, GL_NEAREST);

  // glPixelTransferi(GL_MAP_COLOR,GL_TRUE);
  glPixelMapfv(GL_PIXEL_MAP_I_TO_R,256,mpRed);
  glPixelMapfv(GL_PIXEL_MAP_I_TO_G,256,mpGreen);
  glPixelMapfv(GL_PIXEL_MAP_I_TO_B,256,mpBlue);
  glPixelMapfv(GL_PIXEL_MAP_I_TO_A,256,mpAlpha);

  // glTexImage2D(GL_TEXTURE_2D,0,4,sliceRes,sliceRes,0,GL_COLOR_INDEX,GL_UNSIGNED_BYTE,slice);
  if (imaging.getType() == imaging.DISTANCE)
    {
    int sz = 3*sliceRes*sliceRes;
    unsigned char *clrSlice = new unsigned char[sz];
    for (int i=0;i<sz;)
      {
      char ch = *((char *)slice);
      unsigned char r,g,b;
      if (ch < 0)
        {
        r = g = -ch;
        b = r+g;
        }
      else if (ch > 0)
        {
        g = b = ch;
        r = g + b;
        }
      else
        {
        r = g = b = 255;
        }

      clrSlice[i++] = r;
      clrSlice[i++] = g;
      clrSlice[i++] = b;

      slice++;
      }
    //
    //}
    glTexImage2D(GL_TEXTURE_2D,0,4,sliceRes,sliceRes,0,GL_RGB,GL_UNSIGNED_BYTE,clrSlice);

    delete clrSlice;
    }
  else
    {
    glTexImage2D(GL_TEXTURE_2D,0,4,sliceRes,sliceRes,0,GL_LUMINANCE,GL_UNSIGNED_BYTE,slice);
    }

  // Draw the appropriate slice
  glColor3d(1,1,1);
  glBegin(GL_QUADS);

  glTexCoord2d(0,0);
  glVertex3d(x[0][0][0],x[1][0][0],x[2][0][0]);
  glTexCoord2d(0,trimT);
  glVertex3d(x[0][0][1],x[1][0][1],x[2][0][1]);
  glTexCoord2d(trimS,trimT);
  glVertex3d(x[0][1][1],x[1][1][1],x[2][1][1]);
  glTexCoord2d(trimS,0);
  glVertex3d(x[0][1][0],x[1][1][0],x[2][1][0]);

  glEnd();

  // Draw the box around the image
  glDisable(GL_TEXTURE_2D);
  glScaled(image->size(0),image->size(1),image->size(2));
  glTranslated(0.5,0.5,0.5);
  glColor3d(0.15,0.15,0.35);
  glutWireCube(1);

  glPopMatrix();
  glPopAttrib();
}

// ImageBoxRenderer ibRenderer;
ImageBoxRenderer sliceRenderer;

/*
class LevelSetRenderer : public GLDisplayListRenderer {
private:
public:
    LevelSetRenderer();
    void build();
};

LevelSetRenderer::LevelSetRenderer() {
    setVisible(false);
}

void LevelSetRenderer::build() {
    
    if(imaging.getType() != ImagingSystem::DISTANCE)
        return;

    DistanceTransform *distMap = (DistanceTransform*)imaging.getImage();

    glPushAttrib(GL_LIGHTING | GL_POINT_BIT);
    glDisable(GL_LIGHTING);

    glColor3d(0.8,0.4,0.25);

    // Set coordinate transform
    glPushMatrix();
    
    glTranslated(-0.5,-0.5,-0.5);
    SMLMatrix4f IST = distMap->IS;
    IST.Transpose();

    glMultMatrixf(IST.GetData());
    glPointSize(1.0);

    glBegin(GL_POINTS);
    
    // Go through all the voxels and add to a point list
    for(int cz=0;cz<distMap->cube.size(2);cz++) {
        for(int cy=0;cy<distMap->cube.size(1);cy++) {
            for(int cx=0;cx<distMap->cube.size(0);cx++) {
                DataCube<short> &cube = distMap->cube(cx,cy,cz);
                int bx = cx * cube.size(2);
                int by = cy * cube.size(1);
                int bz = cz * cube.size(2);

                for(int z=0;z<cube.size(2);z+=2) {
                    for(int y=0;y<cube.size(1);y+=2) {
                        for(int x=0;x<cube.size(0);x+=2) {
                            if(cube(x,y,z) < 2)
                                glVertex3i(bx+x,by+y,bz+z);
                        }
                    }
                }
            }
        }
    }

    glEnd();

    glPopMatrix();
    glPopAttrib();
}
*/

#include "implicit.h"

class LevelSetRenderer : public GLDisplayListRenderer
  {
private:
  static AbstractImage3D *img;
  static double impfunDistance(double x,double y,double z);
  static double impfunBinary(double x,double y,double z);
  static int  triCallback(int v1,int v2,int v3,IMP_VERTICES vtx);
  static LevelSetRenderer *current;

  // Vertex and triangle data
  vector<int> tri;
  IMP_VERTICES vtx;

public:
  LevelSetRenderer();
  void buildList();
  void computeLevelSet();
  };

AbstractImage3D *LevelSetRenderer::img = NULL;
LevelSetRenderer *LevelSetRenderer::current = NULL;

double LevelSetRenderer::impfunDistance(double x,double y,double z) {
  return img->interpolateVoxel(x,y,z);
}

double LevelSetRenderer::impfunBinary(double x,double y,double z) {
  return img->interpolateVoxel(x,y,z) - 128;
}

int LevelSetRenderer::triCallback(int v1,int v2,int v3,IMP_VERTICES vtx) {
  current->vtx = vtx;
  current->tri.push_back(v1);
  current->tri.push_back(v2);
  current->tri.push_back(v3);

  //SMLVec3f x(vtx.ptr[v1].position.x,vtx.ptr[v1].position.y,vtx.ptr[v1].position.z),ix;
  //img->SI.TransformPoint(x,ix);

  // cout << ix.x << "\t" <<  ix.y << "\t" <<  ix.z << "\t" << img->interpolateVoxel(x.x,x.y,x.z) << endl;

  // return (vtx.count > 10000) ? 0 : 1;
  return 1;
}

void LevelSetRenderer::computeLevelSet() {
  img = imaging.getImage();
  SMLVec3f start;
  double (*function)(double,double,double) = NULL;

  // Clear the triangle array
  tri.clear();
  vtx.ptr = NULL;
  vtx.count = 0;

  if (img)
    {
    current = this;
    if (imaging.getType() == imaging.DISTANCE)
      {
      // Find a starting point for image match
      for (int z=0;z<img->size(2);z++)
        {
        for (int y=0;y<img->size(1);y++)
          {
          for (int x=0;x<img->size(0)-1;x++)
            {
            if (img->getVoxelNBCFloat(x,y,z) + img->getVoxelNBCFloat(x+1,y,z) == 0.0)
              {
              SMLVec3f istart(x+0.5,y,z);
              start = TransformPoint(img->IS,istart);

              cout << "Starting at voxel " << x << " " << y << " " << z << endl;
              cout << img->interpolateVoxel(start[0],start[1],start[2]) << endl;

              function = impfunDistance;
              z = img->size(2);
              y = img->size(1);
              x = img->size(0);

              }
            }
          }
        }
      }
    else if (imaging.getType() == imaging.BINARY)
      {
      // Find a starting point for image match
      for (int z=0;z<img->size(2);z++)
        {
        for (int y=0;y<img->size(1);y++)
          {
          for (int x=0;x<img->size(0)-1;x++)
            {
            if (fabs((double)(img->getVoxelNBC(x+1,y,z) - img->getVoxelNBC(x,y,z))) > 128.0)
              {
              SMLVec3f istart(x,y,z);
              start = TransformPoint(img->IS,istart);
              function = impfunBinary;
              break;
              }
            }
          }
        }
      }

    if (function)
      {
      // Get the size of the cube
      SMLVec3f ivox(1,1,1);
      SMLVec3f vox = TransformVector(img->IS,ivox);
      // float dmin = 2.7 * min(min(vox.x,vox.y),vox.z);
      float dmin = 1.0 / 256.0f;

      // Run the polygonizer
      polygonize(function,dmin,(int)(1.0 / dmin),start[0],start[1],start[2],triCallback,NOTET);
      }
    }
}

LevelSetRenderer::LevelSetRenderer() {
  setVisible(false);
}

void LevelSetRenderer::buildList() {
  // Compute the level set
  computeLevelSet();

  // Must have some triangles
  if (tri.size() == 0)
    return;

  // Set client state
  glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);
  glEnableClientState(GL_VERTEX_ARRAY);

  // Input the vertex arrays
  glVertexPointer(3,GL_DOUBLE,sizeof(IMP_VERTEX),&vtx.ptr->position);

  // Start the magic
  glNewList(dl,GL_COMPILE);

  // Push attributes
  glPushAttrib(GL_LIGHTING | GL_POLYGON_BIT | GL_TRANSFORM_BIT);
  glDisable(GL_LIGHTING);

  // Move back
  glPushMatrix();
  glTranslated(-0.5,-0.5,-0.5);

  // Set the properties
  glColor3d(0.8,0.4,0.25);
  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

  // Draw the elements
  glBegin(GL_TRIANGLES);
  for (int i=0;i<tri.size();i++)
    glArrayElement(tri[i]);
  glEnd();

  // Pop matrix
  glPopMatrix();

  // Pop the attributes
  glPopAttrib();

  // End the list
  glEndList();

  // Pop client att
  glPopClientAttrib();
}

/*  
    if(imaging.getType() != ImagingSystem::DISTANCE)
        return;

    DistanceTransform *distMap = (DistanceTransform*)imaging.getImage();

    glPushAttrib(GL_LIGHTING | GL_POINT_BIT);
    glDisable(GL_LIGHTING);

    glColor3d(0.8,0.4,0.25);

    // Set coordinate transform
    glPushMatrix();
    
    glTranslated(-0.5,-0.5,-0.5);
    SMLMatrix4f IST = distMap->IS;
    IST.Transpose();

    glMultMatrixf(IST.GetData());
    glPointSize(1.0);

    glBegin(GL_POINTS);
    
    // Go through all the voxels and add to a point list
    for(int cz=0;cz<distMap->cube.size(2);cz++) {
        for(int cy=0;cy<distMap->cube.size(1);cy++) {
            for(int cx=0;cx<distMap->cube.size(0);cx++) {
                DataCube<short> &cube = distMap->cube(cx,cy,cz);
                int bx = cx * cube.size(2);
                int by = cy * cube.size(1);
                int bz = cz * cube.size(2);

                for(int z=0;z<cube.size(2);z+=2) {
                    for(int y=0;y<cube.size(1);y+=2) {
                        for(int x=0;x<cube.size(0);x+=2) {
                            if(cube(x,y,z) < 2)
                                glVertex3i(bx+x,by+y,bz+z);
                        }
                    }
                }
            }
        }
    }

    glEnd();

    glPopMatrix();
    glPopAttrib();
}   */



LevelSetRenderer levelSetRenderer;


class EdgeNetRenderer : public GLDisplayListRenderer
  {
public:
  void loadNetData(const char *netFile, SMLVec3f voxelSize);
  void build();

private:
  struct Edge
    {
    int v[2];
    };

  // Collection of edges
  vector<Edge> e;
  vector<SMLVec3f> v;
  };

void EdgeNetRenderer::build() {
  if (e.size())
    {
    glPushAttrib(GL_LIGHTING_BIT);

    glDisable(GL_LIGHTING);
    glColor3d(1,1,1);

    glBegin(GL_LINES);
    for (int i=0;i<e.size();i++)
      {
      glVertex(v[e[i].v[0]]); 
      glVertex(v[e[i].v[1]]); 
      }
    glEnd();

    glPopAttrib();
    }
}

void EdgeNetRenderer::loadNetData(const char *netFile, SMLVec3f voxelSize) {
  FILE *f = fopen(netFile,"rt");
  char buffer[256];

  SMLVec3f mx,mn;

  while (fgets(buffer,256,f))
    {
    vector<string> tmp;
    char *p = strtok(buffer,"{},\t\n\r ");
    while (p)
      {
      string s = p;
      tmp.push_back(s);
      p = strtok(NULL,"{},\t\n\r");

      }
    if (tmp.size() > 3)
      {
      SMLVec3f x;
      x[0] = atof(tmp[0].c_str()) * voxelSize[0];
      x[1] = atof(tmp[1].c_str()) * voxelSize[1];
      x[2] = atof(tmp[2].c_str()) * voxelSize[2];

      if (v.size())
        {
        mx[0] = vnl_math_max(x[0],mx[0]);
        mx[1] = vnl_math_max(x[1],mx[1]);
        mx[2] = vnl_math_max(x[2],mx[2]);
        mn[0] = vnl_math_min(x[0],mn[0]);
        mn[1] = vnl_math_min(x[1],mn[1]);
        mn[2] = vnl_math_min(x[2],mn[2]);
        }
      else
        {
        mx = mn = x;
        }

      for (int i=3;i<tmp.size();i++)
        {
        int ev = atoi(tmp[i].c_str());
        if (ev > v.size())
          {
          Edge edge;
          edge.v[0] = v.size();
          edge.v[1] = ev;
          e.push_back(edge);
          }
        }

      v.push_back(x);
      }
    }

  fclose(f);

  // Scale and center the object
  SMLVec3f span = mx - mn;
  SMLVec3f center = mn + span * 0.5f;

  float maxSpan = vnl_math_max(vnl_math_max(span[0],span[1]),span[2]);

  for (int i=0;i<v.size();i++)
    {
    v[i] = (v[i] - center) / (maxSpan);
    //v[i].x = (v[i].x - mn.x) / (maxSpan);
    //v[i].y = (v[i].y - mn.y) / (maxSpan);
    //v[i].z = (v[i].z - mn.z) / (maxSpan);
    }

  invalidate();
}


EdgeNetRenderer edgeNetRenderer;

/*************************************************************
 Compute a medial point at u,v
 ************************************************************/
void computeMedialPoint(float u,float v,MedialPoint &M) {

  // Patch number
  SplineGridDefinition *grid = splineData->getGrid();
  int iPatch = grid->patch(0,grid->getIndexForUV(0,u));
  int jPatch = grid->patch(1,grid->getIndexForUV(1,v));

  // Basis weights
  SMLVec4f Wu[4],Wv[4];

  // Compute the weights
  spline->basisJet(0,iPatch+3,u,Wu);
  spline->basisJet(1,jPatch+3,v,Wv);

  // Interpolate the spline up to the second order
  spline->interpolateMedialPoint02(iPatch,jPatch,Wu,Wv,M);

  // Set u,v
  M.u = u;
  M.v = v;
}

GeodesicRenderer::GeodesicRenderer() {
  hasStart = false;
  hasEnd = false;
}

void GeodesicRenderer::build() {
  if (hasStart && hasEnd)
    {

    // Construct geodesic components
    SplineGeoMan sgm(spline,splineData->getGrid());
    EuclideanMetric em(3);  
    SurfaceGeodesic sg;

    // Compute the geodesic
    double geo[(pnum+2)*2];
    sg.geodes(2,3,&sgm,&em,pnum,0.5,1.0e-4,50,uv,geo);

    // Start the GL stuff
    glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
    glDisable(GL_LIGHTING);
    glPushMatrix();
    glTranslated(-0.5,-0.5,-0.5);
    glColor3d(0.3,1.0,0.3);

    // Paint the geodesic
    glBegin(GL_LINE_STRIP);

    for (int i=0;i<2*pnum+4;)
      {
      float u = (float)geo[i++];
      float v = (float)geo[i++];

      MedialPoint M;
      computeMedialPoint(u,v,M);

      glVertex3d(M.F[0],M.F[1],M.F[2]);
      }

    glEnd();

    glColor3d(0,1,0);
    glBegin(GL_POINTS);

    MedialPoint M;
    computeMedialPoint(uv[SUBS2(2,0,0)],uv[SUBS2(2,0,1)],M);
    glVertex3d(M.F[0],M.F[1],M.F[2]);
    computeMedialPoint(uv[SUBS2(2,1,0)],uv[SUBS2(2,1,1)],M);
    glVertex3d(M.F[0],M.F[1],M.F[2]);

    glEnd();


    // End GL stuff
    glPopMatrix();
    glPopAttrib();
    }
}

void GeodesicRenderer::setStart(float u,float v) {
  uv[SUBS2(2,0,0)] = u;
  uv[SUBS2(2,0,1)] = v;
  hasEnd = false;
  hasStart = true;
  invalidate();
}
void GeodesicRenderer::setEnd(float u,float v) {
  uv[SUBS2(2,1,0)] = u;
  uv[SUBS2(2,1,1)] = v;
  hasEnd = true;
  invalidate();
}
bool GeodesicRenderer::needsStart() {
  return(hasStart==false) || (hasEnd==true);
}




/*************************************************************
 Imaging System Code
 *************************************************************/
ImagingSystem::ImagingSystem() {
  type = NONE;
  bim = NULL;
  dt = NULL;
  gray = NULL;
  mVolume = NULL;
  mDistance = NULL;
  mGradient = NULL;
}

ImagingSystem::~ImagingSystem() {
  discard();
}

void ImagingSystem::discard() {
  switch (type)
    {
    case BINARY : 
      delete bim;
      delete mVolume;
      break;
    case DISTANCE :
      delete dt;
      delete mDistance;
      break;
    case GRAYSCALE :
      delete gray;
      delete mGradient;
      break;
    }
  type = NONE;
}

void ImagingSystem::loadImage(const char *fname,ITypes type,IFiles fileType,ostream &out) {
  switch (type)
    {
    case DISTANCE : 
      loadDistance(fname,fileType,out);
      break;
    case BINARY : 
      loadBinary(fname,fileType,out);
      break;
    case GRAYSCALE : 
      loadGrayscale(fname,fileType,out);
      break;
    };
}

// Load different types of files
void ImagingSystem::loadDistance(const char *fname,IFiles ftype,ostream &out) { 
  if (ftype == USE_EXTENSION)
    {
    if (strstr(fname,".dt") || strstr(fname,".cube"))
      loadDistance(fname,CUBE);
    else
      loadDistance(fname,GIPL);
    }
  else
    {
    discard();
    type = DISTANCE;
    dt = new DistanceTransform();
  
    switch (ftype)
      {

      case CUBE :
        dt->loadFromFile(fname);
        break;

      default:
        dt->loadFromITKReadableFile(fname,out);
        break;
      }

    mDistance = new SplineDistanceMatcher(dt);
    }
}

void ImagingSystem::loadBinary(const char *fname,IFiles ftype,ostream &out) {   
  if (ftype == USE_EXTENSION)
    {
    if (strstr(fname,".cube"))
      loadBinary(fname,CUBE);
    else if (strstr(fname,".byu"))
      loadBinary(fname,BYU);
    else 
      loadBinary(fname,GIPL);
    }
  else
    {
    discard();
    type = BINARY;
    bim = new BinaryImage();

    switch (ftype)
      {
      case CUBE :
        bim->loadFromFile(fname);
        break;
      case BYU :
        bim->loadFromBYUFile(fname,256);
        break;
      case GIPL :
        bim->loadFromITKReadableFile(fname,out);
        break;
      default :
        throw "Unknown binary image format";
      }

    mVolume = new SplineBinaryVolumeMatcher(bim);
    }
}

void ImagingSystem::loadGrayscale(const char *fname,IFiles ftype,ostream &out) {    
  if (ftype == USE_EXTENSION)
    {
    if (strstr(fname,".cube"))
      loadGrayscale(fname,CUBE);
    else loadGrayscale(fname,GIPL);
    }
  else
    {
    discard();
    type = GRAYSCALE;
    gray = new GrayImage();

    switch (ftype)
      {
      case CUBE :
        gray->loadFromFile(fname);
        break;
      case GIPL :
        gray->loadFromITKReadableFile(fname,out);
        break;
      default :
        throw "Unknown binary image format";
      }

    mGradient = new GradientImageMatcher(gray);
    }
}

// Convert the image type to another type
void ImagingSystem::toBinary() {
  if (type == GRAYSCALE)
    {
    bim = new BinaryImage();
    bim->loadFromGray(*gray);
    mVolume = new SplineBinaryVolumeMatcher(bim);
    discard();
    type = BINARY;
    }
  else
    {
    throw "Can not convert this image type to binary";
    }
}

// Convert the image type to another type
void ImagingSystem::toGrayscale() {
  if (type == BINARY)
    {
    gray = new GrayImage();
    gray->loadFromBinary(*bim);
    mGradient = new GradientImageMatcher(gray);
    discard();
    type = GRAYSCALE;
    }
  else
    {
    throw "Can not convert this image type to grayscale";
    }
}

void ImagingSystem::toDistance() {
  if (type == BINARY)
    {
    dt = new DistanceTransform();
    dt->loadFromBinaryImage(*bim,6);
    mDistance = new SplineDistanceMatcher(dt);
    discard();
    type = DISTANCE;
    }
  else
    {
    throw "Can not compute distance transform on current image";
    }
}

// Save the image to a file
void ImagingSystem::saveImage(const char *fname,IFiles ftype) {
  if (ftype != CUBE)
    {
    throw "Can not save to this format!";
    }

  if (type == DISTANCE)
    {
    dt->saveToFile(fname);
    }
  else if (type == BINARY)
    {
    bim->saveToFile(fname);
    }
  else if (type == GRAYSCALE)
    {
    gray->saveToFile(fname);
    }
  else
    {
    throw "Can not save current image type!";
    }
}

AbstractImage3D *ImagingSystem::getImage() {
  switch (type)
    {
    case BINARY : return bim;
    case DISTANCE : return dt;
    case GRAYSCALE : return gray;
    }
  return NULL;
}   

SplineImageMatcher *ImagingSystem::getImageMatch() {
  switch (type)
    {
    case BINARY : return mVolume;
    case DISTANCE : return mDistance;
    case GRAYSCALE : return mGradient;
    }
  return NULL;
}

void ImagingSystem::resetImageMatch() {
  switch (type)
    {
    case BINARY : 
      delete mVolume;
      mVolume = new SplineBinaryVolumeMatcher(bim);
      break;
    case DISTANCE : 
      delete mDistance;
      mDistance = new SplineDistanceMatcher(dt);
      break;
    case GRAYSCALE : 
      delete mGradient;
      mGradient = new GradientImageMatcher(gray);
      break;
    }
}

void onImageChange() {
  modeManager->setMode("null");
  sliceRenderer.invalidate();
  levelSetRenderer.invalidate();  
  imaging.resetImageMatch();
}
/*
string makeImageFileName(const char *fname) {
    ostringstream oss;
    oss << settings->getStringValue("dir.images",".");
    oss << "\\";
    oss << fname;
    return oss.str();
}*/

bool tryUnzipFile(char *fname, ostream &err) {
  // Gzip command
  static string cmdUnzip("gzip -d -f -v ");

  // Get fname length
  int len = strlen(fname);

  // If the filename ends with a .gz, decompress the image
  if (strstr(fname,".gz") == fname + (len-3))
    {
    // Execute the unzip command
    string cmd(cmdUnzip + fname);
    system(cmd.c_str());        

    // Chop off the extension
    fname[len-3] = 0;

    // True - means the file was unzipped
    return true;
    }

  // Not a gzip file
  return false;
}

void rezipFile(char *fname) {
  // Gzip command
  static string cmdZip("gzip -f -v ");

  string cmd(cmdZip + fname);
  system(cmd.c_str());
  fname[strlen(fname)] = '.';
}

bool loadImage(const char *fname, ImagingSystem::ITypes type, ostream &err) {
  char fnmod[1024];

  // Load the transform
  try
    {
    // Image has changed
    onImageChange();

    // Copy fname to a modifiable string
    strcpy(fnmod,fname);

    // Unzip file if needed
    bool zip = tryUnzipFile(fnmod,err);

    // Load the image
    try
      {
      imaging.loadImage(fnmod,type,ImagingSystem::USE_EXTENSION,err);
      }
    catch(...)
      {
      err << "Exception when loading image" << endl;
      return false;
      }

    // Rezip if needed
    if (zip)
      rezipFile(fnmod);

    // Success
    return true;
    }
  catch (char *str)
    {
    err << str << endl;
    return false;
    }
}
bool saveImage(const char *fname, ostream &err) {
  // Load the transform
  try
    {
    // imaging.saveImage(makeImageFileName(fname).c_str());
    imaging.saveImage(fname);
    return true;
    }
  catch (char *str)
    {
    err << str << endl;
    return false;
    }
}

bool saveSpline(const char *fname, ostream &err) {
  try
    {
    Registry reg;
    reg.setIntValue("width",spline->dim(0));
    reg.setIntValue("height",spline->dim(1));
    for (int m=0;m<=spline->dim(0);m++)
      {
      for (int n=0;n<=spline->dim(1);n++)
        {
        for (int d=0;d<4;d++)
          {
          reg.setDoubleValue("control[%d][%d][%d]",spline->getControl(m,n,d),m,n,d);
          }
        }
      }
    reg.writeToFile(fname);
    return true;
    }
  catch (RException exc)
    {
    err << "Exception: " << exc << endl;
    return false;
    }
}

// This is the method to destroy all spline associated gook
void destroySpline() {
  // Go to neutral mode
  modeManager->setMode("null");

  // Remove old renderers (safe in this mode)
  GLDisplayDriver::removeRenderer(rndSpline,GLDisplayDriver::TRAN | GLDisplayDriver::NORM);

  // Delete old things
  delete rndSpline;
  delete splineData;
  delete spline;
}

// What do we do when spline pointer is changed???
void replaceSpline(MSpline *newSpline) {

  // Destroy the old spline
  destroySpline();

  // Build and display a spline
  spline = newSpline;
  splineData = new SplineDataCache(spline,RESLN);
  rndSpline = new BSplineRenderer(spline,splineData,RESLN);
  rndSpline->setSurfaceMode(BSplineRenderer::SEETHRU);

  GLDisplayDriver::addRenderer(rndSpline,GLDisplayDriver::TRAN | GLDisplayDriver::NORM);

  uvImage.reset();
}

void newSpline(int m, int n) {
  // Create a new spline
  MSpline *newSpline = new MSpline(m,n);

  // Update the new spline
  replaceSpline(newSpline);
}

bool loadSpline(const char *fname, ostream &err) {
  try
    {
    Registry reg(fname);
    int w = reg.getIntValue("width",4);
    int h = reg.getIntValue("height",4);

    MSpline *newSpline = new MSpline(w,h);

    for (int m=0;m<=newSpline->dim(0);m++)
      {
      for (int n=0;n<=newSpline->dim(1);n++)
        {
        for (int d=0;d<4;d++)
          {
          double ctl = reg.getDoubleValue("control[%d][%d][%d]",0,m,n,d);
          newSpline->setControl(m,n,d,(float)ctl);
          }
        }
      }

    // Update the new spline
    replaceSpline(newSpline);

    return true;
    }
  catch (RException exc)
    {
    err << "Exception: " << exc << endl;
    return false;
    }
}

/** 
 * Scale the whole spline by a range
 */
void scaleSpline(float sx,float sy,float sz,float sr) {
  int i,j;

  // Find the center of mass of the spline
  SMLVec3f Z(0.0f);
  for (i=0;i<=spline->dim(0);i++)
    {
    for (j=0;j<=spline->dim(1);j++)
      {
      Z += spline->getControl(i,j);
      }
    }
  Z /= (i * j);

  // Transform each point
  for (i=0;i<=spline->dim(0);i++)
    {
    for (j=0;j<=spline->dim(1);j++)
      {
      SMLVec4f C = spline->getControl4f(i,j);
      C[0] = (C[0]-Z[0]) * sx + Z[0];
      C[1] = (C[1]-Z[1]) * sy + Z[1];
      C[2] = (C[2]-Z[2]) * sz + Z[2];
      C[3] = C[3] * sr;
      spline->setControl(i,j,C);
      }
    }
}

/**
 Rotate the spline around a vector
 */
void rotateSpline(float vx,float vy,float vz,float angle) {
  int i,j;

  // Find the center of mass of the spline
  SMLVec3f Z(0.0f);
  for (i=0;i<=spline->dim(0);i++)
    {
    for (j=0;j<=spline->dim(1);j++)
      {
      Z += spline->getControl(i,j);
      }
    }
  Z /= (i * j);

  // Create a rotation matrix
  SMLVec3d axis = (angle*M_PI/180) * SMLVec3d(vx,vy,vz).normalize();
  SMLMatrix3f MR = vnl_matops::d2f(vnl_rotation_matrix(axis).transpose());  

  // Transform each point
  for (i=0;i<=spline->dim(0);i++)
    {
    for (j=0;j<=spline->dim(1);j++)
      {
      SMLVec3f C = spline->getControl(i,j);
      SMLVec3f C1 = MR * (C-Z);
      spline->setControl(i,j,C1+Z);
      }
    }

}

void shiftSpline(float dx,float dy,float dz) {
  // Transform each point
  for (int i=0;i<=spline->dim(0);i++)
    {
    for (int j=0;j<=spline->dim(1);j++)
      {
      SMLVec3f C = spline->getControl(i,j);
      spline->setControl(i,j,0,C[0]+dx);
      spline->setControl(i,j,1,C[1]+dy);
      spline->setControl(i,j,2,C[2]+dz);
      }
    }
}



/**
 * Integrate the spline match accross the image
 */
double matchSpline(ostream &out) {
  // Make sure there's an image and a matcher
  SplineImageMatcher *match = imaging.getImageMatch();    
  if (!match)
    return 0;
  match->setSpline(spline,splineData);

  // Check elapsed time
  double tStart = clock();

  // Refresh data
  float lkhd = match->getMatch();

  // Print out match text
  match->printMatchInfo(out);

  // Print out some prior terms as well
  out << "Priors" << endl;

  //ControlPointConstraint cpc;
  //out << "      control penalty    : " << cpc.priorProbabilityFull(spline,splineData) << endl;

  //CurvaturePrior cpr;   
  //out << "      mean curviness     : " << cpr.priorProbabilityFull(spline,splineData) << endl;

  //MRepConstraint mrc(&settings->getSubFolder("optimizer.constraints"));
  //out << "      constraint penalty : " << mrc.priorProbabilityFull(spline,splineData) << endl;

  MRepFastConstraint mrfc(spline);
  out << "       fast m-rep penalty : " << mrfc.priorProbabilityFull(splineData) << endl;

  MRepRadiusConstraint mrrc(spline);
  out << "           radius penalty : " << mrrc.priorProbabilityFull(splineData) << endl;

  MarkovControlPrior mcp(spline);
  out << "          control penalty : " << mcp.priorProbabilityFull(splineData) << endl;

  MarkovCurvaturePrior mcvp(spline);
  out << "        curvature penalty : " << mcvp.priorProbabilityFull(splineData) << endl;

  // Get the elapsed time
  tStart = 1000 * (clock() - tStart) / CLOCKS_PER_SEC;
  out << "      elapsed time (ms)  : " << tStart << endl;

  return lkhd;
}

/*
void testCurvatures(MSpline *spline,int iPatch,int jPatch,int i,int j,ostream &out) {
    float u = splineData->getGrid()->uv(0,iPatch,i);
    float v = splineData->getGrid()->uv(1,jPatch,j);

    SMLVec4f Wu[4],Wv[4],Wue[4],Wve[4],Wuee[4],Wvee[4];
    float eps = 0.001;

    spline->basisJet(0,iPatch+3,u,Wu);
    spline->basisJet(1,jPatch+3,v,Wv);
    spline->basisJet(0,iPatch+3,u+eps,Wue);
    spline->basisJet(1,jPatch+3,v+eps,Wve);
    spline->basisJet(0,iPatch+3,u+eps*2,Wuee);
    spline->basisJet(1,jPatch+3,v+eps*2,Wvee);
    
    // Get the first medial point
    MedialPoint m0 = splineData->patch(iPatch,jPatch)->MP(i,j);
    MedialPoint M0,Mu,Mv,Muv,Muu,Mvv;
    
    spline->interpolateMedialPoint02(iPatch,jPatch,Wu,Wv,M0);
    spline->interpolateMedialPoint02(iPatch,jPatch,Wue,Wv,Mu);
    spline->interpolateMedialPoint02(iPatch,jPatch,Wu,Wve,Mv);
    spline->interpolateMedialPoint02(iPatch,jPatch,Wuee,Wv,Muu);
    spline->interpolateMedialPoint02(iPatch,jPatch,Wu,Wvee,Mvv);
    spline->interpolateMedialPoint02(iPatch,jPatch,Wue,Wve,Muv);
    spline->interpolateBoundaryO2(M0,false);
    spline->interpolateBoundaryO2(Mu,false);
    spline->interpolateBoundaryO2(Mv,false);
    spline->interpolateBoundaryO2(Muv,false);
    spline->interpolateBoundaryO2(Muu,false);
    spline->interpolateBoundaryO2(Mvv,false);

    for(int d=0;d<2;d++) {
        M0.bp[d].Nu = (Mu.bp[d].N - M0.bp[d].N) / eps;
        M0.bp[d].Nv = (Mv.bp[d].N - M0.bp[d].N) / eps;
        M0.bp[d].Xu = (Mu.bp[d].X - M0.bp[d].X) / eps;
        M0.bp[d].Xv = (Mv.bp[d].X - M0.bp[d].X) / eps;

        SMLVec3f Xuu = ((Muu.bp[d].X - Mu.bp[d].X)/eps - (Mu.bp[d].X - M0.bp[d].X)/eps) / eps;
        SMLVec3f Xvv = ((Mvv.bp[d].X - Mv.bp[d].X)/eps - (Mv.bp[d].X - M0.bp[d].X)/eps) / eps;
        SMLVec3f Xuv = ((Muv.bp[d].X - Mu.bp[d].X)/eps - (Mv.bp[d].X - M0.bp[d].X)/eps) / eps;

        float L = Xuu.Dot(M0.bp[d].N);
        float M = Xuv.Dot(M0.bp[d].N);
        float N = Xvv.Dot(M0.bp[d].N);

        // Compare everything
        out << "Nu" << d << ":\t" << (M0.bp[d].Nu-m0.bp[d].Nu).LengthSquared() << endl;
        out << "Nv" << d << ":\t" << (M0.bp[d].Nv-m0.bp[d].Nv).LengthSquared() << endl;
        out << "Xu" << d << ":\t" << (M0.bp[d].Xu-m0.bp[d].Xu).LengthSquared() << endl;
        out << "Xv" << d << ":\t" << (M0.bp[d].Xv-m0.bp[d].Xv).LengthSquared() << endl;
    }
}
*/


void callShellCommand(string command,ostream &sout) {
  string sys = command + " 1> .\\doscmd.mspline.tmp 2> .\\doscmd.mspline2.tmp";
  cout << sys << endl;
  system(sys.c_str());

  char c;

  ifstream ifs2(".\\doscmd.mspline2.tmp");
  while (ifs2.get(c))
    sout.put(c);
  ifs2.close();

  ifstream ifs1(".\\doscmd.mspline.tmp");
  while (ifs1.get(c))
    sout.put(c);
  ifs1.close();

  system("erase .\\doscmd.mspline.tmp");
  system("erase .\\doscmd.mspline2.tmp");
}

void listRegFolder(string folder,ostream &sout) {
  // List the registry settings
  int n = settings->getKeyArraySize();
  int i;

  if (n < 2)
    {
    sout << "Empty Folder" << endl;
    }
  else
    {
    char **keys = new char*[n];
    settings->getFolderKeys(keys);
    for (i=0;keys[i];i++)
      {
      sout << keys[i] << ": " << "SubFolder" << endl;
      }       
    settings->getValueKeys(keys);
    for (i=0;keys[i];i++)
      {
      sout << keys[i] << ": " << settings->getStringValue(keys[i],"undefined") << endl;
      }
    delete keys;
    }   
}

void exportMatlab() {


}

class TrianglePatch
  {
public:
  vector<SMLVec3f> X,N;
  vector<int> T;
  };

class SplineTriangleData
  {
public:
  Array2D<TrianglePatch> M,B[2];
  };

void exportTriangleData(SplineTriangleData &T) {
  T.M.resize(splineData->patch.width(),splineData->patch.height());
  T.B[0].resize(splineData->patch.width(),splineData->patch.height());
  T.B[1].resize(splineData->patch.width(),splineData->patch.height());

  for (int i=0;i<splineData->patch.width();i++)
    {
    for (int j=0;j<splineData->patch.height();j++)
      {

      PatchDataCache *pdc = splineData->patch(i,j);
      TrianglePatch &TM = T.M(i,j);
      TrianglePatch &TB0 = T.B[0](i,j);
      TrianglePatch &TB1 = T.B[1](i,j);


      // Number of triangles
      int nTri = (pdc->MP.width()-1) * (pdc->MP.height()-1) * 2 + pdc->idxTrimStripIdx / 3;
      int nVtx = (pdc->MP.width() * pdc->MP.height()) + pdc->getCrestSize();

      // Index mapping
      Aux2D<int> idxMap(pdc->MP.width(),pdc->MP.height(),pdc->MP.auxSize());

      // Save the vertices
      for (int z=0;z<nVtx;z++)
        {
        if (pdc->in(z) & 0x01)
          {
          idxMap(z) = TM.X.size();
          TM.X.push_back(pdc->MP(z).X());
          TM.N.push_back(pdc->MP(z).N());
          TB0.X.push_back(pdc->MP(z).bp[0].X);
          TB0.N.push_back(pdc->MP(z).bp[0].N);
          TB1.X.push_back(pdc->MP(z).bp[1].X);
          TB1.N.push_back(pdc->MP(z).bp[1].N);                    
          }
        else
          {
          idxMap(z) = -1;
          }
        }

      for (int v=0;v<pdc->MP.height()-1;v++)
        {
        for (int u=0;u<pdc->MP.width()-1;u++)
          {
          if (pdc->in(u,v) & pdc->in(u+1,v) & pdc->in(u,v+1) & pdc->in(u+1,v+1) & 0x01)
            {
            TM.T.push_back(idxMap(u,v));
            TM.T.push_back(idxMap(u,v+1));
            TM.T.push_back(idxMap(u+1,v+1));
            TM.T.push_back(idxMap(u,v));
            TM.T.push_back(idxMap(u+1,v+1));
            TM.T.push_back(idxMap(u+1,v));

            TB0.T.push_back(idxMap(u,v));
            TB0.T.push_back(idxMap(u,v+1));
            TB0.T.push_back(idxMap(u+1,v+1));
            TB0.T.push_back(idxMap(u,v));
            TB0.T.push_back(idxMap(u+1,v+1));
            TB0.T.push_back(idxMap(u+1,v));

            TB1.T.push_back(idxMap(u,v));
            TB1.T.push_back(idxMap(u+1,v+1));
            TB1.T.push_back(idxMap(u,v+1));
            TB1.T.push_back(idxMap(u,v));
            TB1.T.push_back(idxMap(u+1,v));
            TB1.T.push_back(idxMap(u+1,v+1));
            }
          }
        }
      for (int w=0;w<pdc->idxTrimStripIdx;w+=3)
        {
        TM.T.push_back(idxMap(pdc->idxTrimStrip[w]));
        TM.T.push_back(idxMap(pdc->idxTrimStrip[w+2]));
        TM.T.push_back(idxMap(pdc->idxTrimStrip[w+1]));
        TB0.T.push_back(idxMap(pdc->idxTrimStrip[w]));
        TB0.T.push_back(idxMap(pdc->idxTrimStrip[w+1]));
        TB0.T.push_back(idxMap(pdc->idxTrimStrip[w+2]));
        TB1.T.push_back(idxMap(pdc->idxTrimStrip[w]));
        TB1.T.push_back(idxMap(pdc->idxTrimStrip[w+2]));
        TB1.T.push_back(idxMap(pdc->idxTrimStrip[w+1]));
        }
      }
    }
}

void fprintPovRay(FILE *f,TrianglePatch &TP) {
  int z;

  if (TP.T.size() > 0)
    {
    fprintf(f,"mesh2 {\n");
    fprintf(f,"   vertex_vectors {\n");
    fprintf(f,"      %d,\n",TP.X.size());
    for (z=0;z<TP.X.size();z++)
      {
      SMLVec3f X = TP.X[z];
      fprintf(f,"      <%lg,%lg,%lg>%s\n",X[0],X[1],X[2],(z==TP.X.size()-1 ? "" : ","));
      }
    fprintf(f,"   }\n");

    fprintf(f,"   normal_vectors {\n");
    fprintf(f,"      %d,\n",TP.N.size());
    for (z=0;z<TP.N.size();z++)
      {
      SMLVec3f N = TP.N[z];
      fprintf(f,"      <%lg,%lg,%lg>%s\n",N[0],N[1],N[2],(z==TP.N.size()-1 ? "" : ","));
      }
    fprintf(f,"   }\n");

    fprintf(f,"   face_indices {\n");
    fprintf(f,"      %d,\n",TP.T.size()/3);
    int idx = 0;
    for (z=0;z<TP.T.size();z+=3)
      {
      fprintf(f,"      <%d,%d,%d>%s\n",TP.T[z],TP.T[z+1],TP.T[z+2],(z == TP.T.size()-1) ? "" : ",");
      }
    fprintf(f,"   }\n");
    fprintf(f,"}\n\n");
    }
  else
    {
    fprintf(f,"intersection {\n   sphere {<1,1,1> 0.1}\n   sphere {<0,0,0> 0.1}\n}\n");
    }
}

void exportPovRay(const char *file) {
  SplineTriangleData T;
  exportTriangleData(T);

  FILE *f = fopen(file,"wt");

  if (f == NULL)
    throw "Can not open file for writing";

  // Print out patch arrays
  fprintf(f,"#declare patchBnd = array[%d][%d][2]\n",splineData->patch.width(),splineData->patch.height());
  fprintf(f,"#declare patchMed = array[%d][%d]\n\n",splineData->patch.width(),splineData->patch.height());

  // Print each patch
  for (int i=0;i<splineData->patch.width();i++)
    {
    for (int j=0;j<splineData->patch.height();j++)
      {
      // Medial patch
      fprintf(f,"#declare patchMed[%d][%d] = \n",i,j);
      fprintPovRay(f,T.M(i,j));

      // Boundary patches
      for (int d=0;d<2;d++)
        {
        fprintf(f,"#declare patchBnd[%d][%d][%d] = \n",i,j,d);
        fprintPovRay(f,T.B[d](i,j));
        }
      }
    }

  fclose(f);
}

int fprintWaveFront(FILE *f,TrianglePatch &TP,int vFirst=1) {
  int z=0;
  if (TP.T.size() > 0)
    {
    for (z=0;z<TP.X.size();z++)
      {
      SMLVec3f X = TP.X[z];
      SMLVec3f N = TP.N[z];
      fprintf(f,"v %lg %lg %lg\n",X[0],X[1],X[2]);
      fprintf(f,"vn %lg %lg %lg\n",N[0],N[1],N[2]);
      }       

    for (z=0;z<TP.T.size();z+=3)
      {
      int v[3];
      v[0] = vFirst + TP.T[z];
      v[1] = vFirst + TP.T[z+1];
      v[2] = vFirst + TP.T[z+2];
      fprintf(f,"f %d//%d %d//%d %d//%d\n",v[0],v[0],v[1],v[1],v[2],v[2]);
      }
    }

  return vFirst + TP.T.size();
}



/**
 * Export to an OBJ file
 */
void exportWaveFront(const char *file) {
  SplineTriangleData T;
  exportTriangleData(T);

  FILE *f = fopen(file,"wt");

  // Print each patch
  int idx = 1;
  for (int i=0;i<splineData->patch.width();i++)
    {
    for (int j=0;j<splineData->patch.height();j++)
      {
      // Medial patch
      fprintf(f,"g patchMed%02d%02d\n",i,j);
      idx = fprintWaveFront(f,T.M(i,j),idx);

      // Boundary patches
      for (int d=0;d<2;d++)
        {
        fprintf(f,"g patchBnd%02d%02d_%d\n",i,j,d);
        idx = fprintWaveFront(f,T.B[d](i,j),idx);
        }
      }
    }
}


struct CPEffect
  {
  // Derivative of Y w.r.t each parameter of each control point
  // First indexed by control point, last index x,y,z,w;
  SMLVec3f dYdP[2][4][4][4];

  // Precomputed terms for the derivatives

  // Image gradient
  SMLVec3f G[2];
  };

class FlowTester : public GLDisplayListRenderer, public BoundaryMeasure
  {
private:
  Array2D<SMLVec4f> grad,spl;
  bool on;

  // This huge structure represents all points
  Aux2D< Aux2D <CPEffect> *> DY;

  // Current control index (0..3)
  int iCtl,jCtl,iDim;

  // The current image
  AbstractImage3D *img;

public:
  FlowTester() {
    on = false;
  }

  void build();
  void update(ostream &);
  void apply();
  void computeControlGradient(ostream &out);

  float computeBoundaryMeasure(const MedialPoint &mp,int side);
  float computeCrestBoundaryMeasure(const MedialPoint &mp);
  };
FlowTester *flowTester;


void FlowTester::build() {
  if (on)
    {
    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glColor3d(1,0,0);

    glPushMatrix();
    glTranslated(-0.5,-0.5,-0.5);

    glBegin(GL_LINES);

    for (int i=0;i<spl.width();i++)
      {
      for (int j=0;j<spl.height();j++)
        {
        glVertex3d(spl(i,j)[0],spl(i,j)[1],spl(i,j)[2]);
        glVertex3d(spl(i,j)[0]+grad(i,j)[0],spl(i,j)[1]+grad(i,j)[1],spl(i,j)[2]+grad(i,j)[2]);
        }
      }

    glEnd();

    glPopMatrix();

    glPopAttrib();
    }
}



inline SMLVec3f operator *(float f,SMLVec3f vec) {
  return vec * f;
}

// Compute for each point in the spline 
void FlowTester::computeControlGradient(ostream &out) {
  double tStart = clock();

  if ((img = imaging.getImage())== NULL)
    {
    on = false;
    return;
    }

  int m = spline->dim(0);
  int n = spline->dim(1);
  DY.resize(m-2,n-2,0);

  splineData->refreshMedialAllPatches(splineData->getMaxLevel());
  splineData->refreshBoundaryAllPatches(splineData->getMaxLevel());

  for (int jPatch=0;jPatch<splineData->patch.height();jPatch++)
    {
    for (int iPatch=0;iPatch<splineData->patch.width();iPatch++)
      {

      PatchDataCache *pdc = splineData->patch(iPatch,jPatch);
      DY(iPatch,jPatch) = new Aux2D<CPEffect>(pdc->MP.width(),pdc->MP.height(),pdc->getCrestSize());
      Aux2D<CPEffect> &DYP = *DY(iPatch,jPatch);

      for (int k=0;k<DYP.size();k++)
        {
        // Ignore the points outside 
        if (!(pdc->in(k) & 0x01))
          continue;

        CPEffect &cpe = DYP(k);
        MedialPoint &mp = pdc->MP(k);

        // Set the data pointer
        mp.data = &cpe;

        // Create the matrix of alpha, beta, gamma for surrounding control points
        float alpha[4][4],beta[4][4],gamma[4][4];
        SMLVec4f Wu[4],Wv[4];

        // This is unnecessary and should be stored with the off-grid points
        spline->basisJet(0,iPatch+3,mp.u,Wu);
        spline->basisJet(1,jPatch+3,mp.v,Wv);

        // Compute the contribution of each control point and its derivative
        for (int s=0;s<4;s++)
          {
          for (int t=0;t<4;t++)
            {
            alpha[s][t] = Wu[0][s] * Wv[0][t];
            beta[s][t] = Wu[1][s] * Wv[0][t];
            gamma[s][t] = Wu[0][s] * Wv[1][t];
            }
          }

        // Access parts of the medial point
        SMLVec3f &X   = mp.X();
        SMLVec3f &Xu  = mp.Xu();
        SMLVec3f &Xv  = mp.Xv();
        SMLVec3f &N   = mp.N();
        float &R   = mp.R();
        float &Ru  = mp.Ru();
        float &Rv  = mp.Rv();
        float &E = mp.IE, &F = mp.IF, &G = mp.IG;
        float EGF2 = (E*G-F*F);
        float _EGF2 = 1.0 / EGF2;
        float _nrmN = sqrt(_EGF2);

        SMLVec3f B = -mp.GradR;
        float nrmB2 = B.squared_magnitude();
        float nrmB = sqrt(nrmB2);
        float _nrmB = 1.0 / nrmB;
        float sinTh = sqrt(1 - nrmB2);
        float _sinTh = 1.0 / sinTh;

        // Compute the components of the Bp vector 
        for (int p=0;p<3;p++)
          {
          // Corresponding unit vector 
          SMLVec3f Ek(0,0,0);
          Ek[p] = 1.0f;

          // Some precomputed terms
          float EXU = Xu[p];
          float EXV = Xv[p];
          float EXURv = EXU*Rv;
          float EXVRu = EXV*Ru;

          // Coefficients of beta and gamma in Bp
          SMLVec3f CBPb  = - _EGF2 * (
                                     (2 * EXURv - EXVRu) * Xv
                                     - (EXV * Rv) * Xu 
                                     + (G * Ru - F * Rv) * Ek
                                     - (2 * (EXV * F - EXU * G)) * B);

          SMLVec3f CBPg  = - _EGF2 * (
                                     (2 * EXVRu - EXURv) * Xu
                                     - (EXU * Ru) * Xv 
                                     + (E * Rv - F * Ru) * Ek
                                     - (2 * (EXU * F - EXV * E)) * B);

          // Coefficients of beta and gamma in Np
          SMLVec3f XuxEk = vnl_cross_3d(Xu,Ek);
          SMLVec3f EkxXv = vnl_cross_3d(Ek,Xv);
          SMLVec3f CNPg  = _nrmN * (XuxEk - N * dot_product(N,XuxEk));
          SMLVec3f CNPb  = _nrmN * (EkxXv - N * dot_product(N,EkxXv));

          for (int d=0;d<2;d++)
            {
            BoundaryPoint &bp = mp.bp[d];
            float D = 2 * d - 1;

            // Coefficients for beta and gamma in the result
            SMLVec3f Cg = R * (CBPg + D * _sinTh * ((1-nrmB2) * CNPg - dot_product(B,CBPg)*N));
            SMLVec3f Cb = R * (CBPb + D * _sinTh * ((1-nrmB2) * CNPb - dot_product(B,CBPb)*N));
            SMLVec3f Ca = Ek;

            // Store the effects
            for (int s=0;s<4;s++)
              {
              for (int t=0;t<4;t++)
                {
                cpe.dYdP[d][s][t][p] = 
                Ca * alpha[s][t] + 
                Cb * beta[s][t] + 
                Cg * gamma[s][t];
                }
              }

            SMLVec3f T = cpe.dYdP[d][0][0][p];
            if (T.magnitude() > 4)
              {
              cout << "error " << endl;
              }
            }
          }

        // Now do this for the radial controls
        SMLVec3f CBPg  = _EGF2 * (F*Xu - E*Xv);
        SMLVec3f CBPb  = _EGF2 * (F*Xv - G*Xu);

        for (int d=0;d<2;d++)
          {
          BoundaryPoint &bp = mp.bp[d];
          float D = 2 * d - 1;

          // Coefficients for beta and gamma in the result
          SMLVec3f Cg = R * (CBPg - D * _sinTh * (dot_product(B,CBPg)*N));
          SMLVec3f Cb = R * (CBPb - D * _sinTh * (dot_product(B,CBPb)*N));
          SMLVec3f Ca = B + (D * sinTh) * N;

          // Store the effects
          for (int s=0;s<4;s++)
            {
            for (int t=0;t<4;t++)
              {
              cpe.dYdP[d][s][t][3] = 
              Ca * alpha[s][t] + 
              Cb * beta[s][t] + 
              Cg * gamma[s][t];
              }
            }

          // Store the image gradient
          img->interpolateVoxelGradient(bp.X[0],bp.X[1],bp.X[2],cpe.G[d].data_block());
          }
        }
      }
    }

  // Initialize the gradient arrays
  grad.resize(m+1,n+1);
  spl.resize(m+1,n+1);

  // Go through all the control points 
  for (int jControl=0;jControl<=n;jControl++)
    {
    int jPatchStart = vnl_math_max(jControl - 3,0);
    int jPatchEnd = vnl_math_min(jControl,n-3);

    for (int iControl=0;iControl<=m;iControl++)
      {
      int iPatchStart = vnl_math_max(iControl - 3,0);
      int iPatchEnd = vnl_math_min(iControl,m-3);

      for (iDim=0;iDim<4;iDim++)
        {
        // Get the area of the measure
        float area = 0, totArea = 0;
        float measure = 0;

        // Go through all affected patches and integrate our measurement
        jCtl = jControl-jPatchStart;
        for (int jPatch=jPatchStart;jPatch<=jPatchEnd;jPatch++)
          {
          iCtl = iControl-iPatchStart;
          for (int iPatch=iPatchStart;iPatch<=iPatchEnd;iPatch++)
            {
            PatchDataCache *pdc = splineData->patch(iPatch,jPatch);
            measure += pdc->integrateBoundaryMeasure(*this,area);
            totArea += area;
            iCtl--;
            }
          jCtl--;
          }

        // Set the change
        measure *= 0.0001;
        grad(iControl,jControl)[iDim] = totArea > 0.0f ? -measure / totArea : 0.0f;
        }

      spl(iControl,jControl) = spline->getControl4f(iControl,jControl);
      }
    }

  out << "Time elapsed: " << (clock()-tStart) << " ms." << endl;
  on = true;
  invalidate();
}

float FlowTester::computeBoundaryMeasure(const MedialPoint &mp,int side) {
  // Get the effect thingy
  CPEffect *cpe = (CPEffect*)mp.data;

  // Get the displacement vector 
  SMLVec3f dYdp = cpe->dYdP[side][iCtl][jCtl][iDim];
  return dot_product(dYdp,cpe->G[side]);


  // SMLVec3f X = mp.bp[side].X;
  // SMLVec3f G;
  // img->interpolateVoxelGradient(X.x,X.y,X.z,G.data());

  // Return the point change in volume
  //SMLVec3f N = mp.bp[side].N;
  //return dYdp.Dot(N);

  // float dIdp = dYdp.Dot(G);
  // return dIdp;
}

float FlowTester::computeCrestBoundaryMeasure(const MedialPoint &mp) {
  return 0;
}


void FlowTester::update(ostream &out) {
  int i,j;

  SplineImageMatcher *match = imaging.getImageMatch();
  MRepConstraint mrc(spline);

  static float factor = 1.0;

  double tStart = clock();

  if (match)
    {
    match->setSpline(spline,splineData);

    int w = spline->dim(0)+1;
    int h = spline->dim(1)+1;
    grad.resize(w,h);
    spl.resize(w,h);

    // Max grad length
    float maxLen = 0;
    float maxR = 0;

    cout << "=================================================================";

    // Save the current function value
    float f = match->getMatch();// + mrc.priorProbabilityFull(spline,splineData);


    for (i=0;i<w;i++)
      {
      for (j=0;j<h;j++)
        {
        spl(i,j) = spline->getControl4f(i,j);
        for (int k=0;k<4;k++)
          {
          SMLVec4f cp = spl(i,j);
          cp[k] += 0.001;
          spline->setControl(i,j,cp);

          float f1 = match->getMatch();// + mrc.priorProbabilityFull(spline,splineData);
          float df = (f1 - f) / 0.001;
          grad(i,j)[k] = df;
          }
        spline->setControl(i,j,spl(i,j));

        float len = sqrt(grad(i,j)[0]*grad(i,j)[0] + grad(i,j)[1]*grad(i,j)[1] + grad(i,j)[2]*grad(i,j)[2]);

        cout << i << "\t" << j << "\t" << grad(i,j)[0] << "\t" << grad(i,j)[1] 
        << "\t" << grad(i,j)[2] << "\t" << grad(i,j)[3] << endl;

        //grad(i,j).x /= len;
        //grad(i,j).y /= len;
        //grad(i,j).z /= len;
        //grad(i,j).w = 1.0;

        //maxLen = maxLen > len ? maxLen : len;
        //maxR = maxR > fabs(grad(i,j).w) ? maxR : fabs(grad(i,j).w);

        }
      }

    factor *= 0.95;

    for (i=0;i<w;i++)
      {
      for (j=0;j<h;j++)
        {
        // grad(i,j).x *= - 0.05 / maxLen;
        // grad(i,j).y *= - 0.05 / maxLen;
        // grad(i,j).z *= - 0.05 / maxLen;
        // grad(i,j).w *= - 0.05 / maxR;
        grad(i,j) *= -0.0001 * factor;              
        }
      }
    on = true;      
    }
  else
    {
    on = false;
    }
  invalidate();

  out << "Time elapsed: " << (clock()-tStart) << " ms." << endl;
}

void FlowTester::apply() {
  for (int i=0;i<spl.width();i++)
    {
    for (int j=0;j<spl.height();j++)
      {
      spl(i,j) += grad(i,j);
      grad(i,j).fill(0);
      spline->setControl(i,j,spl(i,j));
      }
    }
  invalidate();
}

/**
 * Check that my image gradient code is valid
 */
void testImageDerivative() {
  AbstractImage3D *img = imaging.getImage();
  if (!img)
    return;

  // Constants
  const SMLVec3f X0(0.53,0.47,0.44);
  const float eps = 0.001;

  // Compute analytical derivative    
  SMLVec3f G2;
  for (int k=0;k<3;k++)
    {
    SMLVec3f X1 = X0;
    X1[k] += eps;
    G2[k] = (img->interpolateVoxel(X1[0],X1[1],X1[2]) - img->interpolateVoxel(X0[0],X0[1],X0[2])) / eps;
    }

  // Compute the actual derivative
  SMLVec3f G;
  img->interpolateVoxelGradient(X0[0],X0[1],X0[2],G.data_block());

  // Compute the difference
  SMLVec3f delta = G - G2;
  cout << "Difference is " << delta.magnitude() << endl;
}

/**
 * Test that the derivative computation code I use is valid
 */
void testDerivatives() {

  const int ip0 = 2;
  const int jp0 = 2;
  const int ic0 = 2;
  const int jc0 = 4;
  const int u0 = 13;
  const int v0 = 11;
  const float eps = 0.01;
  const int k0 = 3;
  const int d0 = 1;
/*
    const int ip0 = 0;
    const int jp0 = 0;
    const int ic0 = 0;
    const int jc0 = 0;
    const int u0 = 13;
    const int v0 = 8;
    const float eps = 0.01;
    const int k0 = 0;
    const int d0 = 0;
*/

  // Make sure the data is all computed
  splineData->refreshMedialAllPatches(splineData->getMaxLevel());
  splineData->refreshBoundaryAllPatches(splineData->getMaxLevel());

  // Make a copy of the spline
  MSpline s2(spline->dim(0),spline->dim(1));
  for (int jControl=0;jControl<=spline->dim(1);jControl++)
    {
    for (int iControl=0;iControl<=spline->dim(0);iControl++)
      {
      s2.setControl(iControl,jControl,spline->getControl4f(iControl,jControl));
      }
    }
  SplineDataCache sdc2(&s2,splineData->getMaxLevel());

  // Update the spline
  s2.setControl(ic0,jc0,k0,s2.getControl(ic0,jc0,k0)+eps);
  sdc2.refreshMedialAllPatches(sdc2.getMaxLevel());
  sdc2.refreshBoundaryAllPatches(sdc2.getMaxLevel());

  // Compute the difference vector (analytical derivative)
  PatchDataCache *pdc2 = sdc2.patch(ip0,jp0);
  PatchDataCache *pdc = splineData->patch(ip0,jp0);

  SMLVec3f dydp2 = (pdc2->MP(u0,v0).bp[d0].X - pdc->MP(u0,v0).bp[d0].X) / eps;
  SMLVec3f dB2 = (pdc2->MP(u0,v0).GradR - pdc->MP(u0,v0).GradR) / eps;

  // Now compute the exact derivative using my method
  MedialPoint &mp = pdc->MP(u0,v0);
  BoundaryPoint &bp = mp.bp[d0];

  // Get the jet
  SMLVec4f Wu[4],Wv[4];
  spline->basisJet(0,ip0+3,mp.u,Wu);
  spline->basisJet(1,jp0+3,mp.v,Wv);

  // Compute the contribution of each control point and its derivative
  float alpha = Wu[0][ic0-ip0] * Wv[0][jc0-jp0];
  float beta = Wu[1][ic0-ip0] * Wv[0][jc0-jp0];
  float gamma = Wu[0][ic0-ip0] * Wv[1][jc0-jp0];

  // Check alpha, beta, gamma
  float alpha2 = (pdc2->MP(u0,v0).X() - pdc->MP(u0,v0).X())[k0] / eps;
  float beta2 = (pdc2->MP(u0,v0).Xu() - pdc->MP(u0,v0).Xu())[k0] / eps;
  float gamma2 = (pdc2->MP(u0,v0).Xv() - pdc->MP(u0,v0).Xv())[k0] / eps;

  // Access parts of the medial point
  SMLVec3f &X   = mp.X();
  SMLVec3f &Xu  = mp.Xu();
  SMLVec3f &Xv  = mp.Xv();
  SMLVec3f &N   = mp.N();
  float &R   = mp.R();
  float &Ru  = mp.Ru();
  float &Rv  = mp.Rv();
  float &E = mp.IE, &F = mp.IF, &G = mp.IG;
  float EGF2 = (E*G-F*F);
  float _EGF2 = 1.0 / EGF2;
  float _nrmN = sqrt(_EGF2);

  SMLVec3f B =  - _EGF2 * ((E*Xv-F*Xu)*Rv + (G*Xu-F*Xv)*Ru);
  float nrmB2 = B.squared_magnitude();
  float nrmB = sqrt(nrmB2);
  float _nrmB = 1.0 / nrmB;
  float sinTh = sqrt(1 - nrmB2);
  float _sinTh = 1.0 / sinTh;
  float D = 2 * d0 - 1;

  SMLVec3f dydp;

  // Compute the components of the Bp vector 
  if (k0 < 3)
    {
    // Corresponding unit vector 
    SMLVec3f Ek(0,0,0);
    Ek[k0] = 1.0f;

    // Some precomputed terms
    float EXU = Xu[k0];
    float EXV = Xv[k0];

    // Coefficients of beta and gamma in Bp
    SMLVec3f CBPb  = - _EGF2 * (
                               (2 * EXU * Rv - EXV * Ru) * Xv
                               - (EXV * Rv) * Xu 
                               + (G * Ru - F * Rv) * Ek
                               + (2 * (EXV * F - EXU * G)) * mp.GradR);

    SMLVec3f CBPg  = - _EGF2 * (
                               (2 * EXV * Ru - EXU * Rv) * Xu
                               - (EXU * Ru) * Xv 
                               + (E * Rv - F * Ru) * Ek
                               + (2 * (EXU * F - EXV * E)) * mp.GradR);

    //SMLVec3f CBPgg = _EGF2 * (-E*B - Xu*Ru);
    //SMLVec3f CBPbg = _EGF2 * (Xu*Rv + Xv*Ru + 2*F*B);
    //SMLVec3f CBPbb = _EGF2 * (-G*B - Xv*Rv);

    SMLVec3f dB = CBPb * beta + CBPg * gamma;

    // Coefficients of beta and gamma in Np
    SMLVec3f XuxEk = vnl_cross_3d(Xu,Ek);
    SMLVec3f EkxXv = vnl_cross_3d(Ek,Xv);
    SMLVec3f CNPg  = _nrmN * (XuxEk - N * dot_product(N,XuxEk));
    SMLVec3f CNPb  = _nrmN * (EkxXv - N * dot_product(N,EkxXv));

    // Coefficients for beta and gamma in the result
    SMLVec3f Cg = R * (CBPg + D * _sinTh * ((1-nrmB2) * CNPg - dot_product(B,CBPg)*N));
    SMLVec3f Cb = R * (CBPb + D * _sinTh * ((1-nrmB2) * CNPb - dot_product(B,CBPb)*N));
    SMLVec3f Ca = Ek;

    dydp = Ca * alpha + Cb * beta + Cg * gamma;
    }
  else
    {
    // Now do this for the radial controls
    SMLVec3f CBPg  = _EGF2 * (F*Xu - E*Xv);
    SMLVec3f CBPb  = _EGF2 * (F*Xv - G*Xu);

    // Coefficients for beta and gamma in the result
    SMLVec3f Cg = R * (CBPg - D * _sinTh * (dot_product(B,CBPg)*N));
    SMLVec3f Cb = R * (CBPb - D * _sinTh * (dot_product(B,CBPb)*N));
    SMLVec3f Ca = B + (D * sinTh) * N;

    // Store the effects
    dydp = Ca * alpha + Cb * beta + Cg * gamma;
    }

  // Compare the two vectors
  SMLVec3f delta = dydp - dydp2;
  cout << "Difference is " << delta.magnitude() << endl;

  // Compute difference in image match at the two positions
  AbstractImage3D *img = imaging.getImage();
  if (img)
    {
    SMLVec3f X1 = pdc2->MP(u0,v0).bp[d0].X;
    SMLVec3f X0 = pdc->MP(u0,v0).bp[d0].X;

    float di2 = (img->interpolateVoxel(X1[0],X1[1],X1[2]) - img->interpolateVoxel(X0[0],X0[1],X0[2])) / eps;

    SMLVec3f GI;
    img->interpolateVoxelGradient(X0[0],X0[1],X0[2],GI.data_block());

    float di = dot_product(GI,dydp);
    cout << di << endl;
    }
}



bool matlabMode = false;
/*
void processMatlabCommand(string cmd, ostream &sout) {
    if(cmd == "quit") {
        sout << "Exiting MATLAB mode" << endl;
        matlabMode = false;
    }
    else if(cmd == "exp" || cmd=="export") {
        for(int i=0;i<splineData->patch.width();i++) {
            for(int j=0;j<splineData->patch.height();j++) {
                char name[12];
                PatchDataCache *pdc = splineData->patch(i,j);
                
                sprintf(name,"x%02d%02d",i,j);
                matlab->exportMatrix(name,&(pdc->MP(0,0).F.x),pdc->MP.height(),pdc->MP.width(),sizeof(MedialPoint));
                sprintf(name,"y%02d%02d",i,j);
                matlab->exportMatrix(name,&(pdc->MP(0,0).F.y),pdc->MP.height(),pdc->MP.width(),sizeof(MedialPoint));
                sprintf(name,"z%02d%02d",i,j);
                matlab->exportMatrix(name,&(pdc->MP(0,0).F.z),pdc->MP.height(),pdc->MP.width(),sizeof(MedialPoint));
            }
        }
            
            
    }
    else {
        matlab->command(cmd.c_str(),sout);
    }
}
*/
class CommandString
  {
public:
  CommandString(const char *chrCommand) {
    command = chrCommand;
    istringstream iss(command);
    while (!iss.eof())
      {
      string arg;
      iss >> arg;
      args.push_back(arg);
      }
    cursor = 0;
  }

  CommandString &operator >> (string &s) {
    s = (cursor < args.size()) ? args[cursor++] : "";
    return *this;
  }

  CommandString &operator >> (int &i) {
    string s;
    (*this) >> s;

    if (s.length())
      i = atoi(s.c_str());

    return *this;
  }

  CommandString &operator >> (double &d) {
    string s;
    (*this) >> s;

    if (s.length())
      d = atof(s.c_str());

    return *this;
  }

  CommandString &operator >> (float &d) {
    string s;
    (*this) >> s;

    if (s.length())
      d = atof(s.c_str());

    return *this;
  }

  // Reutn the number of remaining args
  int left() {
    return args.size() - cursor;
  }

  // Reset the arguments
  void reset() {
    cursor = 0;
  }

  // Get the whole command
  string str() {
    return command;
  }

private:
  vector<string> args;
  int cursor;
  string command;
  };

void processCommand(const char *chrCommand, ostream &sout) {
  try
    {
    //string src(chrCommand);
    //istringstream iss(src);
    CommandString iss(chrCommand);

    // Get the command start
    string cmd;
    iss >> cmd;

    // Are we in matlab mode?
    if (matlabMode)
      {
      sout << "MATLAB support disabled in this version!" << endl;
      // processMatlabCommand(cmd,sout);
      }

    // What is the command
    else if (cmd == "help")
      {
      // DOS command
      sout << "No help is available yet.  Please send money to pauly@cs.unc.edu" << endl;
      }

    else if (cmd == "load")
      {
      iss >> cmd; 

      if (cmd == "bin")
        {
        iss >> cmd;
        if (loadImage(cmd.c_str(),ImagingSystem::BINARY,sout))
          {
          sout << "Loaded image " << cmd << endl;
          }
        }

      else if (cmd == "gray" || cmd=="grey")
        {
        iss >> cmd;
        if (loadImage(cmd.c_str(),ImagingSystem::GRAYSCALE,sout))
          {
          sout << "Loaded image " << cmd << endl;
          }
        }

      else if (cmd == "distance" || cmd == "dt")
        {
        iss >> cmd;
        if (loadImage(cmd.c_str(),ImagingSystem::DISTANCE,sout))
          {
          sout << "Loaded distance transform " << cmd << endl;
          }
        }

      else if (cmd == "spline" || cmd == "sp")
        {
        iss >> cmd;
        if (loadSpline(cmd.c_str(),sout))
          {
          sout << "Loaded spline " << cmd << endl;
          }
        }

      else if (cmd == "edge")
        {
        iss >> cmd;
        SMLVec3f voxSize(1.0f,1.0f,1.0f);

        if (iss.left() >= 3)
          {
          iss >> voxSize[0];
          iss >> voxSize[1];
          iss >> voxSize[2];
          }

        edgeNetRenderer.loadNetData(cmd.c_str(),voxSize);
        edgeNetRenderer.setVisible(true);
        }
      }

    else if (cmd == "cd")
      {
      char buff[256];
      sout << iss.str().substr(3).c_str() << endl;
      sout << CHDIR(iss.str().substr(6).c_str());
      sout << " : Current directory" << endl << GETCWD(buff,256) << endl;
      }

    else if (cmd == "save")
      {
      iss >> cmd; 

      if (cmd == "image" || cmd == "im")
        {
        iss >> cmd;
        if (saveImage(cmd.c_str(),sout))
          {
          sout << "Saved image " << cmd << endl;
          }
        }

      else if (cmd == "spline" || cmd == "sp")
        {
        iss >> cmd;
        if (saveSpline(cmd.c_str(),sout))
          {
          sout << "Saved spline " << cmd << endl;
          }
        }

      else if (cmd == "pov" || cmd == "povray")
        {
        iss >> cmd;
        exportPovRay(cmd.c_str());
        sout << "Exported POV-Ray File" << cmd << endl;             
        }

      else if (cmd == "obj" || cmd == "wave")
        {
        iss >> cmd;
        exportWaveFront(cmd.c_str());
        sout << "Exported Wavefront File" << cmd << endl;               
        }
      }

    else if (cmd == "distanceMap" || cmd == "distmap")
      {
      try
        {
        onImageChange();
        imaging.toDistance();
        sout << "Image conveted to distance transform." << endl;
        }
      catch (string exc)
        {
        sout << exc << endl;
        }
      }

    else if (cmd == "convert")
      {
      iss >> cmd;
      if (cmd == "image")
        {
        iss >> cmd;
        if (cmd == "binary" || cmd == "bin")
          {
          imaging.toBinary();
          onImageChange();
          }
        else if(cmd == "grey" || cmd == "gray")
          {
          imaging.toGrayscale();
          onImageChange();
          }
        }
      }

    else if (cmd == "mask")
      {
      if (imaging.getType() == imaging.DISTANCE)
        {
        DistanceTransform *dt = (DistanceTransform*)imaging.getImage();
        onImageChange();
        iss >> cmd;

        float amtAdd = 1000.0f;
        float amtScale = 0.0f;
        iss >> amtAdd;
        iss >> amtScale;

        if (cmd == "inside" || cmd == "in")
          {
          dt->maskInside(amtAdd,amtScale);
          }
        else if (cmd == "outside" || cmd == "out")
          {
          dt->maskOutside(amtAdd,amtScale);
          }
        else
          {
          sout << "specify 'in' or 'out'." << endl;
          }
        }
      else
        {
        sout << "Image not a distance transform." << endl;
        }                       
      }

    else if (cmd == "blur")
      {
      float scale;
      iss >> scale;


      if (imaging.getType() == imaging.GRAYSCALE)
        {
        SMLVec3f sigma = imaging.getImage()->getVoxelSize() * scale;
        imaging.getImage()->gaussianBlur(sigma.data_block());    
        sout << "Image has been blurred" << endl;
        }
      else
        {
        sout << "Operation undefined for this image type" << endl;
        }
      }

    else if (cmd == "list" || cmd == "ls")
      {
      iss >> cmd; 
      if (cmd == "images" || cmd == "i")
        {
        ostringstream oss;
        oss << "dir ";
        // oss << makeImageFileName("*.raw*");
        oss << "*.raw*";
        callShellCommand(oss.str().c_str(),sout);
        }
      else if (cmd == "s" || cmd == "spl" || cmd == "splines" || cmd == "spline")
        {
        callShellCommand("dir /B *.spl*",sout);
        }
      }

    else if (cmd == "show")
      {
      iss >> cmd; 
      if (cmd == "slice")
        {
        iss >> cmd; 
        if (cmd == "off")
          {
          sliceRenderer.setVisible(false);
          }
        else if (cmd == "pause")
          {
          sliceRenderer.sliceMoving = !sliceRenderer.sliceMoving;
          char axis = 'X' + sliceRenderer.dim;
          sout << axis << " axis slice " << sliceRenderer.sliceNo << endl;
          }
        else if (cmd == "last")
          {
          sliceRenderer.sliceMoving = false;
          sliceRenderer.lastSlice();
          char axis = 'X' + sliceRenderer.dim;
          sout << axis << " axis slice " << sliceRenderer.sliceNo << endl;
          }
        else if (cmd == "next")
          {
          sliceRenderer.sliceMoving = false;
          sliceRenderer.nextSlice();
          char axis = 'X' + sliceRenderer.dim;
          sout << axis << " axis slice " << sliceRenderer.sliceNo << endl;
          }
        else if (cmd == "axis")
          {
          sliceRenderer.sliceNo = 0;
          sliceRenderer.dim = (sliceRenderer.dim+1)%3;
          sliceRenderer.setVisible(true);
          sliceRenderer.invalidate();
          }
        else
          {
          sliceRenderer.sliceNo = 0;
          sliceRenderer.dim = (cmd=="xy") ? 2 : (cmd=="xz") ? 1 : 0;
          sliceRenderer.setVisible(true);
          sliceRenderer.sliceMoving = true;
          sliceRenderer.invalidate();
          }
        }
      else if (cmd == "cloud")
        {
        levelSetRenderer.setVisible(!levelSetRenderer.isVisible());
        }
      else if (cmd == "medial")
        {
        iss >> cmd;
        if (cmd == "seethrough" || cmd == "seethru")
          {
          int mode = rndSpline->getSurfaceMode();
          mode = (mode & (!BSplineRenderer::SURFACE)) | BSplineRenderer::SEETHRU;
          rndSpline->setSurfaceMode(mode);
          }
        if (cmd == "solid" || cmd == "flat")
          {
          int mode = rndSpline->getSurfaceMode();
          mode = (mode & (!BSplineRenderer::SEETHRU)) | BSplineRenderer::SURFACE;
          rndSpline->setSurfaceMode(mode);
          }
        }
      else if (cmd == "boundary")
        {
        iss >> cmd;
        if (cmd == "top")
          {
          rndSpline->setSurfaceMode(rndSpline->getSurfaceMode() ^ rndSpline->IBSURFACE_LEFT);
          }
        if (cmd == "bot" || cmd == "bottom")
          {
          rndSpline->setSurfaceMode(rndSpline->getSurfaceMode() ^ rndSpline->IBSURFACE_RIGHT);
          }
        }
      else if (cmd == "control")
        {
        rndSpline->setSurfaceMode(rndSpline->getSurfaceMode() ^ rndSpline->CONTROL);
        }
      else if (cmd == "map")
        {
        iss >> cmd;
        if (cmd == "gauss")
          {
          rndSpline->setColorMapMode(BSplineRenderer::CMAP_GAUSS);
          }
        else if (cmd == "off")
          {
          rndSpline->setColorMapMode(BSplineRenderer::CMAP_NONE);
          }
        else if (cmd == "next")
          {
          rndSpline->setColorMapMode((BSplineRenderer::ColorModes)((rndSpline->getColorMapMode() + 1) % BSplineRenderer::CMAP_COUNT));
          }
        }
      else if (cmd == "edge")
        {
        edgeNetRenderer.setVisible(!edgeNetRenderer.isVisible());
        }

      glutPostRedisplay();
      }

    else if (cmd == "new")
      {
      iss >> cmd; 
      if (cmd == "spline")
        {
        int m, n;
        iss >> m;
        iss >> n;

        newSpline(m,n);
        }
      }

    else if (cmd == "scale")
      {
      iss >> cmd; 
      if (cmd == "spline")
        {
        float sx, sy, sz, sr;
        iss >> sx;
        iss >> sy;
        iss >> sz;
        iss >> sr;

        scaleSpline(sx,sy,sz,sr);
        }
      else if (cmd == "image")
        {
        float sx, sy, sz;
        iss >> sx;
        iss >> sy;
        iss >> sz;

        if (imaging.getImage())
          {
          imaging.getImage()->setVoxelSize(sx,sy,sz);
          onImageChange();
          }
        }
      }

    else if (cmd == "rotate")
      {
      iss >> cmd; 
      if (cmd == "spline")
        {
        float vx, vy, vz, angle;
        iss >> vx;
        iss >> vy;
        iss >> vz;
        iss >> angle;

        rotateSpline(vx,vy,vz,angle);
        }
      }

    else if (cmd == "shift")
      {
      iss >> cmd; 
      if (cmd == "spline")
        {
        float sx, sy, sz;
        iss >> sx;
        iss >> sy;
        iss >> sz;

        shiftSpline(sx,sy,sz);
        }
      }

    else if (cmd == "crop")
      {
      iss >> cmd;
      if (cmd == "image")
        {
        float pad;
        int resample;
        iss >> pad;
        iss >> resample;
        if (imaging.getType() == imaging.BINARY)
          {
          BinaryImage *bim = (BinaryImage*)imaging.getImage();
          bim->autoCrop(pad,resample);
          onImageChange();
          }
        else if (imaging.getType() == imaging.GRAYSCALE)
          {
          int bkg = 0;
          iss >> bkg;
          GrayImage *im = (GrayImage*)imaging.getImage();

          im->autoCrop(pad,resample,bkg);
          onImageChange();
          }
        }
      }

    else if (cmd == "threshold")
      {
      iss >> cmd;
      if (cmd == "image")
        {
        if (imaging.getType() == imaging.BINARY)
          {
          BinaryImage *bim = (BinaryImage*)imaging.getImage();
          bim->threshold(0);
          onImageChange();
          }
        else if (imaging.getType() == imaging.GRAYSCALE)
          {
          GrayImage *im = (GrayImage*)imaging.getImage();

          int intensity = 256*128;
          iss >> intensity;

          im->threshold(intensity);
          onImageChange();
          }
        }
      }

    else if (cmd == "match" || cmd == "m")
      {
      matchSpline(sout);
      }

    else if (cmd == "optimize" || cmd == "opt")
      {

      // Optimize mode
      if (modeManager->getActiveMode() == "optimize")
        {
        modeManager->setMode("null");
        }
      else
        {
        modeManager->setMode("optimize");
        }
      }

    else if (cmd == "rigid")
      {

      // Optimize mode
      if (modeManager->getActiveMode() == "rigid")
        {
        modeManager->setMode("null");
        }
      else
        {
        modeManager->setMode("rigid");
        }
      }

    else if (cmd == "set")
      {
      iss >> cmd;

      if (cmd == "set")
        {
        list<string> lst;
        settings->collectAllKeys(lst);
        for (list<string>::iterator it=lst.begin();it!=lst.end();it++)
          sout << *it << " = " << settings->getStringValue(it->c_str(),"undefined") << endl;
        }
      else if (settings->isFolder(cmd.c_str()))
        {
        Registry *reg = &settings->getSubFolder(cmd.c_str());
        list<string> lst;
        reg->collectAllKeys(lst);
        for (list<string>::iterator it=lst.begin();it!=lst.end();it++)
          sout << *it << " = " << reg->getStringValue(it->c_str(),"undefined") << endl;
        }
      else
        {
        string value;
        iss >> value;
        settings->setStringValue(cmd.c_str(),value.c_str());
        sout << cmd << " = " << settings->getStringValue(cmd.c_str(),"undefined") << endl;
        }
      }

    else if (cmd == "align")
      {
      
        if (imaging.getType() == imaging.BINARY)
          {
          SplineBinaryVolumeMatcher *match = (SplineBinaryVolumeMatcher *)imaging.getImageMatch();    
          match->setSpline(spline,splineData);
          match->alignSplineByMoments();
          }
      }

    else if (cmd == "test")
      {
      iss >> cmd;


      if (cmd == "flow")
        {
        testImageDerivative();
        testDerivatives();
        flowTester->update(sout);
        }

      if (cmd == "fastflow")
        {
        flowTester->computeControlGradient(sout);
        }

      if (cmd == "flowmove")
        {
        flowTester->apply();
        }

      if (cmd == "gl")
        {
        sout << "GL  Vendor     : " << glGetString(GL_VENDOR) << endl;
        sout << "GL  Renderer   : " << glGetString(GL_RENDERER) << endl;
        sout << "GL  Version    : " << glGetString(GL_VERSION) << endl;
        sout << "GL  Extensions : " << endl << glGetString(GL_EXTENSIONS) << endl;
        sout << "GLU Version    : " << gluGetString(GLU_VERSION) << endl;
        sout << "GLU Extensions : " << endl << gluGetString(GLU_EXTENSIONS) << endl;
        }

      else if (cmd == "trilerp")
        {
        AbstractImage3D *img = imaging.getImage();
        if (img)
          {
          vector<SMLVec3f> r;
          for (int i=0;i<10000;i++)
            {
            r.push_back(SMLVec3f(rand(1.0),rand(1.0),rand(1.0)));
            }
          double t0 = clock();
          for (int j=0;j<10000;j++)
            {
            for (int i=0;i<10000;i++)
              {
              SMLVec3f &z = r[i];
              img->interpolateVoxel(z[0],z[1],z[2]);
              }
            }
          t0 = clock() - t0;
          sout << "ten million interpolations in " << t0 << " ms." << endl;
          }
        }

      else if (cmd == "crash")
        {
        try
          {
          char *p = NULL;
          strcpy(p,"crash");
          }
        catch (...)
          {
          sout << "crash" << endl;
          }
        }

      else if (cmd == "paste")
        {
#ifdef _MSC_VER
        if (OpenClipboard(NULL))
          {
          HGLOBAL hglb = GetClipboardData(CF_TEXT);
          char *text = (char *)hglb;

          sout << text << endl;
          CloseClipboard();
          }
#endif // _MSC_VER
        }
      }

    else if (cmd == "sleep")
      {
      int ms = 0;
      iss >> ms;
      sout << "Halting script execution for " << ms << " ms." << endl;
      scriptProcessor.sleep(ms);
      }

    else if (cmd == "quit")
      {
      exitProcedure(0);           
      }

    else if (cmd == "matlab")
      {
      /*
      if(!matlab) {
          matlab = new Matlab();
      }
      if(matlab->connected()) {
          matlabMode = true;
          sout << "Entering MATLAB mode" << endl;
      }
      else {
          sout << "MATLAB not available" << endl;
      }*/
      sout << "MATLAB not available" << endl;
      }

    else if (cmd.substr(0,1)=="!")
      {
      callShellCommand(iss.str().substr(1),sout);
      }

    }
  catch (char *exc)
    {
    sout << "Invalid command syntax: " << exc << endl << "Type 'help' for some help!" << endl;
    }
}

void makeKeyBindings() {
  int ctl = GLUT_ACTIVE_CTRL;
  int alt = GLUT_ACTIVE_ALT;

  globalEventHandler.addKeyBinding('0',alt,"show medial seethrough");
  globalEventHandler.addKeyBinding('1',alt,"show boundary top");
  globalEventHandler.addKeyBinding('2',alt,"show boundary bottom");   
  globalEventHandler.addKeyBinding('c',alt,"show control");   
  globalEventHandler.addKeyBinding('g',alt,"show map next");  
  globalEventHandler.addKeyBinding('f',alt,"test flow");  
  globalEventHandler.addKeyBinding('a',alt,"test flowmove");  
  globalEventHandler.addKeyBinding('s',alt,"test fastflow");  
  globalEventHandler.addKeyBinding('p',alt,"show slice pause");   

  globalEventHandler.addKeyBinding('m',alt,"match");  
  globalEventHandler.addKeyBinding('o',alt,"optimize");   
  globalEventHandler.addKeyBinding(GLUT_KEY_F10,alt,"quit");

  globalEventHandler.addKeyBinding(GLUT_KEY_PAGE_UP,0,"show slice next");
  globalEventHandler.addKeyBinding(GLUT_KEY_PAGE_DOWN,0,"show slice last");
  globalEventHandler.addKeyBinding(GLUT_KEY_END,0,"show slice axis");
}

void exitProcedure(int mode) {
  // Clean up code
  cout << "Cleaning Up\n";

  // Destroy the spline and associated stuff
  destroySpline();

  // Shut down the display driver
  GLDisplayDriver::shutdown();

  // Save and destroy the registry
  // settings->writeToFile("mspline.ini");
  // delete settings;
  settingsAll.writeToFile("mspline.ini");

  // Close matlab
  //if(matlab)
  //    delete matlab;

  // Exit
  exit(mode);
}

#include <Powell.h>
#include <problems.h>

void testOptimization() {
  Vector center(10);
  center.setAll(1.0);
  DeJongF4Pbm pr(center);
  //SixHumpCamelBackPbm pr;

  Vector guess(10);
  for (int i=0;i<guess.size();i++)
    guess(i) = rand(10.0) - 5.0;

  int iter = 0;
  PowellMethod pm(&pr,guess);
  ConjugateGradientMethod cgm(pr,guess);
  while (!pm.isFinished())
    {
    pm.performIteration();
    cgm.performIteration();
    cout << (++iter) << "\t" << pm.getBestEverValue() << "\t" << cgm.getBestEverValue() << "\t";

    pm.getBestEverX().t().print();
    }
}

bool ScriptProcessor::handleIdle() {
  // Only proceed if we are not in waiting mode
  if (tWaitEnd < clock())
    {
    // Make sure there are more commands
    if (commands.size() == 0)
      {
      if(flagInteractive)
        GLDisplayDriver::removeListener(this,GLDisplayDriver::IDLE);
      return true;
      }

    // Pop a command from the stack
    string command = commands.front();
    commands.pop_front();
    executeCommand(command);

    // Automatically sleep 3 seconds
    if (tWaitEnd < clock())
      {
      sleep(250);
      }
    }
  return false;
}

string getenvs(string var) {
  char *ev = getenv(var.c_str());
  string rtn;
  rtn.assign(ev ? ev : "");
  return rtn;
}


// Script commands can use environment variables as '{'
void ScriptProcessor::loadScript(const char *fname) {
  ifstream ifs(fname);
  commands.clear();
  string command;

  while (getline(ifs,command))
    {
    string newcmd,evar;
    bool evmode = false;

    for (int i=0;i<command.length();i++)
      {
      char ch = command.at(i);            
      if (evmode)
        {
        if(ch == '$')
          {
          newcmd += '$';
          evmode = false;
          }
        if(ch == '{')
          {
          evar.assign("");
          }
        else if(ch == '}')
          {
          newcmd += getenvs(evar);
          evmode = false;
          }
        else
          {
          evar += ch;
          }
        }
      else
        {
        if (ch == '$')
          {
          evmode = true;
          }
        else if (ch == '#' && i==0)
          {
          break;
          }
        else
          {
          newcmd += ch;
          }
        }
      }

    if (evmode && evar.length() > 0)
      {
      newcmd += getenvs(evar);
      }

    // Add command to script
    if (newcmd.length() > 0)
      {
      commands.push_back(newcmd);
      }
    }

  ifs.close();
}


// Processing commands from script files

int usage()
  {
  const char *usage = 
    "PROGRAM: cmreps3d \n"
    "  usage: cmreps3d [options] \n"
    "options: \n"
    "  -s, --script FILE      Script file to execute upon starting the program\n"
    "  -i, --init FILE        Location of the initialization file (curr. dir.)\n"
    "  -n, --nogui            Run program in non-interactive mode\n"
    "  -h, --help             Bring up this help message\n";

    cout << usage;
    return -1;
  }

int main(int argc,char *argv[]) {

  // Set memory debugging
  // _CrtSetDbgFlag(_CrtSetDbgFlag(_CRTDBG_REPORT_FLAG) | _CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF);
  // _CrtDumpMemoryLeaks();
  float *test = (float*)_aligned_malloc(sizeof(float) * 4, 256);
  __m128 junk;
  _mm_store_ps(test,junk);

  // Read the command line parameters
  string sInitFile("mspline.ini"), sScriptFile;

  // Read the command line parameters and see if there is a script to execute
  for (int ar=1;ar<argc;ar++)
    {
    string arg = argv[ar];
    if (arg == "--script" || arg == "-s")
      { sScriptFile = argv[++ar]; }
    else if(arg == "--init" || arg == "-i")
      { sInitFile = argv[++ar]; }
    else if(arg == "--nogui" || arg == "-n")
      { flagInteractive = false; }
    else 
      { return usage(); }
    }

  // Load the initialization settings
  try
    {
    settingsAll.readFromFile(sInitFile.c_str());
    }
  catch (RException exc)
    {
    cerr << "Failed to read initialization settings from " << sInitFile << endl;
    }
  settings = &settingsAll.getSubFolder("user");
  settingsUI = &settingsAll.getSubFolder("ui");
  settings->setFlagAddIfNotFound(true);

  // Load the script file
  if(sScriptFile != "")
    scriptProcessor.loadScript(sScriptFile.c_str());

  // Build and display a spline
  spline = new MSpline(9,6);
  splineData = new SplineDataCache(spline,RESLN);
  rndSpline = new BSplineRenderer(spline,splineData,RESLN);
  rndSpline->setSurfaceMode(rndSpline->SEETHRU);

  // If in GUI mode, start all the interactive stuff
  if( flagInteractive ) 
    {
    // Initialize display system
    GLDisplayDriver::init(argc,argv);

    // Add command buttons
    initModes();

    // Add a global event listener at the back of the food chain
    GLDisplayDriver::addListener(&globalEventHandler,GLDisplayDriver::SPECIAL | GLDisplayDriver::KEYS);
    makeKeyBindings();

    GLDisplayDriver::addRenderer(new LightRenderer());
    GLDisplayDriver::addRenderer(rndSpline);

  #ifndef COLOR_SET_1
    GLDisplayDriver::addRenderer(new StarfieldRenderer());
  #endif

    GLDisplayDriver::addRenderer(rndSpline,GLDisplayDriver::TRAN);
    GLDisplayDriver::addRenderer(new FrameRateCountRenderer(),GLDisplayDriver::EYE);
    //GLDisplayDriver::addRenderer(&sliceRenderer,GLDisplayDriver::EYE);
    GLDisplayDriver::addRenderer(&sliceRenderer);
    GLDisplayDriver::addRenderer(&levelSetRenderer);
    GLDisplayDriver::addRenderer(&edgeNetRenderer);

    GLDisplayDriver::addRenderer(&uvImage,GLDisplayDriver::EYE);

    // Add a flow tester
    flowTester = new FlowTester();
    GLDisplayDriver::addRenderer(flowTester);
    GLDisplayDriver::addRenderer(&rndFadeStatus,GLDisplayDriver::EYE);
    rndFadeStatus.setAnchor(true,5,-24);    
    rndFadeStatus.addLine("Welcome to MSpline.  Press F1 for help!");

    // Create the command shell
    fntCourier12->printf("test");
    commandShell = new CommandShell();
    GLDisplayDriver::addRenderer(commandShell,GLDisplayDriver::EYE);
    commandShell->setVisible(false);

    // Load the script parameter
    if(sScriptFile != "")
      {
      // Add the script processor as the idle function
      GLDisplayDriver::addListener(&scriptProcessor,GLDisplayDriver::IDLE);

      // Make the command shell visible
      commandShell->setVisible(true);
      GLDisplayDriver::addListener(commandShell,GLDisplayDriver::KEYS | GLDisplayDriver::BUTTON | GLDisplayDriver::SPECIAL); 
      }

    // Start GLUT
    glutMainLoop();
    }
  else if(sScriptFile == "")
    {
    cerr << "Can not run in non-interactive mode without a script" << endl;
    return -1;
    }
  else
    {
    // Add command buttons
    initModes();

    bool scriptDone = false;
    while(!scriptDone) 
      {
      GLDisplayDriver::runBackgroundProcesses();
      scriptDone = scriptProcessor.handleIdle();
      }
    }

  return 0;
}

