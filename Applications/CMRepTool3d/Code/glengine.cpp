#include "glengine.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cstdarg>
#include "colors.h"

#ifdef _MSC_VER
#include <windows.h>
#endif

using namespace std;

// This method is slower than the above because it does not use lookup tables,
// but it is more precise and scientifically sound
double getGaussianRnd(double mean,double sigma) {
  const double factor = 1.0 / RAND_MAX;
  const double m_PI = acos(-1.0);

  double u1 = factor * ::rand();
  double u2 = factor * ::rand();
  return mean + 2 * sigma * sqrt(-0.5*log(u1)) * cos(2*m_PI*u2);
}

double rand(double max) {
  const double factor = 1.0 / RAND_MAX;
  return max * factor * rand();
}



/*********************************************************
 Display Driver Code
  *******************************************************/
// int GLDisplayDriver::dlAxes;
int GLDisplayDriver::width;
int GLDisplayDriver::height;
int GLDisplayDriver::dragX;
int GLDisplayDriver::dragY;
int GLDisplayDriver::dragButton;
clock_t GLDisplayDriver::tLastRefresh;

SMLVec3f GLDisplayDriver::center(0.0,0,0);
float GLDisplayDriver::scale = 1.0f;

GLDisplayDriver::Modes GLDisplayDriver::displayMode = GLDisplayDriver::NORM;

list <GLRenderer *> GLDisplayDriver::rnd[GLDisplayDriver::INUMMODES];
list<GLEventListener *> GLDisplayDriver::listeners[GLDisplayDriver::INUMEVENTS];

// FracViewer event listener
class AGVEventHandler : public GLEventListener
  {
  virtual bool handleButton(int button,int state,int x,int y) {
    agvHandleButton(button,state,x,y);
    return true;
  }

  virtual bool handleMotion(int x,int y) {
    agvHandleMotion(x,y);
    return true;
  }

  virtual bool handleKeys(unsigned char key, int x, int y) {
    agvHandleKeys(key,x,y);
    return true;
  }
  };

AGVEventHandler agvEventHandler;

void GLDisplayDriver::reshape(int w,int h) {
  width = w;
  height = h;

  glViewport(0,0,w,h);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0, (GLdouble)w/h, 0.01, 100000);
  glPushMatrix();

  glMatrixMode(GL_MODELVIEW);
  glFlush();
};

extern GLfloat EyeAz,EyeDist,EyeEl;

void GLDisplayDriver::worldTransform() {
  // Set up the projection component: apply the agv transform
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glPushMatrix();  
  agvViewTransform();

  // Draw the model
  glMatrixMode(GL_MODELVIEW); 
}

void GLDisplayDriver::drawMode(IModes idx) {
  displayMode = (Modes)(1 << idx);

  list<GLRenderer*>::iterator it;
  for (it=rnd[idx].begin();it!=rnd[idx].end();it++)
    if ((*it)->isVisible())
      (*it)->onDraw();
}

void GLDisplayDriver::draw() {
  // Save the last refresh time
  tLastRefresh = clock();

  // Set up the projection component: apply the agv transform
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glPushMatrix();
  agvViewTransform();

  // Enter model view matrix mode
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Draw the preliminary mode
  drawMode(IPRE);

  // Clear the contents
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Draw the objects in the untransformed space
  drawMode(IUNIT);

  // Push space matrix
  glPushMatrix();
  glScaled(1.0f / scale, 1.0f / scale, 1.0f / scale);
  glTranslated(-center[0],-center[1],-center[2]);

  // Draw all renderers
  drawMode(INORM);

  // Make the depth buffer read only to render transparent objects
  if (rnd[1].size())
    {
    glEnable(GL_BLEND);
    //glDepthMask(GL_FALSE);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Draw all see through renderers
    drawMode(ITRAN);

    //glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
    }

  // Pop space matrix
  glPopMatrix();

  // Transform the view into eye coordinates
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0.0,width,0.0,height);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  drawMode(IEYE);

  // Flush the display
  glFlush();
  glutSwapBuffers();
}

void GLDisplayDriver::handleButton(int button,int state,int x,int y) {
  // Save the drag coordinates
  if (state == GLUT_DOWN)
    {
    dragButton = button;
    dragX = x;
    dragY = y;
    }

  list<GLEventListener*>::iterator it;
  for (it=listeners[IBUTTON].begin();it!=listeners[IBUTTON].end();it++)
    {
    GLEventListener *lst = *it;
    if (lst->handleButton(button,state,x,y))
      return;
    }
}

void GLDisplayDriver::handleMotion(int x,int y) {

  list<GLEventListener*>::iterator it;
  for (it=listeners[IMOTION].begin();it!=listeners[IMOTION].end();it++)
    {
    GLEventListener *lst = *it;
    if (lst->handleMotion(x,y))
      break;
    }

  // Update the drag coordinates
  dragX = x;
  dragY = y;
}

void GLDisplayDriver::handlePassiveMotion(int x,int y) {
  list<GLEventListener*>::iterator it;
  for (it=listeners[IPASSIVE].begin();it!=listeners[IPASSIVE].end();it++)
    {
    GLEventListener *lst = *it;
    if (lst->handlePassiveMotion(x,y))
      return;
    }
}

void GLDisplayDriver::handleKeys(unsigned char key,int x,int y) {
  list<GLEventListener*>::iterator it;
  for (it=listeners[IKEYS].begin();it!=listeners[IKEYS].end();it++)
    {
    GLEventListener *lst = *it;
    if (lst->handleKeys(key,x,y))
      return;
    }
}

void GLDisplayDriver::handleSpecial(int key,int x,int y) {
  list<GLEventListener*>::iterator it;
  for (it=listeners[ISPECIAL].begin();it!=listeners[ISPECIAL].end();it++)
    {
    GLEventListener *lst = *it;
    if (lst->handleSpecial(key,x,y))
      return;
    }
}

void GLDisplayDriver::runBackgroundProcesses() 
{
  // If there are any idlers, call them
  list<GLEventListener*>::iterator it;
  for (it=listeners[IIDLE].begin();it!=listeners[IIDLE].end();it++)
    {
    GLEventListener *lst = *it;
    if (lst->handleIdle())
      break;
    }
}

void GLDisplayDriver::handleIdle() 
{
  // Do backround processing
  runBackgroundProcesses();
  
  // Make sure we refresh the screen regularly (30 frames per second);
  if (clock() - tLastRefresh > 20)
    {
    // Call the virtual trackball's idle routine
    if (agvMoving)
      {
      agvMove();
      }
    glutPostRedisplay();
    }
}

void GLDisplayDriver::addRenderer(GLRenderer *renderer, int renderPoints) {
  int bit = 1;
  for (int i=0;i<INUMMODES;i++)
    {
    if (renderPoints & bit)
      rnd[i].push_back(renderer);
    bit = bit << 1;
    }
}

void GLDisplayDriver::removeRenderer(GLRenderer *renderer, int renderPoints) {
  int bit = 1;
  for (int i=0;i<INUMMODES;i++)
    {
    if (renderPoints & bit)
      {
      list<GLRenderer*>::iterator it = find(rnd[i].begin(),rnd[i].end(),renderer);
      if (it != rnd[i].end())
        rnd[i].erase(it);
      }
    bit = bit << 1;
    }
}

void GLDisplayDriver::addListener(GLEventListener *listener,int events) {
  for (int i=0;i<INUMEVENTS;i++)
    {
    if (events & (1 << i))
      {
      listeners[i].push_front(listener);
      }
    }

  // Only attach a passive listener if needed
  if (events & PASSIVE)
    {
    glutPassiveMotionFunc(GLDisplayDriver::handlePassiveMotion);
    }
}

void GLDisplayDriver::removeListener(GLEventListener *listener,int events) {
  list<GLEventListener*>::iterator it;
  for (int i=0;i<INUMEVENTS;i++)
    {
    if (events & (1 << i))
      {
      it = find(listeners[i].begin(),listeners[i].end(),listener);
      if (it != listeners[i].end())
        listeners[i].erase(it);
      }
    }

  // Only attach a passive listener if needed
  if (listeners[IPASSIVE].size() == 0)
    {
    glutPassiveMotionFunc(NULL);
    }
}

void GLDisplayDriver::init(int argc,char *argv[]) {
  // Initalize glut
  // win3d = new Fl_Window(20,20,800,800,"3D Display");
  // win3d->show();                       // glut will die unless parent window visible
  // win3d->begin();                      // this will cause Glut window to be a child
  // glutInitWindowSize(1024, 768);
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowSize(800, 600);
  glutInitWindowPosition(0,0);        // place it inside parent window
  glutCreateWindow("3D Display");
  // glutFullScreen();

  // win3d->resizable(glut_window);
  // win3d->end();

  // agvInit(1); /* 1 cause we don't have our own idle */

  // Last refresh time
  tLastRefresh = clock();

  // Set the listener functions
  glutMouseFunc(GLDisplayDriver::handleButton);
  glutMotionFunc(GLDisplayDriver::handleMotion);
  glutKeyboardFunc(GLDisplayDriver::handleKeys);
  glutSpecialFunc(GLDisplayDriver::handleSpecial);
  glutIdleFunc(GLDisplayDriver::handleIdle);
  glutReshapeFunc(GLDisplayDriver::reshape);
  glutDisplayFunc(GLDisplayDriver::draw);

  // Initialize spinner
  addListener(&agvEventHandler,BUTTON | MOTION | KEYS);
  agvInit(0); 

  // Make a display list for the axes
  // dlAxes = glGenLists(1);
  // agvMakeAxesList(dlAxes);

  // Run GL initialization code - establish a scene
  GLfloat light_ambient[] = { 0.4f, 0.4f, 0.4f, 1.0f};
  GLfloat light_diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f};
  GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f};
  GLfloat light_position[] = { 0.0f, 4.0f, 0.0f, 1.0f};

  GLfloat light1_ambient[] = { 0.3f, 0.3f, 0.3f, 1.0f};
  GLfloat light1_diffuse[] = { 0.3f, 0.3f, 0.3f, 1.0f};
  GLfloat light1_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f};
  GLfloat light1_position[] = { -1.0f, -1.0f, -2.0f, 0.0f};

  GLfloat lmodel_ambient[] = { 1.0f, 1.0f, 1.0f, 1.0f};

  glEnable(GL_LIGHTING);  
  glEnable(GL_LIGHT0);

  glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  // glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
/*
    glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
    glLightfv(GL_LIGHT1, GL_SPECULAR, light1_specular);
    glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
*/

  // glEnable(GL_LIGHT1);

  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);

  // glEnable(GL_NORMALIZE);

  glShadeModel(GL_SMOOTH);

  // glEnable(GL_COLOR_MATERIAL);
  // glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);   

  glClearColor(clrBack[0],clrBack[1],clrBack[2],1);

  //glEnableClientState(GL_NORMAL_ARRAY);
  //glEnableClientState(GL_VERTEX_ARRAY);

  // glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);


}

void GLDisplayDriver::shutdown() {
}

void GLDisplayDriver::addCommandIcon(GLCommandIcon *icon) {
  addRenderer(icon,EYE);
  addListener(icon,BUTTON);
}

void StarfieldRenderer::buildList() {
  // Array of point coordinates
  int d = 3;
  int dim[] = {0,500,600,650};
  int i;

  // Allocate arrays
  GLfloat *vtx = new GLfloat[dim[d]*3];
  GLbyte *clr = new GLbyte[dim[d]*3];

  // Different sizes    
  for (i=0;i<d;i++)
    {
    for (int j=dim[i];j<dim[i+1];j++)
      {
      // Graylevel
      double gLev = getGaussianRnd(0.6,0.2);

      // Tint - red, blue or yellow
      double rTint = getGaussianRnd(gLev,0.05);
      double yTint = getGaussianRnd(gLev,0.05);
      double bTint = getGaussianRnd(gLev,0.05);
      double r=gLev,g=gLev,b=gLev;

      if (rTint > yTint && rTint > bTint)
        {
        r = rTint;
        }
      else if (yTint > rTint && yTint > bTint)
        {
        r = yTint;
        g = yTint;
        }
      else
        {
        b = bTint;
        }

      double x = rand(100) - 50;
      double y = getGaussianRnd(0,15);
      double z = rand(100) - 50;

      clr[j*3] = (GLbyte)(r*255);
      clr[j*3+1] = (GLbyte)(g*255);
      clr[j*3+2] = (GLbyte)(b*255);
      vtx[j*3] = (GLfloat)x;
      vtx[j*3+1] = (GLfloat)y;
      vtx[j*3+2] = (GLfloat)z;
      }

    glEnd();
    }

  glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);

  glVertexPointer(3,GL_FLOAT,0,vtx);
  glColorPointer(3,GL_UNSIGNED_BYTE,0,clr);

  glNewList(dl,GL_COMPILE);

  glPushAttrib(GL_LIGHTING_BIT);
  glDisable(GL_LIGHTING);

  for (i=0;i<d;i++)
    {
    glPointSize(i+1);
    glDrawArrays(GL_POINTS,dim[i],dim[i+1]-dim[i]);
    }

  glPopAttrib();
  glEndList();

  glPopClientAttrib();

  delete vtx;
  delete clr;
}


struct tgaHeader
  {
  unsigned char blah[12];
  short width;
  short height;
  short blah1;
  };


unsigned int loadTransparentMipMap(const char *name,const GLColor &trans) {
  // Error code
  const unsigned int FAIL = (unsigned int) -1;
  
  // Store tex info   
  glPushAttrib(GL_TEXTURE_BIT);

  // Allocate a texture name
  GLenum tName = GLTextureFactory::factory.getTexture();

  tgaHeader header;

  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D,tName);
  glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER, GL_NEAREST_MIPMAP_NEAREST);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST);
  glPixelStorei(GL_UNPACK_ALIGNMENT,1);
  glPixelStorei(GL_UNPACK_ROW_LENGTH,0);
  glPixelStorei(GL_UNPACK_SKIP_ROWS,0);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS,0);

  FILE *f = fopen(name,"rb");
  if (f==NULL)
    return FAIL;

  fread(&header,18,1,f);

  if (header.width != header.height)
    {
    fclose(f);
    return FAIL;
    }

  int r = header.width;

  // Read the bytes
  unsigned char *array = new unsigned char[r*r*3];
  unsigned char *rgba = new unsigned char[r*r*4];
  fread(array,1,r*r*3,f);

  int ind1=0,ind2=0;
  int tR = (int)(255 * trans.fv[0]);
  int tG = (int)(255 * trans.fv[1]);
  int tB = (int)(255 * trans.fv[2]);

  for (int i=0;i<r*r;i++)
    {
    int r = array[ind1+2];
    int g = array[ind1+1];
    int b = array[ind1];

    // Luminance
    int l = (r > b) ? r : b;
    l = (l > g) ? l : g;


    rgba[ind2] = r;
    rgba[ind2+1] = g;
    rgba[ind2+2] = b;
    rgba[ind2+3] = l == 0 ? 0 : 255;

    ind1+=3;
    ind2+=4;
    }

  gluBuild2DMipmaps(GL_TEXTURE_2D,4,r,r,GL_RGBA,GL_UNSIGNED_BYTE,rgba);

  delete array;
  delete rgba;

  glPopAttrib();

  return tName;
}


unsigned int loadMipMap(const char *name,int startRes,int endRes,GLuint intFormat) {
  // Store tex info   
  glPushAttrib(GL_TEXTURE_BIT);

  // Allocate a texture name
  GLenum tName = GLTextureFactory::factory.getTexture();

  char fname[256];
  tgaHeader header;

  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D,tName);
  glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER, GL_NEAREST_MIPMAP_NEAREST);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST);
  glPixelStorei(GL_UNPACK_ALIGNMENT,1);
  glPixelStorei(GL_UNPACK_ROW_LENGTH,0);
  glPixelStorei(GL_UNPACK_SKIP_ROWS,0);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS,0);

  int frame = 0;
  for (int r=endRes;r>=startRes;r/=2)
    {
    sprintf(fname,"%s_%d.tga",name,r);
    FILE *f = fopen(fname,"rb");
    if (f==NULL)
      continue;

    fread(&header,18,1,f);
    if (header.width!=r || header.height!=r)
      {
      fclose(f);
      continue;
      }

    // Read the bytes
    unsigned char *array = new unsigned char[r*r*3];
    fread(array,1,r*r*3,f);

    // Create a bitmap
    if (frame==0)
      gluBuild2DMipmaps(GL_TEXTURE_2D,intFormat,r,r,GL_BGR_EXT,GL_UNSIGNED_BYTE,array);
    else
      glTexImage2D(GL_TEXTURE_2D,frame,intFormat,r,r,0,GL_BGR_EXT,GL_UNSIGNED_BYTE,array);
    delete array;

    frame++;
    }

  glPopAttrib();

  cout << "Texture " << tName << 
  "\tIsResident = " << GLTextureFactory::factory.isResident(tName) << 
  "\tMipmap = " << name << "\n";

  return tName;
}

/*******************
 * Texture Factory
 ***************/
GLTextureFactory::GLTextureFactory() {

}

GLenum GLTextureFactory::getTexture() {
  GLenum i;
  glGenTextures(1,&i);
  names.push_back(i);

  return i;
}

bool GLTextureFactory::isResident(GLenum tex) {
  GLint i;
  glGetTexParameteriv(tex,GL_TEXTURE_RESIDENT,&i);
  return i != 0;
}

void GLTextureFactory::apply2D(GLenum tex) {
  glBindTexture(GL_TEXTURE_2D,tex);
}

void GLTextureFactory::release(GLenum tex) {
  glDeleteTextures(1,&tex);
  names.remove(tex);
}

GLTextureFactory::~GLTextureFactory() {
  list<GLenum>::iterator i;
  for (i=names.begin();i!=names.end();i++)
    {
    GLenum tex = *i;
    glDeleteTextures(1,&tex);
    }
}

GLTextureFactory GLTextureFactory::factory;

/*******************
 * Icon
 ***************/
GLIcon::GLIcon(const char *fname) {
  file = fname;
  tn = 0;
}


GLIcon::~GLIcon() {
  glDeleteTextures(1,&tn);        
};

void GLIcon::onDraw() {
  if (tn == 0)
    tn = loadTransparentMipMap(file.c_str(),GLColor(0,0,0));

  glPushAttrib(GL_TEXTURE_BIT | GL_COLOR_BUFFER_BIT);

  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D,tn);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glBegin(GL_QUADS);

  glTexCoord2d(0.0,1.0);
  glVertex2d(0.0,0.0);
  glTexCoord2d(0.0,0.0);
  glVertex2d(0.0,1.0);
  glTexCoord2d(1.0,0.0);
  glVertex2d(1.0,1.0);
  glTexCoord2d(1.0,1.0);
  glVertex2d(1.0,0.0);

  glEnd();

  glPopAttrib();
};


GLCommandIcon::GLCommandIcon(const char *fnActive,const char *fnPassive,float x,float y,float w,float h) :
icPassive(fnPassive),icActive(fnActive)
{
  this->x = x;
  this->y = y;
  this->w = w;
  this->h = h;
  clickState = false;
}

bool GLCommandIcon::handleButton(int button,int state,int x,int y) {
  // float x1 = 2.0 * x / GLDisplayDriver::width - 1.0;
  // float y1 = 1.0 - 2.0 * y / GLDisplayDriver::height;
  float x1 = x;
  float y1 = GLDisplayDriver::height - y;
  float x0 = this->x >= 0 ? this->x : GLDisplayDriver::width + this->x;
  float y0 = this->y >= 0 ? GLDisplayDriver::height - this->y : - this->y;

  //cout << x << ", " << y << "\n";
  //cout << x0 << ", " << y0 << "\n";

  if (state==GLUT_DOWN)
    {
    if (x0 <= x1 && x0+w >= x1 && y0 <= y1 && y0+h >= y1)
      {
      // We are in the icon square
      cout << "hit!\n";

      clickState = !clickState;
      // invalidate();
      onStateChange();

      glutPostRedisplay();

      return true;
      }
    }

  return false;
}

void GLCommandIcon::onDraw() 
{
  glPushAttrib(GL_LIGHTING_BIT);
  glPushMatrix();

  glDisable(GL_LIGHTING);

  glTranslated(x >= 0 ? x : GLDisplayDriver::width + x,
               y >= 0 ? GLDisplayDriver::height-y : -y,0);
  glScaled(w,h,1);

  //glTranslated(400,400,0);
  //glScaled(200,200,1);

  glColor3d(1,1,1);

  if (clickState)
    {
    icPassive.onDraw();
    }
  else
    {
    icActive.onDraw();
    }

  glPopMatrix();
  glPopAttrib();
}

/**
 * Write text to the screen
 */
void stroke_output(GLfloat x, GLfloat y, GLfloat height,char *format,...)
{
  va_list args;
  char buffer[200], *p;

  va_start(args, format);
  vsprintf(buffer, format, args);
  va_end(args);

  glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT | GL_COLOR_BUFFER_BIT);
  glDisable(GL_LIGHTING);

  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glLineWidth(0.5);


  glPushMatrix();
  glTranslatef(x, y, 0.0f);
  glScaled(height/120.0, height/120.0, 0.0);
  for (p = buffer; *p; p++)
    glutStrokeCharacter(GLUT_STROKE_ROMAN, *p);
  glPopMatrix();

  glPopAttrib();
}

FrameRateCountRenderer::FrameRateCountRenderer() {
  for (int i=0;i<SPAN;i++)
    t[i] = clock();
  tIndex=0;
}

void FrameRateCountRenderer::onDraw() {
  // Compute the frame rate
  clock_t t0 = clock();
  if (t0 - t[tIndex] > 0)
    {
    int fr = (SPAN * CLOCKS_PER_SEC) / (t0 - t[tIndex]);
    stroke_output(GLDisplayDriver::width-40,25,14,"%d",fr);
    }
  t[tIndex] = t0;
  tIndex = (tIndex+1)%SPAN;
}

void DefaultLightRenderer::onDraw() {
  // Transform the view into eye coordinates
  glPushMatrix();

  glLoadIdentity();
  glRotatef(-EyeAz, 0, 1, 0);
  glRotatef(-EyeEl, 1, 0, 0);
  glTranslatef(0, 0, EyeDist);

  // Run GL initialization code - establish a scene
  GLfloat light_position[] = { 0.0f, 0.0f, 100.0f, 1.0f};
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);

  glPopMatrix();
}

const GLColor clrWhite(1);
const GLColor clrBlack(0);
const GLColor clrLightGray(0.75);
const GLColor clrGray(0.5);
const GLColor clrDarkGray(0.25);
const GLColor clrRed(1,0,0);
const GLColor clrGreen(0,1,0);
const GLColor clrBlue(0,0,1);


GLFont::GLFont(string title,int height,bool bold,bool italics) {
  this->title = string(title);
  this->height = height;
  this->bold = bold;
  this->italics = italics;
  dlBase = -1;
}

void GLFont::print(const char *text,unsigned int length) {

#ifdef _MSC_VER
  
  // Build the font on demand
  if (dlBase < 0)
    build();

  glPushAttrib(GL_LIST_BIT);                          // Pushes The Display List Bits

  glListBase(dlBase - 32);                                // Sets The Base Character to 32
  glCallLists(length, GL_UNSIGNED_BYTE, text);    // Draws The Display List Text

  glPopAttrib();                                      // Pops The Display List Bits

#else
  while(*text)
    {
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,*text);
    text++;
    }
#endif
}

void GLFont::printf(const char *fmt, ...)                   // Custom GL "Print" Routine
{
  char        text[1024];                             // Holds Our String
  va_list     ap;                                     // Pointer To List Of Arguments

  if (fmt == NULL)                                    // If There's No Text
    return;                                         // Do Nothing

  va_start(ap, fmt);                                  // Parses The String For Variables
  int length = vsprintf(text, fmt, ap);                           // And Converts Symbols To Actual Numbers
  va_end(ap);                                         // Results Are Stored In Text

  print(text,length);
}

void GLFont::build() {

#ifdef _MSC_VER
  HFONT   font;                                       // Windows Font ID

  dlBase = glGenLists(128);                               // Storage For 96 Characters

  font = CreateFont(  -height,                            // Height Of Font
                      0,                              // Width Of Font
                      0,                              // Angle Of Escapement
                      0,                              // Orientation Angle
                      bold?FW_BOLD:FW_NORMAL,         // Font Weight
                      italics,                        // Italic
                      FALSE,                          // Underline
                      FALSE,                          // Strikeout
                      ANSI_CHARSET,                   // Character Set Identifier
                      OUT_TT_PRECIS,                  // Output Precision
                      CLIP_DEFAULT_PRECIS,            // Clipping Precision
                      ANTIALIASED_QUALITY,            // Output Quality
                      FF_DONTCARE|DEFAULT_PITCH,      // Family And Pitch
                      title.c_str());                 // Font Name

  HDC hDC = wglGetCurrentDC();
  SelectObject(hDC, font);                            // Selects The Font We Want

  // Compute the widths
  GetCharWidth32(hDC,0,255,asciiCharWidths);

  cout << wglUseFontBitmaps(hDC, 32, 128, dlBase) << endl;            // Builds 96 Characters Starting At Character 32    
  cout << GetLastError() << endl;

#else // _MSC_VER

#endif // _MSC_VER
}

int GLFont::getHeight() {
  return height;  
}

int GLFont::getStringWidth(const char *text) {
  const char *p = text;
  int w = 0;
  while (*p)
    {
    w+=asciiCharWidths[*p];
    p++;
    }
  return w;
}

GLFont::~GLFont()
{
  glDeleteLists(dlBase, 128);                         // Delete All 96 Characters
}
