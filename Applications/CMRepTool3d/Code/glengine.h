#ifndef GL_ENGINE
#define GL_ENGINE

//#include <Fl\gl.h>
//#include <Fl\glut.h>
//#include <Fl\fl_file_chooser.h>
// #include <GL\gl.h>
// #include <GL\glu.h>
#include <GL/glut.h>
#include <smlmath.h>

#include <vector>
#include <list>
#include <string>
#include <ctime>

#include "fracviewer.h"

// A renderer
class GLRenderer {
private:
    bool visible;
public:
    // Called when the renderer is to be drawn
    virtual void onDraw() {};

    GLRenderer() {
        visible = true;
    }
    virtual ~GLRenderer() {}

    virtual void setVisible(bool visible) {
        this->visible = visible;
    }

    virtual bool isVisible() {
        return visible;
    }
};

// A display list renderer
class GLDisplayListRenderer : public GLRenderer {
protected:
    int dl;
    bool valid;

    GLDisplayListRenderer() {
        dl = -1;
        valid = false;
    }

    virtual ~GLDisplayListRenderer() {
        if(dl != -1)
            glDeleteLists(dl,1);
    }


    virtual void buildList() {
        glNewList(dl,GL_COMPILE);
        build();
        glEndList();
    }

    virtual void build() {};

public:
    // This method called to draw in GL
    virtual void onDraw() {
        if(dl == -1) {
            dl = glGenLists(1);
        }
        
        if(!valid) {
            buildList();
            valid = true;
        }

        glCallList(dl);
    }

    // Invalidate list for next display operation
    void invalidate() {
        valid = false;
    }
};

/**
 * Event listener classes 
 */
class GLEventListener {
protected:
    virtual bool handleButton(int button,int state,int x,int y) { return false; };
    virtual bool handleMotion(int x,int y) { return false; };
    virtual bool handlePassiveMotion(int x,int y) { return false; };
    virtual bool handleKeys(unsigned char key, int x, int y) { return false; };
    virtual bool handleSpecial(int key, int x, int y) { return false; };
    virtual bool handleIdle() { return false; };

    friend class GLDisplayDriver;
};

class GLColor {
public:
    GLfloat fv[4];
    
    GLColor(double r,double g,double b,double a=1.0) {
        fv[0] = r;
        fv[1] = g;
        fv[2] = b;
        fv[3] = a;
    }

    GLColor(double gLevel,double a=1.0) {
        fv[0] = gLevel;
        fv[1] = gLevel;
        fv[2] = gLevel;
        fv[3] = a;
    }

    GLColor() {
        fv[0] = fv[1] = fv[2] = 0;
        fv[3] = 1;
    }
};

extern const GLColor clrWhite;
extern const GLColor clrBlack;
extern const GLColor clrLightGray;
extern const GLColor clrGray;
extern const GLColor clrDarkGray;
extern const GLColor clrRed;
extern const GLColor clrGreen;
extern const GLColor clrBlue;

class GLIcon : public GLRenderer {
private:
    unsigned int tn;
    std::string file;
public:
    // Load the icon
    GLIcon(const char *fname);
    ~GLIcon();

    // Display the icon on the unit square
    void onDraw();
};

/**
 * This is a GL command icon
 */ 
class GLCommandIcon : public GLEventListener,public GLRenderer {
public:
    GLCommandIcon(const char *fnActive,const char *fnPassive,float x,float y,float w,float h);
    
    virtual void onStateChange() {};

    virtual void setState(bool state) {
        clickState = state;
        // invalidate();
    }

    void onDraw();

protected:
    // Passive and active states for the icon
    GLIcon icPassive,icActive;

    // Position and size of the icon in range -1 to 1
    float x,y,w,h;

    // State of the control
    bool clickState;

    // Handle the button press event
    bool handleButton(int button,int state,int x,int y);
};



class GLTextureFactory {
private:
    std::list<GLenum> names;
public:
    GLTextureFactory();
    GLenum getTexture();
    
    void apply2D(GLenum tex);
    void release(GLenum tex);
    bool isResident(GLenum tex);

    ~GLTextureFactory();

    static GLTextureFactory factory;
};

/*
class GLRenderGroup : public GLRenderer {
public:
    void onDraw();

    void addRenderer(GLRenderer *rnd);
    void removeRenderer(GLRenderer *rnd);

    void addTransparentRenderer(GLRenderer *rnd);
    void removeTransparentRenderer(GLRenderer *rnd);
    
protected:
    // Override this method to get a transform
    virtual void transform() {};


};
*/

// A display driver
class GLDisplayDriver {
private:
    // static int dlAxes;

    // Internal versions of the enums
    enum IModes {INORM = 0, ITRAN, IEYE, IPRE, INUMMODES};
    enum IEvents {IBUTTON = 0, IKEYS, IMOTION, IPASSIVE, ISPECIAL, IIDLE, INUMEVENTS};

    // List of renderers that are rendered in world coordinates
    static std::list <GLRenderer *> rnd[INUMMODES];

    // Event listener lists (for add/remove)
    static std::list<GLEventListener*> listeners[INUMEVENTS];

    // Draw all renderers in a given mode
    static void drawMode(IModes idx);

    // GLUT response functions
    static void reshape(int w,int h);
    static void draw();

    // Event handler routines
    static void handleButton(int button,int state,int x,int y);
    static void handleMotion(int x,int y);
    static void handlePassiveMotion(int x,int y);
    static void handleKeys(unsigned char key,int x,int y);
    static void handleSpecial(int key,int x,int y);
    static void handleIdle();

public:
    // Display size
    static int width,height;

    // Drag start coordinates
    static int dragX, dragY, dragButton;

    // Last refresh time
    static clock_t tLastRefresh;

    static void init(int argc,char *argv[]);
    static void shutdown();

    // This method adds a control icon to the display driver
    static void addCommandIcon(GLCommandIcon *icon);

    // Add renderer to the display driver's list
    static void addRenderer(GLRenderer *renderer,int points = NORM);
    static void removeRenderer(GLRenderer *renderer,int points = NORM);
    
    // Add and remove button listeners
    static void addListener(GLEventListener *listener,int events);
    static void removeListener(GLEventListener *listener,int events);

    // This method prepares a buffer for drawing
    static void worldTransform();

    // Center of the world
    static SMLVec3f center;

    // An enumeration of display modes
    enum Modes {NORM = 1, TRAN = 2, EYE = 4, PRE = 8};
    static Modes displayMode;

    // Enumeration of event listener types
    enum Events {BUTTON = 1, KEYS = 2, MOTION = 4, PASSIVE = 8, SPECIAL = 16, IDLE = 32};
};

// Star field renderer
class StarfieldRenderer : public GLDisplayListRenderer {
    void buildList();
};

// Frame rate counter
class FrameRateCountRenderer : public GLRenderer {
private:
    enum {SPAN=10};
    
    clock_t t[SPAN];
    int tIndex;
public: 
    void onDraw();
    FrameRateCountRenderer();
};

void stroke_output(GLfloat x, GLfloat y, GLfloat height,char *format,...);

/***********************************
 * A utility class to draw fonts
 ***********************************/
class GLFont {
private:
    std::string title;
    int height;
    bool bold,italics;
    
    int dlBase;
    void build();

    // Widths of characters
    int asciiCharWidths[256];

public:
    GLFont(std::string title,int height,bool bold=false,bool italics=false);
    ~GLFont();

    // Get the number of pixels used to represent a string
    int getStringWidth(const char *text);
    int getHeight();

    void print(const char *text, unsigned int length);
    void printf(const char *fmt, ...);
};

class GLMaterial {
    int dl;
public:
    GLColor ambient,diffuse,specular,emission;
    GLfloat shininess;
    GLint face;

    GLMaterial(int fce = GL_FRONT_AND_BACK,
             const GLColor &a = clrBlack,const GLColor &d = clrBlack,
             const GLColor &s = clrBlack,double sh=0,
             const GLColor &e = clrBlack) 
    {
        ambient = a;
        diffuse = d;
        specular = s;
        emission = e;
        shininess = sh;
        face = fce;
        dl = -1;
    }

    void apply() const {
        glMaterialfv(face,GL_AMBIENT,ambient.fv);
        glMaterialfv(face,GL_DIFFUSE,diffuse.fv);
        glMaterialfv(face,GL_SPECULAR,specular.fv);
        glMaterialfv(face,GL_EMISSION,emission.fv);
        glMaterialf(face,GL_SHININESS,shininess);
    }

    void buildDL() {
        dl = glGenLists(1);
        glNewList(dl,GL_COMPILE);
        apply();
        glEndList();
    }
    
    // Gets the id for the material
    int getDL() {
        if(dl < 0) 
            buildDL();
        return dl;
    }
};


#endif //GL_ENGINE