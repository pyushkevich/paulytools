#ifndef GLUI_H
#define GLUI_H

// My includes
#include "mspline.h"
#include "optimization.h"
// #include "ui.h"

// C++ includes
#include <iostream>
#include <ctime>
#include <map>
#include <algorithm>
#include <string>
#include <registry.h>

// Open GL includes
#include "glengine.h"

#ifndef M_PI
extern const double M_PI;
#endif

// Extenal classes for less inclusion
class BSplineRenderer;
class SplineOptDriver;
class SplineImageMatcher;

/**
 * A collection of index arrays at a specific level
 */
class SplineIndexArraySet {
private:
    //int dPatchNew;
    //int wQuadArray;
    //int hQuadArray;

    // The quad strip arrays, [row][col] 
    // int **idxQuad;

    // Array by patches, then by index
    Array2D< Array2D<unsigned int> > idxPatchQuad;

    // The start and count pairs.  This is necessary because some patches may have 'holes' where 
    // the shape 

public:
    SplineIndexArraySet(int level,int maxLevel,int m,int n);
    ~SplineIndexArraySet();
    
    // int *getQuadStripArray(int row); 
    unsigned int *getPatchQuadStripArray(int iPatch,int jPatch,int row) {
        return idxPatchQuad(iPatch,jPatch).row(row);
    }
};

/*********************************************************
 BSpline Renderer - Renders Everything for a spline
  *******************************************************/
class BSplineRenderer : public GLRenderer {
public:
    // Constructor for this renderer
    BSplineRenderer(DynamicBSpline2D *spline,SplineDataCache *dataCache, int level);
    ~BSplineRenderer();

    // Find a UV pair at a location
    void findSampleUV(int x,int y,int &u,int &v);

    // Find a control point at a position
    void findControlPoint(int x,int y,int &ci,int &cj);

    enum Modes {
        SURFACE = 1,WIRE = 2,SEETHRU = 4, CONTROL = 8, IBSURFACE_LEFT = 16, IBSURFACE_RIGHT = 32, TRIM_CURVE = 64
    };

    // Color maps
    enum ColorModes {
        CMAP_NONE = 0, CMAP_GAUSS = 1, CMAP_KOEN = 2, CMAP_FLETCH = 3, CMAP_IMAGE = 4, CMAP_COUNT = 5
    };

    // Set the current display mode as a combination of mode enums
    void setSurfaceMode(int mode);

    // Set the color map mode
    void setColorMapMode(ColorModes mode);

    // Get the color map mode
    ColorModes getColorMapMode() {
        return colorMapMode;
    }

    // Ge the current surface mode
    int getSurfaceMode() {
        return surfMode;
    }

    // Materials for surface rendering in different modes
    GLMaterial matFlatSurface[3],matSeeThruSurface[3],matBoundary[3];

    // Functions to build the different display lists
    void buildWire(int dl);
    void buildControl(int dl);
    void buildUVIndex(int dl);
    void buildControlIndex(int dl);

private:
    // A spline that we are rendering
    DynamicBSpline2D *spline;

    // A spline grid definition
    SplineGridDefinition *splineGrid;

    // A spline cache of data
    SplineDataCache *dataCache;

    // Index arrays for each level
    SplineIndexArraySet **arraySet;

    // Maximum and minimum levels
    int minLevel,maxLevel;

    // Patch dimensions
    int pw,ph;

    // Grid dimensions
    int gw,gh;

    // Sizes of the patch for min, max levels
    int dPatchMin,dPatchMax;

    // A general on-draw method
    void onDraw();

    // Control point array from the spline
    Array2D<SMLVec3f> C;

    // Color index array for UV lookup
    Array2D< Array2D<SMLVec3f> > UV;
    
    // Current surface mode
    int surfMode;

    // Current color map mode
    ColorModes colorMapMode;

    void makeMaterials();

    // Create a color patch based on the color patch mode
    void makeBoundaryColorPatch(PatchDataCache *pdc,int side, int iPatch, int jPatch);

    // Boundary color maps
    Array2D< Aux2D<SMLVec3f> > cmBoundary[2];

    void buildUVIndexPatch(int dl,int iPatch,int jPatch);
    void buildSurfacePatch(int dl,int iPatch,int jPatch,GLMaterial materials[]);
    void buildBoundaryPatch(int dl,int iPatch,int jPatch,int side,GLMaterial materials[]);
    void buildPatchTrimCurve(int dl,int iPatch,int jPatch);

    /**
    * This is a renderer used to display just a single patch
    */
    class PatchRenderer : public GLRenderer {
    protected:
        // The display list for the patch
        int dl;
        
        // Timestamp from when this list was last built
        int ts;
    public:     
        // The parent renderer
        BSplineRenderer *parent;
        
        // The index of the patch
        int iPatch,jPatch;
        SplineDataCache::Types type;
        
        // Build method
        virtual void buildList() = 0;

        // INvalidate
        void invalidate() {
            ts = -1;
        }
        
        // Draw method: check if the data is dirty, then render it
        void onDraw() {
            if(parent->dataCache->getPatchTimeStamp(iPatch,jPatch,parent->minLevel,type) > ts) {
                if(dl < 0)
                    dl = glGenLists(1);
                buildList();
                ts = parent->dataCache->getPatchTimeStamp(iPatch,jPatch,parent->minLevel,type);
            }
            glCallList(dl);
        }

        // Init method because I don't want to mess with contstructors
        void init(BSplineRenderer *parent,int iPatch,int jPatch,SplineDataCache::Types type) {
            this->parent = parent;
            this->iPatch = iPatch;
            this->jPatch = jPatch;
            this->type = type;
            dl = -1;
            ts = -2;
        }
    };

    /**
    * This is a renderer used to display just a single patch
    */
    class Renderer : public GLRenderer {
    protected:
        // The display list for the patch
        int dl;
        
        // Timestamp from when this list was last built
        int ts;

    public:     
        // The parent renderer
        BSplineRenderer *parent;
        
        // Subclassed for list building
        virtual void buildList() = 0;
        
        // Draw method: check if the data is dirty, then render it
        void onDraw() {
            if(parent->dataCache->getTimeStamp(parent->minLevel,SplineDataCache::MEDIAL) > ts) {
                if(dl < 0)
                    dl = glGenLists(1);
                buildList();
                ts = parent->dataCache->getTimeStamp(parent->minLevel,SplineDataCache::MEDIAL);
            }
            glCallList(dl);
        }

        Renderer() {
            dl = -1;
            ts = -2;
        }
    };

    // Wirefrace renderer for the surface
    class Wire : public /*BSplineRenderer::*/Renderer {
    public:
        void buildList() {
            parent->buildWire(dl);
        }
    } rndWire;
    
    // UV Index renderer for the surface
    class ControlIndex : public /*BSplineRenderer::*/Renderer {
    public:
        void buildList() {
            parent->buildControlIndex(dl);
        }
    } rndControlIndex;
    
    // UV Index renderer for the surface
    class Control : public /*BSplineRenderer::*/Renderer {
    public:
        void buildList() {
            parent->buildControl(dl);
        }
    } rndControl;   

    // UV Index renderer for the surface
    class PatchFlat : public /*BSplineRenderer::*/PatchRenderer {
    public:
        void buildList() {
            parent->buildSurfacePatch(dl,iPatch,jPatch,parent->matFlatSurface);
        }
    };
    Array2D<PatchFlat> rndPatchFlat;

    // UV Index renderer for the surface
    class PatchSeeThrough : public /*BSplineRenderer::*/PatchRenderer {
    public:
        void buildList() {
            parent->buildSurfacePatch(dl,iPatch,jPatch,parent->matSeeThruSurface);
        }
    };
    Array2D<PatchSeeThrough> rndPatchSeeThrough;

    // UV Index renderer for the surface
    class PatchTrimCurve : public /*BSplineRenderer::*/PatchRenderer {
    public:
        void buildList() {
            parent->buildPatchTrimCurve(dl,iPatch,jPatch);
        }
    };
    Array2D<PatchTrimCurve> rndPatchTrimCurve;

    // UV Index renderer for the surface
    class PatchBoundary : public /*BSplineRenderer::*/PatchRenderer {
    public:
        int side;
        void buildList() {
            parent->buildBoundaryPatch(dl,iPatch,jPatch,side,parent->matBoundary);
        }
    };
    Array2D<PatchBoundary> rndPatchBoundary[2];

    // UV Index renderer for the surface
    class PatchUVIndex : public /*BSplineRenderer::*/PatchRenderer {
    public:
        void buildList() {
            parent->buildUVIndexPatch(dl,iPatch,jPatch);
        }
    };
    Array2D<PatchUVIndex> rndPatchUVIndex;
    
    friend class Renderer;
    friend class PatchRenderer;
    friend class Wire;
    friend class PatchUVIndex;
    friend class PatchBoundary;
    friend class PatchTrimCurve;
    friend class PatchSeeThrough;
    friend class PatchFlat;
    friend class Control;
    friend class ControlIndex;
};


/**
 * An output window for communicating with the user
 */
class GLOutputWindow : public GLRenderer {
protected:
    list<string> lines;
    // string title;

    // Text color
    GLColor color;

    // Base position of y
    bool anchorTop;
    int anchorX,anchorY;

public:
    // Font to use for display
    GLFont *font;

    // Initialize
    GLOutputWindow(GLFont *font,GLColor color);

    void setAnchor(bool top,int x,int y);

    void clear();
    void addLine(const char *fmt,...);
    
    void onDraw();
};

// This is a command icon for selection mode
class ModeIcon : public GLCommandIcon {
private:
    string mode;
    static list<ModeIcon *> allIcons;
public:
    ModeIcon(string img1,string img2,string modeName,int icIndex);
    virtual ~ModeIcon();

    void onStateChange();
};

/**
 * A mode handler is used to run the system in a specific mode.  The tool
 * handler recieves system events and has a display routine that is used to
 * give feedback to the user
 */
class ModeHandler : public GLEventListener
{
public:
    virtual void start() {};
    virtual void stop() {};
};

/**
 * A general purpose point movement mode handler
 */
class PointMovementMH : public ModeHandler, public GLRenderer {
protected:
    // Selected medial point
    MedialPoint mp;

    // Normalized frame of the selected point
    SMLVec3f XU,XV;

    // Starting point for tracking motion
    SMLVec2f motionStart;

    // Projection, Model and Viewport for motion tracking
    double MP[16],MM[16];
    int MV[4];

    // Have we clicked the right button
    bool clicked,showInstructions;

    // Different modes
    enum Modes {NORMAL,TANGENT,RADIUS} mode;

    // Data about the transformed point, its original position
    SMLVec3f T,T0;

    // For radius manipulation
    float Rnew;

    // Text output window
    GLOutputWindow rndTextOut;

    // Draw method
    virtual void drawSelectedPoint();
    virtual void drawText() = 0;

    virtual void testMousePosition(int x,int y) = 0;
    virtual void computeFrame() = 0;
    virtual bool havePoint() = 0;
    virtual void applyMovement() = 0;

    virtual void doHandleMotion(int x,int y);

    // Schedule to test mouse position on the next draw event
    int mptX,mptY;
    bool doPassiveOnDraw,doMotionOnDraw;
    
    void scheduleMousePositionTest(int x,int y);
    void scheduleMouseMotion(int x,int y);

public:
    PointMovementMH();

    virtual void onDraw();

    virtual bool handleButton(int button,int state,int x,int y);
    virtual bool handlePassiveMotion(int x,int y);
    virtual bool handleMotion(int x,int y);
    virtual bool handleKeys(unsigned char key,int x,int y);

    virtual void start();
    virtual void stop();
private:
};

/**
 * A mode handler for surface point selection and manipulation
 */
class SurfaceSelectionMH : public PointMovementMH {
    // Data about the selected point
    int su,sv;

    void drawText();
    void computeFrame();
    void testMousePosition(int x,int y);
    void applyMovement();
    
    bool havePoint() {
        return su >= 0;
    }

    void onDraw();

public:
    SurfaceSelectionMH() {
        su = sv = -1;
    }

    // Allows coordinate system change
    virtual bool handleKeys(unsigned char key,int x,int y);
};

/**
 * A mode handler for surface point selection and manipulation
 */
class ControlPointEditMH : public PointMovementMH {
    // Index of the current control point
    int ci,cj;

    enum Frame {POLYGON=0,LOCAL=1,XYZ=2};

    // The coordinate frame index 
    Frame frameIndex;

    void drawText();
    void computeFrame();
    void applyMovement();
    void testMousePosition(int x,int y);
    
    bool havePoint() {
        return ci >= 0;
    }


public:
    ControlPointEditMH() {
        ci = cj = -1;
        frameIndex = POLYGON;
    }

    // Allows coordinate system change
    virtual bool handleKeys(unsigned char key,int x,int y);
};


// Optimization Mode
class OptimizationMH : public ModeHandler, GLRenderer {
private:
    // Optimization value
    double value,cpStartOptCost,optCost;

    // Current control point index
    int u,v;

    // Time left for current control point
    int cpTimeLeft;

    // Window size
    int window;
    int msPerWindow;

    // Optimization driver
    SplineOptDriver *od;
    SplineImageMatcher *match;

    // Sequence of control points left to optimize
    vector<int> cpList;

    // Make a sequence of control point indices
    void makeControlSequence();

public:
    OptimizationMH();
    void start();
    void stop();
    bool handleIdle();
    void onDraw();
};


class RigidMatchMode : public ModeHandler, GLRenderer {
private:
    SplineRigidProblem *srp;
    SplineImageMatcher *match;
    EvolutionaryStrategy *es;
    SolutionSpace *ss;

public:
    void start();
    void stop();
    bool handleIdle();
    void onDraw();
};
/*
class CurveAxisMatchMode : public ModeHandler, GLRenderer {
private:
    CurveAxisProblem *srp;
    EvolutionaryStrategy *es;
    SolutionSpace *ss;

public:
    void start();
    void stop();
    bool handleIdle();
    void onDraw();
};
*/
/**
 * A mode manager - manages mode handlers
 */
class ModeManager {
private:
    map<string,ModeHandler *> modes;
    ModeHandler *activeMode;
    string activeModeName;
public:
    ModeManager();
    void addMode(string name,ModeHandler *handler);
    void setMode(string name);
    string getActiveMode();
};



/**
 * An output only (status) window
 */
class GLFadingWindow : public GLOutputWindow {
private:
    // Time when text was last changed
    clock_t tLastChange,tView,tFade;
    bool isReset;

public:
    void onDraw();

    GLFadingWindow(GLFont *font,GLColor color);
    void clear();
};

/**
 * Global key binidng
 */
class GlobalKeyBinding {
public:
    int key;
    int state;
    string command;
    string desc;

    GlobalKeyBinding(int key,int state,string command);
};

/**
 * A global event handler
 */
class GlobalEventHandler : public GLEventListener {
private:
    list<GlobalKeyBinding> bind;
public:
    bool handleSpecial(int key,int x,int y);
    bool handleKeys(unsigned char key,int x,int y);
    void addKeyBinding(int key,int state,string command);
};

class MoveableWindow : public GLRenderer, public GLEventListener {
protected:
    // Current position of window
    SMLVec3f position, size;

    // Size of the 'hat'
    int szHat, szHatFont;

    // Drag position
    int dragX,dragY;

    // Window title
    string title;

    // Focus status
    bool hasFocus() {
        return false;
    }

    // Children's draw method
    virtual void drawClient(int width, int height) = 0;

public:
    void onDraw();

    virtual bool handleKeys(unsigned char key,int x,int y);
    virtual bool handleButton(int button,int state,int x,int y);
    virtual bool handleMotion(int x,int y);

    MoveableWindow(int x,int y,int w,int h);
};

/**
 * A command shell window, oh yeah!
 */
class CommandShell : public GLRenderer, public GLEventListener {
private:
    enum {REGHISTSIZE = 50};
    
    vector<string> lines;

    vector<string> history;

    string prompt;

    bool hasFocus,needRedraw;

    // Position and size of the command shell (use z coordinate)
    SMLVec3f pos,size;
    int rows,cols;

    // History index
    int idxHistory;

    // Drag position
    int dragX,dragY;

    // Cursor position
    int cursor;

    // An output window
    GLOutputWindow outwin;

public:
    void onDraw();

    bool handleKeys(unsigned char key,int x,int y);
    bool handleSpecial(int key,int x,int y);
    bool handleButton(int button,int state,int x,int y);
    bool handleMotion(int x,int y);

    void print(const char *text,int length = -1);

    CommandShell();
};


/****************************
 Imaging System - Allows for swapping of different image 
 types in the program
 ****************************************************/
class ImagingSystem {
public:
    // Types of images that the system can have loaded
    enum ITypes {NONE, DISTANCE, BINARY, GRAYSCALE};

    // Types of files that the system knows how to read
    enum IFiles {CUBE, RAW, BYU, GIPL, USE_EXTENSION};

    // Constructor
    ImagingSystem();
    ~ImagingSystem();

    // Load an image of one of the types
    void loadImage(const char *fname,ITypes type,IFiles fileType = CUBE,ostream &out=cout);

    // Convert the image type to another type
    void toBinary();
    void toGrayscale();
    void toDistance();

    // Get a reference to the image
    AbstractImage3D *getImage();

    // Get the image match object
    SplineImageMatcher *getImageMatch();

    // Reset the image/spline matcher object
    void resetImageMatch();

    // Get the image type
    ITypes getType() {
        return type;
    }

    // Save the image to a file
    void saveImage(const char *fname,IFiles type = CUBE);

private:
    // Current image type
    ITypes type;

    // Load different types of files
    void loadDistance(const char *fname,IFiles type = CUBE,ostream &out=cout);
    void loadBinary(const char *fname,IFiles type = USE_EXTENSION,ostream &out=cout);
    void loadGrayscale(const char *fname,IFiles type = USE_EXTENSION,ostream &out=cout);

    // Pointers to different image types
    BinaryImage *bim;
    DistanceTransform *dt;
    GrayImage *gray;

    // Pointers to different matchers
    SplineBinaryVolumeMatcher *mVolume;
    SplineDistanceMatcher *mDistance;
    // ProfileMatcher *mProfile;
    GradientImageMatcher *mGradient;
    
    // Discard the current image
    void discard();
};

// A background script processor for batch operation of the program
class ScriptProcessor : public GLEventListener {
public:
    bool handleIdle();

    void loadScript(const char *fname);

    ScriptProcessor() {
        tWaitEnd = 0;
    }

    void sleep(long ms) {
        tWaitEnd = ms + (1000 * clock()) / CLOCKS_PER_SEC;
    }
private:
    list<string> commands;
    clock_t tWaitEnd;
};

/*************************************************************
 Geodesic Renderer
 ************************************************************/
class GeodesicRenderer : public GLDisplayListRenderer {
public:
    void build();
    void setStart(float u,float v);
    void setEnd(float u,float v);
    bool needsStart();
    
    GeodesicRenderer();
private:
    

    bool hasStart,hasEnd;
    enum {pnum = 100};
    double uv[4];
};

// Use this method to execute a command
void executeCommand(string command);
void echoMessage(string message);

/***************************************************************
    Global Variables
 ***************************************************************/
// The spline that can be edited
extern MSpline *spline;

extern SplineDataCache *splineData;

// Spline renderer
extern BSplineRenderer *rndSpline;

// Mode manager
extern ModeManager *modeManager;

// Some defined fonts
extern GLFont *fntCourier18;
extern GLFont *fntCourier10;

// Command shell
extern CommandShell *commandShell;

// Status window
extern GLFadingWindow rndFadeStatus;

// Registry
extern Registry *settings;

// Script processor
extern ScriptProcessor scriptProcessor;

// Imaging subsystem
extern ImagingSystem imaging;

// Geodesic Renderer
extern GeodesicRenderer geodesicRenderer;


/**
 * This subroutine initializes the modes and mode icons
 */
void initModes();

// This is called to exit the program
void exitProcedure(int mode);

// Vertex paint routine
void glVertex(const SMLVec3f &v);

#endif

