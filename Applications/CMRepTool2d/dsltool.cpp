#pragma warning( disable : 4786 )  

#include "misc.h"
#include "tools.h"
#include "dslwin.h"
#include "dsltool.h"
#include "undo.h"
#include "ui.h"
#include "reg.h"
#include "uibind.h"

#include <Fl/fl_file_chooser.H>
#include <Fl/Fl_Color_Chooser.H>
#include <Fl/gl.h>
//#include <Gl/glu.h>

#include <registry.h>
#include "ispace.h"
#include "likehood.h"
#include "readpgm.h"
#include "SplineMrep.h"

// #include <cimage_ops.h>

#include <vector>

using namespace std;

class JunkStream : public ostringstream {
public:
	basic_ostream<char,char_traits<char> >& write(const char_type *s,streamsize n) {
		char *text=new char[n+1];
		cout << "YO! " << text << endl;
		strncpy(text,s,n);
		txtOutput->insert(text);
		delete text;
		return *this;
	}	

	basic_ostream<char,char_traits<char> >& put(char ch) {
		cout << "YOOOOOOOOOOOOOOOOOOOOOOOOOOOOO";
	}
};

JunkStream js;



/***********************************************************
* Global variables for the DSLTool application				  *
***********************************************************/
// The image that is being displayed
// cimage baseImage;
// Array2D<GLbyte> baseBytes;

// Scale space built on top of the image
// ScaleSpace *scaleSpace = NULL;
// CImageSpace *imageSpace = NULL;
ImageSpace *imageSpace = NULL;

// The image match computer
IMatchComputer iMatchComputer;

// Spline match object
SplineSampleImageMatch *iMatchSpline = NULL;

// The model that we will be working with
MGraph model;

// Number of selected primitives in current model
int nSelected;

// Application settings in registry format
Registry regOptions;

// We have an undo buffer
UndoBuffer<MGraph> *undoBuffer = new UndoBuffer<MGraph>(50);

// A spline MRep
SplineObject *spline;
RegularSplineSample *splineSample;

// Which mode are we running in
DisplayModeEnum eDisplayMode = DISPLAY_MREP;

/**********************************************************************
* Customized undo methods
**********************************************************************/
void updateUndoMenus() {
   if(undoBuffer->canUndo()) {
		miUndo->activate();
      bUndo->activate();
   }
   else {
		miUndo->deactivate();
      bUndo->deactivate();
   }
	
   if(undoBuffer->canRedo()) {
		miRedo->activate();
      bRedo->activate();
   }
   else {
      miRedo->deactivate();
      bRedo->deactivate();
   }
}

void pushModel() {
   // Copy the model
   // MGraph copy = model;
	
   // Save the model in the stack
   undoBuffer->push(model);
	
   // Update the undo/redo buttons/menus
   updateUndoMenus();
}

/**********************************************************************
* In-Responce Methods
**********************************************************************/
void onModelOrSelectionChange() {
   if(eDisplayMode == DISPLAY_MREP || eDisplayMode == DISPLAY_PEM) {
		
		// Update the selection box
		computeSelectionBox();
		
		// Update the editor
		updatePrimitiveEditor(&model);
		
		// Recompute the medialness
		if(imageSpace) {
			double medness = iMatchComputer.compute(*model.selection(),*imageSpace);
			outMedialness->value(medness);
		}
	}
	else if(eDisplayMode == DISPLAY_SPLINE) {
		
		// Check the timestamp of the spline
		if(splineSample->tsOverall < spline->tsStructure) {
			delete splineSample;
			splineSample = new ArcSplineSample(spline,20);
			
			if(imageSpace) {
				if(iMatchSpline) delete iMatchSpline;
				iMatchSpline = new SplineSampleImageMatch(splineSample,imageSpace);
			}
		}
		else {
			splineSample->update();
		}

		// Recompute the medialness
		if(imageSpace) {
			double medness = iMatchSpline->integrateImageMatch();
			outMedialness->value(medness);

			// Compute penalties and stuff
			SNGSegment sng(spline,0,0,spline->curves.front()->size()-1);
			SplineMatchProblem smp(&sng,imageSpace);
			cout << "reg prior = " << smp.computeRegularizationPrior() << endl;
		}
	}
}

/*************************************************************
FL and GL Initialization code
***********************************************************/
void initGL() {
	glEnable(GL_NORMALIZE); 
		
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	 
	if(regOptions.getBooleanValue("display.gl.blend",true)) {
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	else {
		glDisable(GL_BLEND);
	}
	
	if(regOptions.getBooleanValue("display.gl.smooth",true)) {
		glEnable(GL_LINE_SMOOTH);
      // glEnable(GL_POLYGON_SMOOTH);
		glLineWidth(regOptions.getDoubleValue("display.gl.lineWidth",2.0));	
	}
	else {
		glDisable(GL_LINE_SMOOTH);
      // glDisable(GL_POLYGON_SMOOTH);
		glLineWidth(1.0);
	}
	
	glFlush();
	
}	

void status(const char *message) {
	cout << message << endl;
}

string modelFileName;

/*
* Routine:		histEq
*              A homegrown image enchancement routine
* Parms:
*  input		Image as array of doubles
*  output		Output image as array of bytes
*  size		Number of pixels in both images
*  histSize	Number of bins in by histogram (try 8 or 16)
*/
void histEq(cimage &input,Array2D<GLbyte> &output,int histSize) {
	double max = 0;
	const int minPix = 0;
   
	
   int w = cimage_xdim(input);
   int h = cimage_ydim(input);
   // Array2D<unsigned char> output(w,h);
	
   for(int i=0;i<w;i++) {
      for(int j=0;j<h;j++) {
         double pix = cimage_get(input,i,j);
         max = (pix>max) ? pix : max;
      }
   }
	
	max *= 1.0001;
	
	int *hist = new int[histSize];
	memset(hist,0,sizeof(int)*histSize);
	
   for(i=0;i<w;i++) {
      for(int j=0;j<h;j++) {
			int bin = (int)((histSize*cimage_get(input,i,j)) / max);
			hist[bin] ++;
      }
	}
	
	int freeCells = 256 - histSize*minPix;
	int *startHist = new int[histSize+1];
	startHist[0]=0;
	
	for(i=0;i<histSize;i++) {
		hist[i] = minPix + (hist[i]*freeCells) / (w*h);
		startHist[i+1] = startHist[i]+hist[i];
	}
	
	for(i=0;i<w;i++) {
      for(int j=0;j<h;j++) {
         double pix = cimage_get(input,i,j);
         int bin = (int)((histSize*pix) / max);
			double frac = (pix*histSize)/max - bin;
			output(i,j) = startHist[bin] + (int)(frac*hist[bin]);
      }
	}
	
	delete hist;
	delete startHist;
}


void imageChanged(cimage newImage) {
	// cimage_destroy(baseImage);
	cimage baseImage = newImage;
   
   int w = cimage_xdim(newImage);
   int h = cimage_ydim(newImage);
	
   cimage imGray = cimage_init(); // cimage_new(CGREY,w,h,1);
   cimage imAHE = cimage_init();
	
   cimage_convert(newImage,imGray,CGREY);
   
   int nreg[] = {1,8,8};
   //cimage_ahe(imGray, imAHE, nreg, 5.0);
	//baseImage = imAHE;
	baseImage = imGray;
   
   cimage_max_min(baseImage);
	
	// Allocate a byte array for GL display
	// baseBytes.resize(cimage_xdim(baseImage),cimage_ydim(baseImage));
   
	double iMin = cimage_imin(baseImage);
	double iRange = cimage_imax(baseImage) - iMin;
	
	/*
	for(int j=0;j<cimage_ydim(baseImage);j++) {
		for(int i=0;i<cimage_xdim(baseImage);i++) {
			baseBytes(i,j) = (unsigned char) (255.0 * (cimage_get(baseImage,i,j)-iMin) / iRange);
		}
	}*/
   
   //histEq(baseImage,baseBytes,16);
	
	// Create a scale space
	if(imageSpace)
		delete imageSpace;
	
   // Create an image space
   imageSpace = new ScaleSpace(baseImage);         

   iMatchSpline = new SplineSampleImageMatch(splineSample,imageSpace);
   //imageSpace = new CImageSpace(baseImage);         
	
   // A little test
   /*
   double b1 = imageSpace->getValue(Vector2D(106,60),4);
   double b2 = imageSpace->getValue(Vector2D(106,60),8);
   double b3 = imageSpace->getValue(Vector2D(106,60),16);
	
	  double g1 = imageSpace->getGradient(Vector2D(106,60),4).twoNorm();
	  double g2 = imageSpace->getGradient(Vector2D(106,60),8).twoNorm();
	  double g3 = imageSpace->getGradient(Vector2D(106,60),16).twoNorm();
   */
	
	// Unless we are in slow mode of display, update the zoom image
	winGL->updateTexture();
	winGL->redraw();
   
   // Compute the first few slices
   winImageProgress->hotspot(inILoadProgress);
   winImageProgress->show();
   double scale = 1.0;
   for(int i=0;i<10;i++) {
      char label[40];
      sprintf(label," s = %lg",scale);
      inILoadProgress->value(i/9.0);
      inILoadProgress->label(label);
      winImageProgress->redraw();
		
      if (Fl::ready()) 
         Fl::check();      
		
      imageSpace->getValue(Vector2D(0,0,true),scale);
      scale *= sqrt(2.0);
		
   }
   winImageProgress->hide();
}
void readImage(char *fname,bool forcePGM=false,bool forceIM=false) {
	cimage temp = cimage_init();
	bool readOK = false;
	
	// Determine the type of the image
	if(strstr(fname,".pgm") || forcePGM) {
		readOK = (PGM_OK == readpgm(fname,temp));
	}
	else if(strstr(fname,".im") != NULL || strstr(fname,".Im") != NULL || forceIM) {
		readOK = 0 != cimage_readusrimage(temp, fname);
	}
	else {
		status("Unknown image extention!");
		return;
	}
	
	if(readOK) {
		// cimage t2 = cimage_init();
		// cimage_multconst(temp,t2,0.25);		
		imageChanged(temp);
		status("Successfully read image file");
	}
	else { 
		cimage_destroy(temp);
		status("Unable to read image file");
	}
}

void readModel(char *fname,char *tagName) {
	try {
		Registry folder(fname);
		modelFileName = fname;
		
		// Generate tag name
		string tag = "model";
		if(tagName) {
			tag += "[";
			tag += tagName;
			tag += "]";
		}
		
		// Save model on undo stack before making changes
		pushModel();
		
		// Load appropriate model
		if(!folder.hasKey((char *)tag.c_str()))
			return;
		
		model.loadFromModelFile(folder.getSubFolder((char *)tag.c_str()));
		
		// Scale the model if it is too big to fit on the 1 by 1 image
		Vector2D ex1,ex2;
		model.getExtents(ex1,ex2);
		
		// Infinity norm is the abs. value of the greater of two components of the vector
		double extent = (ex1-ex2).infinityNorm(); 
		if(extent > 1.5) {
			model.scaleToFit(0.1,0.1,0.8,0.8);
		}
		
		winGL->haveSelection = false;				
		onModelOrSelectionChange();
		winGL->redraw();
		
      // Switch to selection tool
		model.select();
      // uiRunTool(TOOL_SELECT);
		
		status("Successfully read model file");      
	}
	catch(RException exc) {
		cerr << exc;
		status("Unable to read model file");
	}
}


void uiOpenPGM(Fl_Menu_*, void*) {
	const char *flast = regOptions.getStringValue("history.pgmFile","");
	char *fname = fl_file_chooser("Select an Image File", "PGM Files (*.im)", flast);
	if(fname) {
		uiScript.uiCommand("set history.pgmFile \"%s\"",fname);
		uiScript.uiCommand("load pgm \"%s\"",fname);
	}
}

void uiOpenImage(Fl_Menu_*, void*) {
	const char *flast = regOptions.getStringValue("history.imFile","");
	char *fname = fl_file_chooser("Select an Image File", "UNC image Files (*.im)", flast);
	if(fname)  {
		uiScript.uiCommand("set history.imFile \"%s\"",fname);
		uiScript.uiCommand("load image \"%s\"",fname);
	}
}

void uiOpenModel(Fl_Menu_*, void*) {
	const char *flast = regOptions.getStringValue("history.modelFile","");
	char *fname = fl_file_chooser("Select an Model File", "*.mod*",flast);
	if(fname) {
		uiScript.uiCommand("set history.modelFile \"%s\"",fname);
		uiScript.uiCommand("load model \"%s\"",fname);
	}
}

void makeNewModel(int nFigs,int nAtoms,bool keepOld,double rho) {
   // Make sure nothing is selected
	model.deselect();
	
	// Empty the model if needed
   if(!keepOld) {
      model.removeAll();
   }
	
   // Add all figures as needed
	model.beginTopologyEdit();
	for(int i=0;i<nFigs;i++) {
		vector<MAtom> atoms;
		atoms.reserve(nAtoms);		
		for(int j=0;j<nAtoms;j++) {
			MAtom atom(1+3*i,1+3*j,1,-90,90);
			atom.select();
			if(j==0 || j==nAtoms-1)
				atom.endRC(1.2);
			atoms.push_back(atom);
		}
		model.addFigure(atoms);
	}
	model.endTopologyEdit();
	
	// Get the added (selected) atoms and place them in the middle of the current view
	MNodeGroup *ng = model.selection();
	Vector2D wx(winGL->winPosition);
	Vector2D ws(winGL->winSize,winGL->winSize);
	ng->scaleToFit(wx.x,wx.y,ws.x,ws.y);
	ng->translateBy(wx+ws*0.5-ng->getCOG());
	ng->scaleBy(0.6);
	ng->rho(rho);
	
   // Selection has changed
	onModelOrSelectionChange();
	
   status("Successfully created new model");
}

void uiNewModel(int nFigs,int nPrims,bool keepOld,double rho) {
   pushModel();
   makeNewModel(nFigs,nPrims,keepOld,rho);
   winGL->redraw();
}

void uiDeleteSelected(Fl_Menu_*, void*) {
   // This is simple enough
	model.beginTopologyEdit();
	model.removeGroup(*model.selection());
	model.endTopologyEdit();
	
   // Switch to selection tool
   uiRunTool(TOOL_SELECT);
	
   onModelOrSelectionChange();
   winGL->redraw();
}

// Array of tool buttons
vector<Fl_Button *> vToolButtons;

// Array of tool pointers
vector<ToolHandler *>vTools;

extern SplineToolHandler *uiCallbackSplineTool;

// This method runs one of the selected tools
void uiRunTool(ToolEnum id) {
	eToolId = id;

	// Create a new array of tool hanlerrs
	if(vTools.size() == 0) {
		vTools.push_back(new ZoomPanToolHandler(*winGL,true));
		vTools.push_back(new ZoomPanToolHandler(*winGL,false));
		vTools.push_back(new SelectionToolHandler(*winGL));
		vTools.push_back(new TranslateToolHandler(*winGL));
		vTools.push_back(new ScaleRotateToolHandler(*winGL,false));
		vTools.push_back(new ScaleRotateToolHandler(*winGL,true));
		vTools.push_back(new StretchToolHandler(*winGL));
		vTools.push_back(new BendToolHandler(*winGL));
		
      uiCallbackSplineTool = new SplineToolHandler(*winGL);
      vTools.push_back(uiCallbackSplineTool);

		vTools.push_back(new PEMToolHandler(*winGL));
		vTools.push_back(new BranchToolHandler(*winGL));

		// Array of tool buttons
		vToolButtons.push_back(bZoom);
		vToolButtons.push_back(bPan);
		vToolButtons.push_back(bSelect);
		vToolButtons.push_back(bTranslate);
		vToolButtons.push_back(bScale);
		vToolButtons.push_back(bRotate);
		vToolButtons.push_back(bStretch);
		vToolButtons.push_back(bBend);
		vToolButtons.push_back(bSpline);
		vToolButtons.push_back(bPEM);
		vToolButtons.push_back(bBranch);
	}

	// Get the current tool
	ToolHandler *tool = vTools.at(id);

	// Set the current tool button to enabled
	for(unsigned int iTool=0;iTool<vTools.size();iTool++) {
		vToolButtons.at(iTool)->value((iTool == id) ? 1 : 0);
	}

	// Set the context help and the label
	lToolName->label(tool->getLabel());
	lContextHelp->label(tool->getContextHelp());

	// Redraw the main window
	winGL->setTool(vTools[id]);
	winMain->redraw();

	// Other, custom tasks
	switch(id) {
	case TOOL_PEM : 
		// Show the primitive editor
		uiScript.uiCommand("show winPrimitive");
		break;
	};

}

ToolEnum eToolId;

void uiSetZoom(double factor) {
	winGL->winSize = factor;
	winGL->winPosition = Vector2D(0,0);
	
	// winGL->updateImageDL();
	winGL->redraw();
}


double modelThinLastAmount;
void uiStretchMode(Fl_Button *b, void*) {
	pushModel();
	modelThinLastAmount = 0;
	winModelThinner->show();
};

void uiModelThin(double amount) {
	// Thin the model
	// model.selection()->thin(pow(2.0,));
	model.selection()->thin(pow(2.0,amount-modelThinLastAmount));
	modelThinLastAmount = amount;
	
   // Update displays
   onModelOrSelectionChange();
	
	// Repaint
	winGL->redraw();
}

void uiUndo(Fl_Menu_*, void*) {
   // Push model on the undo stack
   model = undoBuffer->undo(model);
	
   // Update displays
   onModelOrSelectionChange();
	
   // Redo the buttons
   updateUndoMenus();
   
   // Repaint
	winGL->redraw();
}

void uiRedo(Fl_Menu_*, void*) {
   // Push model on the undo stack
   model = undoBuffer->redo(model);
	
   // Update displays
   onModelOrSelectionChange();
	
   // Redo the buttons
   updateUndoMenus();
   
   // Repaint
	winGL->redraw();
}

// This method updates all the fields in the primitive editor.
void updatePE(MNode *node) {
   // Enable all controls and set values
   
   /*
   inPWinX->activate();
   inPWinY->activate();
   inPWinR->activate();
   inPWinTA->activate();
   inPWinTO->activate();
   inPWinLeftPol->activate();
   inPWinRightPol->activate();
   inPWinEndcap->activate();
   sliderPWinEndcap->activate();
   inPWinLeftX->activate();
   inPWinLeftY->activate();
   inPWinLeftN->activate();
   inPWinLeftSigma->activate();
   inPWinLeftBNess->activate();
   inPWinRightX->activate();
   inPWinRightY->activate();
   inPWinRightN->activate();
   inPWinRightSigma->activate();
   inPWinRightBNess->activate();   
   winPrimitive->activate();
   */
	
	BAtom lba = node->bAtom(MAtom::LEFT);
	BAtom rba = node->bAtom(MAtom::RIGHT);
	
   inPWinX->value(node->x().x);
   inPWinY->value(node->x().y);
   inPWinR->value(node->r());
   inPWinTA->value(node->faDeg());
   inPWinTO->value(node->oaDeg());
	
   inPWinLeftPol->value(lba.polarity);
   inPWinRightPol->value(rba.polarity);
   inPWinEndPol->value(node->polarity(MAtom::HEAD));
	
   inPWinEndcap->value(node->endRC());
   sliderPWinEndcap->value(node->endRC());
	
   // Also the boundary sites
   inPWinLeftX->value(lba.x.x);
   inPWinLeftY->value(lba.x.y);
   inPWinLeftN->value(lba.slope);
   inPWinLeftSigma->value(lba.scale);   
	
   inPWinRightX->value(rba.x.x);
   inPWinRightY->value(rba.x.y);
   inPWinRightN->value(rba.slope);
   inPWinRightSigma->value(rba.scale);   
	
   // Compute the boundary strength
   if(imageSpace) {		
      double bLeft = iMatchComputer.compute(lba,*imageSpace);
      double bRight = iMatchComputer.compute(rba,*imageSpace);
      inPWinRightBNess->value(bRight);
      inPWinLeftBNess->value(bLeft);
   }
	
   // Set the profiles
   winRightProfile->set(rba);
   winLeftProfile->set(lba);
   winRightGradient->set(rba);
   winLeftGradient->set(lba);
}

void updatePrimitiveEditor(MNode *node) {
   static char str[256];
   sprintf(str,"One primitive selected");
   gPWinSelected->label(str);
   updatePE(node);
   winPrimitive->redraw();
}

// This method updates all the fields in the primitive editor.
void updatePrimitiveEditor(MGraph *model) {
   MNodeGroup *sg = model->selection();
	
	// Print number of selected prims
   static char str[256];
   if(sg->size()==0) 
      sprintf(str,"Zero primitives selected");
   else if(sg->size()==1) 
      sprintf(str,"One primitive selected");
   else
      sprintf(str,"%d primitives selected",sg->size());
   
   gPWinSelected->label(str);  
   
   if(sg->size() != 1) {
      // Disable all controls
      /*
      inPWinX->deactivate();
      inPWinY->deactivate();
      inPWinR->deactivate();
      inPWinTA->deactivate();
      inPWinTO->deactivate();
      inPWinLeftPol->deactivate();
      inPWinRightPol->deactivate();
      inPWinEndcap->deactivate();
      sliderPWinEndcap->deactivate();
      */
		
      inPWinX->value(0);
      inPWinY->value(0);
      inPWinR->value(0);
      inPWinTA->value(0);
      inPWinTO->value(0);
      inPWinLeftPol->value(0);
      inPWinRightPol->value(0);
      inPWinEndPol->value(0);
      inPWinEndcap->value(0);
      sliderPWinEndcap->value(0);
      inPWinLeftX->value(0);
      inPWinLeftY->value(0);
      inPWinLeftN->value(0);
      inPWinLeftSigma->value(0);
      inPWinLeftBNess->value(0);
      inPWinRightX->value(0);
      inPWinRightY->value(0);
      inPWinRightN->value(0);
      inPWinRightSigma->value(0);
      inPWinRightBNess->value(0);
      
      /*
      inPWinLeftX->deactivate();
      inPWinLeftY->deactivate();
      inPWinLeftN->deactivate();
      inPWinLeftSigma->deactivate();
      inPWinLeftBNess->deactivate();
		
		  inPWinRightX->deactivate();
		  inPWinRightY->deactivate();
		  inPWinRightN->deactivate();
		  inPWinRightSigma->deactivate();
		  inPWinRightBNess->deactivate();
      */
		
      // Clear the intensity profile displays
      winLeftProfile->unset();
      winRightProfile->unset();
      winLeftGradient->unset();
      winRightGradient->unset();
      winPrimitive->redraw();
		
      // winPrimitive->deactivate();
      
		
      return;
   }
   else {
		updatePE(sg->node(0));
   }
	
   winPrimitive->redraw();
}

// Set value for all selected atoms
void uiPEditSet(int which) {
	// Get the selection
	MNodeGroup *sg = model.selection();
	
   // Store the model for undo, unless we had just called this method
	static MGraph *lastModel = NULL;
	if(lastModel == NULL || lastModel != undoBuffer->getTail()) {
		pushModel();	
		lastModel = undoBuffer->getTail();
	}
	
	
   // Set the value on all selected primitives
	for(int i=0;i<sg->size();i++) {
		MNode *node = sg->node(i);
		switch(which) {
		case 0 : node->x(Vector2D(inPWinX->value(),node->x().y));break;
		case 1 : node->x(Vector2D(node->x().x,inPWinY->value()));break;
		case 2 : node->r(inPWinR->value());break;
		case 3 : node->faDeg(inPWinTA->value());break;
		case 4 : node->oaDeg(inPWinTO->value());break;
			
			// Set endcap using input
		case 5 : node->endRC(inPWinEndcap->value());
			sliderPWinEndcap->value(inPWinEndcap->value());
			break;
			
			// Set endcap using slider
		case 6 : node->endRC(sliderPWinEndcap->value());
			inPWinEndcap->value(sliderPWinEndcap->value());
			break;
			
		case 7 : node->polarity(inPWinLeftPol->value(),MAtom::LEFT);break;
		case 8 : node->polarity(inPWinRightPol->value(),MAtom::RIGHT);break;
      case 9 : node->polarity(inPWinEndPol->value(),MAtom::HEAD);
			node->polarity(inPWinEndPol->value(),MAtom::TAIL);
			break;
		}
	}
	
   // Refill the window
   onModelOrSelectionChange();
	
   // Redraw
   winGL->redraw();
}


void uiRegister(Fl_Button *b, void*) {
   runRegistration();
   bRegister->value(0);
   winMain->redraw();
}

void uiAssist(Fl_Button *b, void*) {
   // runNodeOptimization();
	if(eDisplayMode==DISPLAY_MREP || eDisplayMode==DISPLAY_PEM) {
		runMultiNodeOpt();
	}
	else if(eDisplayMode==DISPLAY_SPLINE) {
		runSplineMatch(spline,imageSpace);
	}
	
   bAA->value(0);
   winMain->redraw();
}

void uiShowInterpolator(Fl_Menu_ *,void *) {
	MFigure *f = getSingleSelectedFigure();	
	onModelOrSelectionChange();
	winGL->redraw();
	
	winInterpolate->show();
}

void uiSSUpdateMain(Fl_Button *b, void*) {
   // This causes the main window display to be replaced with the image currently in scale space.
}

/******************************************************
* Command processor for scripted commands
******************************************************/
bool scompare(const char *src,const char *dst) {
   int i=0;
   
   if(src==NULL) {
      if(dst == NULL) 
         return true;
      return false;
   }
   else if(dst == NULL) 
      return false;
	
   while(1) {
      if(tolower(src[i]) != tolower(dst[i])) return false;
      if(src[i]==0) return true;
      i++;
   }
}

// Is second string prefix of the first?
bool sprefix(const char *text,const char *prefix) {
   int i=0;
   
   if(prefix == NULL) 
      return true;
	
	if(text == NULL)
		return false;
	
   while(1) {
      if(prefix[i] == 0) return true;
		if(tolower(text[i]) != tolower(prefix[i])) return false;      
      i++;
   }
}

void readSpline(const char *file) {
	try {
		Registry folder;
		folder.readFromFile(file);
		// TODO: spline = new SplineMRep(&folder);
		spline->load(folder);
		
		onModelOrSelectionChange();
		winGL->redraw();
	} catch (...) {
		status("Error loading spline file");
	}
}

void saveSpline(const char *file) {
	try {
		Registry folder;
		// spline->save(&folder);
		spline->save(folder);
		folder.writeToFile(file);
		
		status("Wrote spline to file");


	} catch (...) {
		status("Error writing spline file");
	}
}

void newSpline(int size,double rho) {	
	spline = new SplineObject();
	spline->addCurve(size,rho);

	onModelOrSelectionChange();
	winNewSpline->hide();
	winGL->redraw();
}

void writeSplineAsTable(const char *fname) {
	FILE *f = fopen(fname,"wt");
	if(spline->curves.size()==1) {
		SplineCurve *crv = spline->curves[0];
		for(int i=0;i<crv->size();i++) {
			MySMLVec3f v = crv->getControlVec(i);
			fprintf(f,"%lg %lg %lg\n",v.x,v.y,v.z);
		}
	}
	fclose(f);
}

void commandProcessor(int argc,char **argv,bool fromScript) {
   string command(argv[0]);
   for(int i=1;i<argc;i++) {
      command += " ";
      command += argv[i];
   }
   outScript->add(command.c_str());
	
   if(scompare(argv[0],"set")) {
      // A set command.  Every set command sets an option in the registry.
      regOptions.setStringValue(argv[1],argv[2]);
		
		// If the command is from script, set the widget as well
		if(fromScript)
			uiBind.setWidget(argv[1],argv[2]);
		
      // Redraw the display
		if(sprefix(argv[1],"display.gl")) {
			// Window must revalidate itself
			winGL->valid(false);
			winGL->updateTexture();
         winGL->redraw();
      }
		
      else if(sprefix(argv[1],"imageMatch")) {
			// Image computer can not afford to check the registry every time we call it.  Instead we set all the
			// values on the image computer
			iMatchComputer.readRegistry(regOptions.getSubFolder("imageMatch"));
      }
   }
	
   else if(scompare(argv[0],"show")) {
      // Second argument is window name
      Fl_Window *win = NULL;
      int x = (argc > 2) ? atoi(argv[2]) : -1;
      int y = (argc > 3) ? atoi(argv[3]) : -1;
		
      if(scompare(argv[1],"winMain")) {
         winMain->show();
         winGL->show();
         if (x >= 0 && y >= 0) 
            winMain->position(x,y);
      }
		
      else if(scompare(argv[1],"winScaleSpace")) {
         winScaleSpace->show();
         boxScaleSpace->show();
         if (x >= 0 && y >= 0) 
            winScaleSpace->position(x,y);
      }
		
      else if(scompare(argv[1],"winPrimitive")) {
         winPrimitive->show();
         winLeftProfile->show();
         winRightProfile->show();
         winLeftGradient->show();
         winRightGradient->show();
         if (x >= 0 && y >= 0) 
            winPrimitive->position(x,y);
      }
		
      else if(scompare(argv[1],"winNewModel")) {
         winNewModel->show();
         if (x >= 0 && y >= 0) 
            winNewModel->position(x,y);
      }
		
      else if(scompare(argv[1],"winCoreWizard")) {
         winCoreWizard->show();
         if (x >= 0 && y >= 0) 
            winCoreWizard->position(x,y);
      }
   }
	
   else if(scompare(argv[0],"print")) {
      if(argc==1) {
			list<string> keys;
         regOptions.collectAllKeys(keys);
			
			for(list<string>::iterator i=keys.begin();i!=keys.end();i++) {
				char *key = (char *) i->c_str();
				outScript->add(key,(void *)regOptions.getStringValue(key,NULL));
			}
      }
      else {
         outScript->add(regOptions.getStringValue(argv[1],"Unknown key"));
      }
   }
	
	else if(scompare(argv[0],"load")) {
		// Load a model file
		if(scompare(argv[1],"image")) {
			readImage(argv[2],false,true);
		}
		else if(scompare(argv[1],"pgm")) {
			readImage(argv[2],true,false);
		}
		else if(scompare(argv[1],"model")) {
			readModel(argv[2],(argc>3)?argv[3]:NULL);
		}
		else if(scompare(argv[1],"spline")) {
			readSpline(argv[2]);
		}
	}
	
	else if(scompare(argv[0],"save")) {
		// Load a model file
		if(scompare(argv[1],"model")) {
			// Load the model file
			Registry folder;
			try {
				// Read a model if it's there
				folder.readFromFile(argv[2]);
			} catch(...) {}
			
			try {
				// Generate tag name
				string tag = "model";
				if(argc>3) {
					tag += "[";
					tag += argv[3];
					tag += "]";
				}
				
				model.writeToModelFile(folder.getSubFolder((char *)tag.c_str()));
				folder.writeToFile(argv[2]); 
				
				status("Model file was saved");
			}
			catch(RException exc) {
				cout << exc;
				status("Failed to save model file");
			}
		}
		else if(scompare(argv[1],"matlab")) {
			try {
				model.writeToMatlabFile("argv[2]");
			}
			catch(...) {
				status("Failed to save model file");
			}
		}
		else if(scompare(argv[1],"xml")) {
			try {
				model.writeXML(argv[2]);
			}
			catch(...) {
				status("Failed to save XML file");
			}
		}
		else if(scompare(argv[1],"spline")) {
			if(strstr(argv[2],".dat")) {
				writeSplineAsTable(argv[2]);
			}
			else {
				saveSpline(argv[2]);
			}
		}
	}
	
	else if(scompare(argv[0],"resample")) {
		if(scompare(argv[1],"uniform")) {
			pushModel();
			
			MFigure *f = getSingleSelectedFigure();
			resampleFigureUniform(f);
			
			winInterpolate->hide();
			onModelOrSelectionChange();
			winGL->redraw();
		}
	}
	
	else if(scompare(argv[0],"interpolate")) {
		int nAtoms = 2;
		double rho = 0.25;
		try {
			if(argc > 1)
				nAtoms = atoi(argv[1])+1;
			if(argc > 2)
				rho = atof(argv[2]);
		} catch(...) {
		}
		
		pushModel();
		
		MFigure *f = getSingleSelectedFigure();
		interpolateFigure(f,nAtoms);
		model.selection()->rho(rho);
		
		winInterpolate->hide();
		onModelOrSelectionChange();
		winGL->redraw();
	}
	
	else if(scompare(argv[0],"register")) {
		uiRegister(NULL,NULL);
	}
	
	else if(scompare(argv[0],"optimize")) {
		uiAssist(NULL,NULL);
	}
	
	else if(scompare(argv[0],"execute")) {
		UserInterfaceScript uiTempScript(commandProcessor);
		uiTempScript.load(argv[1]);
		uiTempScript.run(argc-2,&(argv[2]));
	}

	else if(scompare(argv[0],"create") && argc > 1) {
		if(scompare(argv[1],"spline") && argc > 3) {
			try {
				int n = atoi(argv[2]);
				double rho = atof(argv[3]);
				newSpline(n,rho);
			} catch(...) {
				status("Error parsing command!");
			}
		}
	
	}

	else if(scompare(argv[0],"minimize")) {
		winMain->iconize();
		js << "Yo yo yo!" << endl;
	}

	else if(scompare(argv[0],"rho") &&  argc > 1) {
		double rho = atof(argv[1]);
		for(unsigned int iCurve=0;iCurve<spline->curves.size();i++) {
			spline->curves[iCurve]->setRho(rho);
		}

	}

	else if(scompare(argv[0],"mode") &&  argc > 1) {
		int id = atoi(argv[1]);
		uiRunTool((ToolEnum)id);
	}
	else if(scompare(argv[0],"quit")) {
		regOptions.writeToFile("dsltool.ini"); 
		exit(0);
	}
	
	outScript->bottomline(outScript->size());
}

void uiSaveModel(Fl_Menu_*, void*) {
	if(modelFileName.length()==0)
		uiSaveModelAs(NULL,NULL);
	else 
		uiScript.uiCommand("save model %s",modelFileName.c_str());
}

void uiSaveModelXML(Fl_Menu_*, void*) {
	char *fname = fl_file_chooser("Select an Model File", "*.xml",NULL);
	if(fname) {
		uiScript.uiCommand("save xml \"%s\"",fname);
	}
}

void uiSaveModelAs(Fl_Menu_*, void*) {
	const char *flast = regOptions.getStringValue("history.modelFile","");
	char *fname = fl_file_chooser("Select an Model File", "*.model",flast);
	if(fname) {
		uiScript.uiCommand("set history.modelFile \"%s\"",fname);
		modelFileName = fname;
		uiSaveModel(NULL,NULL);
	}
	else {
		status("Model file was not specified");
	}
}

void uiOpenSpline(Fl_Menu_*, void*) {
	const char *flast = regOptions.getStringValue("history.splineFile","");
	char *fname = fl_file_chooser("Select a Spline File", "*.bspline",flast);
	
	if(fname) {
		uiScript.uiCommand("set history.splineFile \"%s\"",fname);
		uiScript.uiCommand("load spline \"%s\"",fname);
	}
	else {
		status("Model file was not specified");
	}
}


void uiSaveSplineAs(Fl_Menu_*, void*) {
	const char *flast = regOptions.getStringValue("history.splineFile","");
	char *fname = fl_file_chooser("Select a Spline File", "*.bspline",flast);
	
	if(fname) {
		uiScript.uiCommand("set history.splineFile \"%s\"",fname);
		uiScript.uiCommand("save spline \"%s\"",fname);
	}
	else {
		status("Model file was not specified");
	}
}

void uiNewSplineCreate(Fl_Return_Button*, void*) {
	uiScript.uiCommand("create spline %d %lf",(int)inNewSplineControls->value(),inNewSplineRho->value());
}

UserInterfaceScript uiScript(commandProcessor);
UIBinder uiBind;

void uiSetCmd(Fl_Widget *widget) {
	string command;
	string key,value;
	
	if(uiBind.readWidget(widget,key,value)) {
		command = "set " + key + " " + value;
		uiScript.uiSliderCommand(widget,command.c_str());
	}
}

void uiInputCommand(Fl_Input*, void*) {
   uiScript.runCommand(inScript->value());
}

void uiColorButton(Fl_Light_Button *,void *) {
}

bool uiPopulateTreeWidgetRec(Fl_ToggleTree *tree,Registry &reg,string prefix) {
	char **p;
	bool needsub = true;

	// Get all the folder keys
	int nKeys = reg.getKeyArraySize();
	char **keys = new char *[nKeys];
	reg.getValueKeys(keys);
	for(p=keys;*p!=NULL;p++) {
		string *fullKey = new string(prefix+*p);
		Fl_ToggleNode *node = new Fl_ToggleNode(*p,0,NULL,fullKey);
		if(needsub) {
			needsub = false;
			tree->add_sub(node);
		}
		else {
			tree->add_next(node);
		}
	}
	reg.getFolderKeys(keys);
	for(p=keys;*p!=NULL;p++) {
		Fl_ToggleNode *node = new Fl_ToggleNode(*p,1);
		if(needsub) {
			needsub = false;
			tree->add_sub(node);
		}
		else {
			tree->add_next(node);
		}		
		tree->close(node);
		bool empty = uiPopulateTreeWidgetRec(tree,reg.getSubFolder(*p),prefix + *p + ".");
		node->can_open(!empty);

	}
	
	if(!needsub)
		tree->traverse_up();

	return needsub;
}

void uiPopulateTreeWidget(Fl_ToggleTree *tree,Registry &reg) {
	
	tree->edit_on_reselect(false);
	if(tree->first()) {
		tree->traverse_start();
		tree->clear();
	}
	tree->traverse_start();

	// Get all the folder keys
	uiPopulateTreeWidgetRec(tree,reg,string(""));
}

void uiTreeSelection(Fl_ToggleTree *tree, void*) {
	Fl_ToggleNode* current = tree->selected();
	if(current) {
		string *fullKey = (string *)current->data();
		if(fullKey) {
			inTreeStringValue->value(regOptions.getStringValue(fullKey->c_str(),NULL));
			bTreeUpdate->activate();
			return;
		}
	}
	inTreeStringValue->value("");
	bTreeUpdate->deactivate();
}

void uiTreeUpdate(Fl_Button *button,void *) {
	Fl_ToggleNode* current = treeSettings->selected();
	if(current) {
		string *fullKey = (string *)current->data();
		if(fullKey) {
			regOptions.setStringValue(fullKey->c_str(),inTreeStringValue->value());
			return;
		}
	}
}

void uiTreeRefresh(Fl_Button *button,void *) {
	uiPopulateTreeWidget(treeSettings,regOptions);
}

/******************************************************
* Main
c:\dev\images\pelvis.im c:\dev\images\pelvis_nice.model
******************************************************/
char *lookatme = (char *)0x029ac5f0;

void main(int argc,char *argv[]) {      
	// initialize the image
	// baseImage = cimage_init();
	
	// Initialize all the FL junk
	make_window();

	// Show main window
	uiScript.runCommand("show winMain 10 32");
	
	// Position all other sig. windows
	winPrimitive->position(655,32);
	winScaleSpace->position(655,530);
	winCoreWizard->position(985,32);
	
	// Read image if there are args
	makeNewModel(1,8,false,0.25);

	// Initialize the spline
/*
	spline = new SplineObject();
	SplineCurve *c1 = spline->addCurve(9,0.25);
	c1->updateControl(4,Vector2D(0.5,0.45),0.1);
	spline->addBranch(c1,4,5);
	splineSample = new SplineSample(spline,20);
*/

	spline = new SplineObject();
	splineSample = new ArcSplineSample(spline,40);
/*
	SplineCurve *c1 = spline->addCurve(14,0.25);
	// c1->updateControl(4,Vector2D(0.5,0.45),0.1);
	spline->addBranch(c1,3,4);
	spline->addBranch(spline->curves[1],7,4);
	// spline->addBranch(spline->curves[2],5,4);
	splineSample = new SplineSample(spline,40);
*/
	// Read default options
	try {
		regOptions.readFromFile("dsltool.ini");
	} catch(...) {
	}

	regOptions.setFlagAddIfNotFound(true);
	uiPopulateTreeWidget(treeSettings,regOptions);
	
	// Set the widgets in the default options file
	list<string> keys;
	regOptions.collectAllKeys(keys);
	
	uiScript.startRecording();
	for(list<string>::iterator i=keys.begin();i!=keys.end();i++) {
		char *key = (char *) i->c_str();
		char *val = (char *) regOptions.getStringValue(key,NULL);
		cout << key << "=" << val << endl;
		uiScript.uiCommand("set %s \"%s\"",key,val);
	}
	uiScript.stopRecording();
	uiScript.run();

	// Check if there has been a script submitted
	for(int r=1;r<argc;r++) {
		if(strstr(argv[r],".model"))
			readModel(argv[r],NULL);
		else if(strstr(argv[r],".im"))
			readImage(argv[r]);
		else if(strstr(argv[r],".pgm"))
			readImage(argv[r]);	
		else if(strstr(argv[r],".script")) {
			UserInterfaceScript uiTempScript(uiScript);
			uiTempScript.load(argv[r]);
			uiTempScript.run(argc-r-1,argv+r+1);
		}
	}
	
	Fl::visual(FL_RGB);
	Fl::run();
	
	// Write the current options
	regOptions.writeToFile("dsltool.ini");   
}



