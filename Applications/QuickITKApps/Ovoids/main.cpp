#include "blobmodl.h"

#include "ui.h"

#include <FL/glut.H>
#include <FL/Fl_Color_Chooser.H>
#include <FL/Fl_File_Chooser.H>
#include "fracviewer.h"
#include "glcode.h"
#include <Registry/registry.h>

/**********************************************************************
 * GLOBAL VARIABLES                         									 *
 **********************************************************************/
Model *model;
ParabolicSurface *psurf;
CausticSurface *csurf;
RidgeSurface *rsurf;

ISurface *currSurface;


/**********************************************************************
 * LOCAL FUNCTION DEFINITIONS               									 *
 **********************************************************************/
void fillChoiceBox();

/**********************************************************************
 * MAIN                                   									 *
 **********************************************************************/
int main (int argc, char *argv[])
{
  // Must do this first, otherwize no GL will work
  initGLDisplay();

  // Create a couple ovoids andplace them in a model
  model = new Model();

  Ovoid *o1 = new Ovoid(2,2);
  //o1->setAffineTransform(2,3,4,-1.9,0,0,30,20,-20);
  o1->setAffineTransform(0.5,0.5,0.5,0,-0.8,0,0,0,0);

  Ovoid *o2 = new Ovoid(2,2);
  //o2->setAffineTransform(1,1,1,1.9,0,0,34,12,12);
  o2->setAffineTransform(0.5,0.5,0.5,0,0.8,0,0,0,0);

  model->ovoids->append(o1);
  model->ovoids->append(o2);

  iSurfList->append(o1);
  iSurfList->append(o2);
  iSurfList->append(model);

  model->recompute();

  psurf = new ParabolicSurface(model);
  iSurfList->append(psurf);

  csurf = new CausticSurface(model);
  iSurfList->append(csurf);

  rsurf = new RidgeSurface(model);
  iSurfList->append(rsurf);

  currSurface = model;

  // Initialize GUI
  createUI();
  fillChoiceBox();

  // Update the radio buttons
  uiCSChange(NULL,NULL);

  winMain->show();

  // Initialize GLUT and container window
  glutMainLoop();

  return 0;
}

/*****************************************************************
  All kinds of UI functions
 *****************************************************************/
void fillInputFields(Ovoid *ovoid) {
  inputXR->value(ovoid->Rx); 
  inputYR->value(ovoid->Ry); 
  inputZR->value(ovoid->Rz); 
  inputXT->value(ovoid->Tx); 
  inputYT->value(ovoid->Ty); 
  inputZT->value(ovoid->Tz); 
  inputXS->value(log(ovoid->Sx)); 
  inputYS->value(log(ovoid->Sy)); 
  inputZS->value(log(ovoid->Sz)); 

  sliderAlpha->value(log(ovoid->getAlpha()) / log(2.0) ); 
  sliderBeta->value( log(ovoid->getBeta() ) / log(2.0) ); 
}


void blobAffineChange() {
  // Compute a matrix
  int o = choiceOvoid->value();
  Ovoid *ovoid = (Ovoid *) model->ovoids->get(o);

  ovoid->setAffineTransform(exp(inputXS->value()),
    exp(inputYS->value()),
    exp(inputZS->value()),
    inputXT->value(),
    inputYT->value(),
    inputZT->value(),
    inputXR->value(),
    inputYR->value(),
    inputZR->value());

  glutPostRedisplay();
}

int currentOvoid = 0;

void ovoidChanged() {
  // All surfaces become invalid
  model->invalidate();
  psurf->invalidate();
  csurf->invalidate();
  rsurf->invalidate();

  // Update the radio buttons
  uiCSChange(NULL,NULL);

  // redisplay
  glutPostRedisplay();
}


void affineSliderChange(Fl_Value_Slider*a, void*b) {
  blobAffineChange();
  ovoidChanged();
}

void fillChoiceBox() {
  Ovoid *ovoid;

  choiceOvoid->clear();
  for(int i=0;i<model->ovoids->getSize();i++)  {
    ovoid = (Ovoid *)model->ovoids->get(i);
    char text[256];
    sprintf(text,"Blob %d, (%lg,%lg)",i,ovoid->getAlpha(),ovoid->getBeta());
    choiceOvoid->add(text);
  }

  choiceOvoid->value(currentOvoid);
}

void parameterSliderChange(Fl_Value_Slider *s, void*junk) {
  int o = choiceOvoid->value();
  Ovoid *ovoid = (Ovoid *) model->ovoids->get(o);

  double alpha = pow(2.0,sliderAlpha->value());
  double beta = pow(2.0,sliderBeta->value());
  ovoid->setParms(alpha,beta);

  blobAffineChange();
  ovoidChanged();
}

void ovoidSelected(Fl_Choice *choice,void *junk) {
  int o = choiceOvoid->value();
  Ovoid *ovoid = (Ovoid *) model->ovoids->get(o);
  Ovoid *lastOvoid = (Ovoid *) model->ovoids->get(currentOvoid);

  if(!ovoid)
    return;

  fillInputFields(ovoid);

  lastOvoid->setColor(0.6,0.6,0.6,1);
  ovoid->setColor(0,0,0.6,1);

  currentOvoid = o;
}

bool showOvoids = true;

void uiEditOvoids(Fl_Button*, void*) {
  fillChoiceBox();
  winCompEditor->show();
  winMain->hide();

  // All surfaces become invalid
  model->invalidate();
  psurf->invalidate();
  csurf->invalidate();
  rsurf->invalidate();

  // Show all ovoids
  for(int i=0;i<model->ovoids->getSize();i++) {
    Ovoid *ovoid = (Ovoid *) model->ovoids->get(i);
    ovoid->setDisplayed(true);
  }

  // Set all fields to current ovoid
  ovoidSelected(NULL,NULL);

  // redisplay
  glutPostRedisplay();
}

void uiShowHideOvoids(Fl_Menu_*, void*) {
  showOvoids = !showOvoids;

  for(int i=0;i<model->ovoids->getSize();i++) {
    Ovoid *ovoid = (Ovoid *) model->ovoids->get(i);
    ovoid->setDisplayed(showOvoids);
  }

  glutPostRedisplay();
}

void blobOKPressed(Fl_Return_Button *button,void *junk) {
  Ovoid *lastOvoid = (Ovoid *) model->ovoids->get(currentOvoid);
  lastOvoid->setColor(0.6,0.6,0.6,1);

  winCompEditor->hide();
  currentOvoid = 0;

  // Recompute the model
  model->setDisplayed(true);

  // Update the radio buttons
  uiCSChange(NULL,NULL);

  // Hide all ovoids if flag says so
  if(!showOvoids) {
    for(int i=0;i<model->ovoids->getSize();i++) {
      Ovoid *ovoid = (Ovoid *) model->ovoids->get(i);
      ovoid->setDisplayed(false);
    }
  }

  // redisplay
  glutPostRedisplay();

  // Show the main window
  winMain->show();
}

void uiDeleteOvoid(Fl_Button*, void*) {
  Ovoid *ovoid = (Ovoid *) model->ovoids->get(currentOvoid);

  // Remove it from model
  model->ovoids->set(currentOvoid,NULL);
  model->ovoids->compact();

  // Remove it from surface list
  iSurfList->set(iSurfList->getIndex(ovoid),NULL);
  iSurfList->compact();

  // Delete it
  delete ovoid;

  // Rebuild the choice box
  fillChoiceBox();
  choiceOvoid->value(0);
  currentOvoid = 0;
  ovoidSelected(NULL,NULL);

  // Ovoid configuratoin changed
  ovoidChanged();
}

void createOvoid(Fl_Button*b, void*j) {
  Ovoid *o = new Ovoid(2,2);
  model->ovoids->append(o);

  // Update the window
  fillChoiceBox();
  choiceOvoid->value(model->ovoids->getSize()-1);
  ovoidSelected(NULL,NULL);

  // Add ovoid to list
  iSurfList->append(o);

  ovoidChanged();
}

void uiRebuildModel(Fl_Button*, void*) {
  model->recompute();
  glutPostRedisplay();
}

void uiCSChange(Fl_Choice*, void*) {
  switch(choiceSurface->value()) {
  case 0 : currSurface = model;break;
  case 1 : currSurface = psurf;break;
  case 2 : currSurface = rsurf;break;
  case 3 : currSurface = csurf;break;

  }

  // Set radio buttons 
  rbCSNotComputed->value(0);
  rbCSNotDisplayed->value(0);
  rbCSWireframe->value(0);
  rbCSSolid->value(0);

  if(!currSurface->computed) 
    rbCSNotComputed->value(1);
  else if(!currSurface->displayed)
    rbCSNotDisplayed->value(1);
  else if(!currSurface->solid)
    rbCSWireframe->value(1);
  else
    rbCSSolid->value(1);

  inputPolygonSize->value(currSurface->getTetraSize());
}	

  void uiCSDisplayMode(bool display,bool solid) {
    if(solid)
      currSurface->drawSolid();
    else
      currSurface->drawWireframe();

    currSurface->setDisplayed(display);

    // Update the radio buttons
    uiCSChange(NULL,NULL);

    glutPostRedisplay();
  }

void uiCSColor(Fl_Button*, void*) {
  double r = currSurface->color[0],g = currSurface->color[1],b = currSurface->color[2];
  fl_color_chooser("Surface Color", r, g, b);
  currSurface->setColor(r,g,b,1);
  glutPostRedisplay();
}

void uiLoadModel(Fl_Menu_*, void*) {
  char *file = fl_file_chooser("Load Model","*.blobs",NULL);
  if(!file)
    return;

  Registry reg;
  reg.readFromFile(file);

  model->read(&reg.getSubFolder("model"));

  delete psurf;
  psurf = new ParabolicSurface(model);

  delete csurf;
  csurf = new CausticSurface(model);

  delete rsurf;
  rsurf = new RidgeSurface(model);

  iSurfList->clear();
  for(int o=0;o<model->ovoids->getSize();o++) {
    iSurfList->append(model->ovoids->get(o));
    iSurfList->append(model);
    iSurfList->append(psurf);
    iSurfList->append(csurf);
    iSurfList->append(rsurf);
  }

  // Update the radio buttons
  model->recompute();
  uiCSChange(NULL,NULL);
  fillChoiceBox();

  glutPostRedisplay();
}

void uiSaveModel(Fl_Menu_*, void*) {
  char *file = fl_file_chooser("Load Model","*.blobs",NULL);
  if(!file)
    return;

  Registry reg;
  model->write(&reg.getSubFolder("model"));
  reg.writeToFile(file);
}

void uiCSPolygonResize(Fl_Value_Input*, void*) {
  currSurface->setTetraSize(inputPolygonSize->value());
  currSurface->setDisplayed(true);
  glutPostRedisplay();
}


void fillCurrentPointInfo() {
  if(currentPoint) {
    outKappa1->value(currentPoint->kappa1);
    outKappa2->value(currentPoint->kappa2);
    outGradMag->value(currentPoint->gradMag);
  }
  else {
    outKappa1->value(0);
    outKappa2->value(0);
    outGradMag->value(0);
  }
}

void uiFlyBy(Fl_Menu_*, void*) {
  agvSwitchMoveMode(FLYING);
}

void uiPolarMovement(Fl_Menu_*, void*) {
  agvSwitchMoveMode(POLAR);
}
