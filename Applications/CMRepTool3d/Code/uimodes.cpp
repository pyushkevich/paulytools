#include "glui.h"

void SurfaceSelectionMH::drawText() {

  rndTextOut.clear();
  rndTextOut.addLine("Surface Editing Mode");

  if (su < 0)
    {
    if (showInstructions)
      {
      rndTextOut.addLine("This mode allows you to push and pull points on the surface");
      rndTextOut.addLine("Move the mouse over the surface");
      }
    else
      {
      rndTextOut.addLine("");
      }
    }
  else if (!clicked)
    {

    if (showInstructions)
      {
      rndTextOut.addLine("Under mouse you see local tangent plane and the normal vector");
      rndTextOut.addLine("Click and hold the RIGHT mouse button to begin deforming the surface");
      }
    else
      {
      // The current point's u,v coordinates
      double u = splineData->getGrid()->uv(0,su);
      double v = splineData->getGrid()->uv(1,sv);

      rndTextOut.addLine("U = %.03g,  V = %.03g",u,v);
      rndTextOut.addLine("X = [%.03g, %.03g, %.03g]",(double)mp.F[0],(double)mp.F[1],(double)mp.F[2]);
      rndTextOut.addLine("R = %.03g",(double)mp.F[3]);
      rndTextOut.addLine("K11 = %.03g,  K12 = %.03g",(double)mp.bp[0].kappa1,(double)mp.bp[0].kappa2);
      rndTextOut.addLine("K21 = %.03g,  K22 = %.03g",(double)mp.bp[1].kappa1,(double)mp.bp[1].kappa2);
      }

    }
  else
    {

    if (showInstructions)
      {
      rndTextOut.addLine("Move the mouse while keeping the right putton pressed."  );
      rndTextOut.addLine("The red line represents the direction in which the surface will be deformed");
      rndTextOut.addLine("Press SPACEBAR to toggle between moving along the normal and in the tangent plane.");
      rndTextOut.addLine("Release the RIGHT mouse button to deform the surface.");
      }
    else
      {
      // The current point's u,v coordinates
      double u = splineData->getGrid()->uv(0,su);
      double v = splineData->getGrid()->uv(1,sv);

      rndTextOut.addLine("U = %.03lg,  V = %.03lg",u,v);
      }
    }

  if (showInstructions)
    {
    rndTextOut.addLine("");
    rndTextOut.addLine("Other Commands:");
    rndTextOut.addLine("   0 : cancel current deformation");
    rndTextOut.addLine("   C : toggle control polygon");
    rndTextOut.addLine("   H : hide these instructions");
    rndTextOut.addLine("   F : flatten the surface");
    rndTextOut.addLine("  ESC : exit the program");
    }
  else
    {
    rndTextOut.addLine("Press 'H' for instructions.");
    }

}

void SurfaceSelectionMH::onDraw() {
  if (GLDisplayDriver::displayMode == GLDisplayDriver::EYE)
    {
    AbstractImage3D *image = imaging.getImage();
    if (image != NULL)
      {
      int prfInt[1024],prfSize;
      float prfPos[1024];

      SMLVec3f XI,UI[2];
      XI = TransformPoint(image->SI,mp.X());
      UI[0] = TransformVector(image->SI,mp.bp[0].N * mp.R());
      UI[1] = TransformVector(image->SI,mp.bp[1].N * mp.R());

      glPushAttrib(GL_LIGHTING_BIT);
      glDisable(GL_LIGHTING);

      // Place somewhere visible
      glPushMatrix();
      glTranslated(600,0,0);
      glScaled(150,50,1);

      // Draw the box around this thing           
      glBegin(GL_LINES);
      glColor3d(0.2,0.7,0.9);
      glVertex2d(-1.6,0.2);
      glVertex2d(1.6,0.2);
      glVertex2d(-1.0,0.1);
      glVertex2d(-1.0,0.9);
      glVertex2d(1.0,0.1);
      glVertex2d(1.0,0.9);
      glVertex2d(0.0,0.0);
      glVertex2d(0.0,1.0);
      glEnd();

      for (int d=0;d<2;d++)
        {
        parseLine3D(image,XI.data_block(),UI[d].data_block(),1.5f,prfInt,prfPos,prfSize);
        float sign = -1.0 + d*2.0;

        glBegin(GL_LINE_STRIP);
        glColor3d(0.8,0.9,0.3);
        //glVertex2d(sign*prfPos[0],0.2);
        for (int i=0;i<prfSize-1;i++)
          {
          glVertex2d(sign*prfPos[i],0.6 * (prfInt[i] / 256.0) + 0.2);
          glVertex2d(sign*prfPos[i+1],0.6 * (prfInt[i] / 256.0) + 0.2);
          }
        //glVertex2d(sign*prfPos[i],0.2);
        glEnd();
        }

      glPopMatrix();
      glPopAttrib();
      }
    }
  else
    {
    PointMovementMH::onDraw();
    }

}

bool SurfaceSelectionMH::handleKeys(unsigned char key,int x,int y) {
  if (PointMovementMH::handleKeys(key,x,y))
    {
    return true;
    }

  float u = splineData->getGrid()->uv(0,su);
  float v = splineData->getGrid()->uv(1,sv);

  if (key=='x' || key=='X')
    {
    if (geodesicRenderer.needsStart())
      {
      // Mark the starting point for geodesic
      geodesicRenderer.setStart(u,v);
      }
    else
      {
      // Mark the starting point for geodesic
      geodesicRenderer.setEnd(u,v);
      }

    return true;
    }
  return false;
}


PointMovementMH::PointMovementMH() : rndTextOut(fntCourier18,GLColor(0.2,0.6,0.7,0.8)) {

}

void PointMovementMH::drawSelectedPoint() {
  glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT | GL_COLOR_BUFFER_BIT);

  glDisable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);

  // Draw a planar patch around the normal
  glPushMatrix();

  // Create a rotation/translation matrix
  SMLMatrix4f M;
  M.set_identity();
  for (int i=0;i<3;i++)
    {
    M.put(2,i,mp.N()[i]);
    M.put(1,i,XV[i]);
    M.put(0,i,XU[i]);
    M.put(3,i,mp.X()[i]);
    }

  glTranslated(-0.5,-0.5,-0.5);
  glMultMatrixf(M.data_block());
  glScaled(0.1,0.1,0.1);

  // Save projection, model and viewport matrices
  glGetDoublev(GL_MODELVIEW_MATRIX,MM);
  glGetDoublev(GL_PROJECTION_MATRIX,MP);
  glGetIntegerv(GL_VIEWPORT,MV);

  // Draw some lines
  static int dl = -1;
  if (dl < 0)
    {
    dl = glGenLists(1);
    glNewList(dl,GL_COMPILE);

    glBegin(GL_LINES);
    for (int u=-5;u<=5;u++)
      {
      if (u==0)
        glColor3d(1,1,0.5);
      else
        glColor3d(0.5,0.5,0.25);

      double U = u * 0.2;

      glVertex3d(U,-1.0,0);
      glVertex3d(U,1.0,0);
      glVertex3d(-1.0,U,0);
      glVertex3d(1.0,U,0);

      }

    glColor3d(1,1,0.5);
    glVertex3d(0,0,0);
    glVertex3d(0,0,1);

    glEnd();

    glEndList();
    }

  if (mode==NORMAL || mode==TANGENT)
    {
    glCallList(dl);

    // Fat lines
    glLineWidth(2);

    // Draw the line to T
    glBegin(GL_LINES);
    glColor3d(1,0,0);
    glVertex3d(0,0,0);
    glVertex(T);
    glEnd();

    // Draw dotted lines
    glLineWidth(1);
    glEnable(GL_LINE_STIPPLE);
    glLineStipple(3,0xAAAA);

    glBegin(GL_LINE_STRIP);
    glVertex3d(0,0,0);
    glVertex3d(T[0],T[1],0);
    glVertex(T);
    glEnd();
    }

  // Draw the fitted sphere
  if (mode==RADIUS && mp.sinTheta2 >= 0)
    {
    // SMLVec3f X1 = (XB[0] - X) * 10 ;
    // SMLVec3f X2 = (XB[1] - X) * 10 ;

    SMLVec3f X1 = (mp.bp[0].X - mp.X()) * 10.0f;
    SMLVec3f X2 = (mp.bp[1].X - mp.X()) * 10.0f;

    // glColor3d(1,0,1);
    glColor3d(0.5,0.5,0.25);
    glutWireSphere(mp.R()*10,16,16);

    if (Rnew != mp.R())
      {
      glutWireSphere(fabs(Rnew)*10,16,16);
      }
    }

  glPopMatrix();  

  glPopAttrib();
}

void PointMovementMH::onDraw() {
  // Handle mouse motion in 'predisplay' mode
  if (GLDisplayDriver::displayMode == GLDisplayDriver::PRE)
    {
    // Handle mouse motion (if any)
    if (doPassiveOnDraw)
      {
      testMousePosition(mptX,mptY);
      doPassiveOnDraw = false;

      // Set the transformed point to zero
      T[0] = T[1] = T[2] = 0.0;
      Rnew = mp.R();

      if (havePoint())
        {
        computeFrame();
        }
      }

    else if (doMotionOnDraw)
      {
      doHandleMotion(mptX,mptY);
      doMotionOnDraw = false;
      }
    }
  else if (GLDisplayDriver::displayMode == GLDisplayDriver::NORM)
    {
    // Draw the text
    glColor3d(0.2,0.5,0.7);
    drawText();

    if (havePoint())
      {
      drawSelectedPoint();
      }
    }
}

void PointMovementMH::scheduleMousePositionTest(int x,int y) {
  mptX = x;
  mptY = y;
  doPassiveOnDraw = true;
  glutPostRedisplay();
}


void PointMovementMH::start() {
  GLDisplayDriver::addRenderer(&rndTextOut,GLDisplayDriver::EYE);
  GLDisplayDriver::addRenderer(this,GLDisplayDriver::NORM | GLDisplayDriver::PRE | GLDisplayDriver::EYE);
  GLDisplayDriver::addListener(this,GLDisplayDriver::PASSIVE);

  GLDisplayDriver::addRenderer(&geodesicRenderer,GLDisplayDriver::NORM);

  // rndSpline->setSurfaceMode(BSplineRenderer::SEETHRU | BSplineRenderer::CONTROL);

  // Motion positions are reset
  mptX = -1,mptY = -1;

  clicked = false;
  mode = NORMAL;
  showInstructions = false;
}

void PointMovementMH::stop() {
  GLDisplayDriver::removeRenderer(&rndTextOut,GLDisplayDriver::EYE);
  GLDisplayDriver::removeRenderer(this,GLDisplayDriver::PRE | GLDisplayDriver::EYE | GLDisplayDriver::NORM);
  GLDisplayDriver::removeListener(this,GLDisplayDriver::PASSIVE);

  GLDisplayDriver::removeRenderer(&geodesicRenderer,GLDisplayDriver::NORM);

  // rndSpline->setSurfaceMode(BSplineRenderer::SURFACE);
}

bool PointMovementMH::handlePassiveMotion(int x,int y) {

  // If the user is just moving the mouse around, track its position on the surface           
  // testMousePosition(x,GLDisplayDriver::height - y);
  scheduleMousePositionTest(x,GLDisplayDriver::height - y);
  return false;
}

void SurfaceSelectionMH::computeFrame() {
  spline->interpolateMedialPoint02(*splineData->getGrid(),su,sv,mp);
  spline->interpolateBoundary(mp);
  // spline->computeCurvatures(mp);
  XU = mp.Xu();
  XV = mp.Xv();
  XU.normalize();
  XV.normalize();
}

void SurfaceSelectionMH::testMousePosition(int x,int y) {
  rndSpline->findSampleUV(x,y,su,sv);
}

void SurfaceSelectionMH::applyMovement() {
  if (mode==RADIUS)
    {
    float Rold;

    spline->interpolateGridPoint(*splineData->getGrid(),su,sv,0,0,3,3,&Rold);

    float dr = Rnew-Rold;
    spline->push(*splineData->getGrid(),su,sv,3,3,&dr);
    }
  else
    {
    // SMLVec3f V = Xu * T.x * 0.1 + Xv * T.y * 0.1 + N * T.z * 0.1;
    SMLVec3f V = (XU * T[0] +XV * T[1] + mp.N() * T[2]) * 0.1f;

    SMLVec3f Xn;

    spline->interpolateGridPoint(*splineData->getGrid(),su,sv,0,0,0,2,Xn.data_block());

    V = V - (Xn - mp.X());

    spline->push(*splineData->getGrid(),su,sv,0,2,V.data_block());
    }
}

bool PointMovementMH::handleKeys(unsigned char key,int x,int y) {
  // cout << "Key: " << key << (glutGetModifiers() & GLUT_ACTIVE_SHIFT != 0) <<  "\n";
  if (key=='t' || key=='T')
    {
    mode = TANGENT;
    motionStart[0] = x;
    motionStart[1] = GLDisplayDriver::height - y;
    T0 = T;
    glutPostRedisplay();
    return true;
    }
  else if (key=='r' || key=='R')
    {
    mode = RADIUS;
    motionStart[0] = x;
    motionStart[1] = GLDisplayDriver::height - y;
    Rnew = mp.R();
    glutPostRedisplay();
    return true;
    }
  else if (key=='n' || key=='N')
    {
    mode = NORMAL;
    motionStart[0] = x;
    motionStart[1] = GLDisplayDriver::height - y;
    T0 = T;
    glutPostRedisplay();
    return true;
    }
/*
    else if(key=='C' || key=='c') {
        rndSpline->setSurfaceMode(rndSpline->getSurfaceMode() ^ BSplineRenderer::CONTROL);
        glutPostRedisplay();
    }
    else if(key=='1') {
        rndSpline->setSurfaceMode(rndSpline->getSurfaceMode() ^ BSplineRenderer::IBSURFACE_LEFT);
        glutPostRedisplay();
    }
    else if(key=='2') {
        rndSpline->setSurfaceMode(rndSpline->getSurfaceMode() ^ BSplineRenderer::IBSURFACE_RIGHT);
        glutPostRedisplay();
    }*/
  else if (key=='H' || key=='h')
    {
    showInstructions = !showInstructions;
    glutPostRedisplay();
    return true;
    }

  else if (key=='F' || key=='f')
    {
    // Flatten the surface
    for (int i=0;i<=spline->dim(0);i++)
      {
      for (int j=0;j<=spline->dim(1);j++)
        {
        spline->setControl(i,j,2,0.5);
        }
      }

    glutPostRedisplay();
    }

  else if (key=='0' || key==')')
    {

    if (clicked)
      {
      T[0] = T[1] = T[2] = 0;
      Rnew = mp.R();
      motionStart[0] = x;
      motionStart[1] = GLDisplayDriver::height - y;
      T0 = T;
      }

    glutPostRedisplay();
    return true;
    }

  else if (key == 'C' || key == 'c')
    {
    GLDisplayDriver::center = mp.X() - SMLVec3f(0.5,0.5,0.5);
    glutPostRedisplay();
    return true;
    }


  return false;
}

void PointMovementMH::scheduleMouseMotion(int x,int y) {
  mptX = x;
  mptY = y;
  doMotionOnDraw = true;
}

void PointMovementMH::doHandleMotion(int x,int y) {
  if (clicked && havePoint())
    {
    if (mode == TANGENT)
      {

      double W0[3],WX1[3],WY1[3];

      // Project the normal vector on viewport plane
      gluProject(0,0,0,MM,MP,MV,W0,W0+1,W0+2);
      gluProject(1,0,0,MM,MP,MV,WX1,WX1+1,WX1+2);
      gluProject(0,1,0,MM,MP,MV,WY1,WY1+1,WY1+2);

      SMLVec2f WX(WX1[0]-W0[0],WX1[1]-W0[1]);
      SMLVec2f WY(WY1[0]-W0[0],WY1[1]-W0[1]);

      SMLVec2f DM = SMLVec2f(x,y) - motionStart;

      float lx = dot_product(WX,DM) / WX.squared_magnitude();
      float ly = dot_product(WY,DM) / WY.squared_magnitude();

      // Append the z value by len
      T[0] = T0[0] + lx;
      T[1] = T0[1] + ly;

      // Warp the surface, what the hell
      applyMovement();
      // rndSpline->updateSurface();

      // glutPostRedisplay();
      // return true;     
      }
    else if (mode==NORMAL)
      {
      double W0[3],W1[3];

      // Project the normal vector on viewport plane
      gluProject(0,0,0,MM,MP,MV,W0,W0+1,W0+2);
      gluProject(0,0,1,MM,MP,MV,W1,W1+1,W1+2);

      SMLVec2f W(W1[0]-W0[0],W1[1]-W0[1]);
      SMLVec2f DM = SMLVec2f(x,y) - motionStart;

      // W.Normalize();
      // DM.Normalize();

      float len = dot_product(W,DM) / W.squared_magnitude();

      // Append the z value by len
      T[2] = T0[2] + len;

      // Warp the surface, what the hell
      applyMovement();
      // rndSpline->updateSurface();

      // glutPostRedisplay();

      // return true;
      }
    else if (mode == RADIUS)
      {
      float dy = 10.0f * (y - motionStart[1]) / GLDisplayDriver::height;
      float dx = 10.0f * (x - motionStart[0]) / GLDisplayDriver::width;
      //Rnew = r * pow(1.5,dy);
      Rnew = mp.R() + dy / 30.0 + dx / 3;
      applyMovement();
      }
    }
}

bool PointMovementMH::handleMotion(int x,int y) {
  y = GLDisplayDriver::height - y;
  if (clicked && havePoint())
    {
    scheduleMouseMotion(x,y);
    glutPostRedisplay();
    return true;
    }

  return false;
}

bool PointMovementMH::handleButton(int button,int state,int x,int y) {

  if (button == GLUT_RIGHT_BUTTON && havePoint())
    {
    if (state == GLUT_DOWN)
      {
      clicked = true;

      motionStart[0] = x;
      motionStart[1] = GLDisplayDriver::height- y;
      T0 = T;

      glutPostRedisplay();

      return true;
      }
    else if (state == GLUT_UP)
      {
      // Apply it!
      applyMovement();

      // Surface changed
      // rndSpline->updateSurface();

      glutPostRedisplay();

      clicked = false;
      return true;
      }
    }

  return false;
}


void ControlPointEditMH::drawText() {

  const char* frameName[] = {"control polygon's","surface local","XYZ"};

  rndTextOut.clear();
  rndTextOut.addLine("Control Point Mode");

  if (!havePoint())
    {
    if (showInstructions)
      {
      rndTextOut.addLine("This mode allows you to move the B-spline control points");
      rndTextOut.addLine("Move the mouse near a control point to begin moving it");
      }
    else
      {
      rndTextOut.addLine("");
      }
    }
  else if (!clicked)
    {

    if (showInstructions)
      {
      rndTextOut.addLine("RIGHT - click the mouse to begin control point movement");
      }
    else
      {
      rndTextOut.addLine("Selected point (%d, %d)",ci,cj);
      rndTextOut.addLine("Using %s coordinate frame",frameName[frameIndex]);
      }

    }
  else
    {

    if (showInstructions)
      {
      rndTextOut.addLine("Move the mouse while keeping the right putton pressed."  );
      rndTextOut.addLine("The red line represents the direction in which the control point will be moved");
      rndTextOut.addLine("Press SPACEBAR to toggle between moving along the normal and in the tangent plane.");
      rndTextOut.addLine("Release the RIGHT mouse button to deform the surface.");
      }
    else
      {
      rndTextOut.addLine("Selected point (%d, %d)",ci,cj);
      rndTextOut.addLine("Using %s coordinate frame",frameName[frameIndex]);
      }
    }

  if (showInstructions)
    {
    rndTextOut.addLine("");
    rndTextOut.addLine("Other Commands:");
    rndTextOut.addLine("");
    rndTextOut.addLine("Other Commands:");
    rndTextOut.addLine("   0 : cancel current deformation");
    rndTextOut.addLine("   C : toggle control polygon");
    rndTextOut.addLine("   H : hide these instructions");
    rndTextOut.addLine("   F : flatten the surface");
    rndTextOut.addLine("   X : toggle coordinate frame [control point/local/xyz]");
    rndTextOut.addLine("  ESC : exit the program");
    }
  else
    {
    rndTextOut.addLine("Press 'H' for instructions.");
    }

}

void ControlPointEditMH::testMousePosition(int x,int y) {
  rndSpline->findControlPoint(x,y,ci,cj);
}

void ControlPointEditMH::applyMovement() {
  if (mode==RADIUS)
    {
    spline->setControl(ci,cj,3,Rnew);   
    }
  else
    {
    SMLVec3f V = mp.X() + (XU * T[0] + XV * T[1] + mp.N() * T[2]) * 0.1f;
    spline->setControl(ci,cj,0,V[0]);
    spline->setControl(ci,cj,1,V[1]);
    spline->setControl(ci,cj,2,V[2]);
    }
}

void ControlPointEditMH::computeFrame() {
  // The position X
  mp.F = spline->getControl4f(ci,cj);
  mp.sinTheta2 = 1;

  if (frameIndex == POLYGON)
    {
    // The initial Xu and Xv vectors
    mp.Fu = ((ci > 0) 
             ? (spline->getControl4f(ci,cj) - spline->getControl4f(ci-1,cj)) * 0.5f 
             : (spline->getControl4f(ci+1,cj) - spline->getControl4f(ci,cj)) * 0.5f) 
            + ((ci < spline->dim(0)) 
               ? (spline->getControl4f(ci+1,cj) - spline->getControl4f(ci,cj)) * 0.5f 
               : (spline->getControl4f(ci,cj) - spline->getControl4f(ci-1,cj)) * 0.5f);

    mp.Fv = ((cj > 0) 
             ? (spline->getControl4f(ci,cj) - spline->getControl4f(ci,cj-1)) * 0.5f 
             : (spline->getControl4f(ci,cj+1) - spline->getControl4f(ci,cj)) * 0.5f) 
            + ((cj < spline->dim(1)) 
               ? (spline->getControl4f(ci,cj+1) - spline->getControl4f(ci,cj)) * 0.5f 
               : (spline->getControl4f(ci,cj) - spline->getControl4f(ci,cj-1)) * 0.5f);

    // The normal vector
    mp.N() = vnl_cross_3d(mp.Xu(),mp.Xv());
    mp.N().normalize();

    // Make an orthogonal frame     
    mp.Xu().normalize();
    mp.Xv() = vnl_cross_3d(mp.N(),mp.Xu());
    XU = mp.Xu();XV = mp.Xv();
    }
  else if (frameIndex == LOCAL)
    {
    // For border points use same frame as for regular points
    int pi = ci > 0 ? (ci < spline->dim(0) ? ci : ci-1) : ci+1;
    int pj = cj > 0 ? (cj < spline->dim(1) ? cj : cj-1) : cj+1;

    // Find the index that corresponds to the knot
    int u = splineData->getGrid()->patchStart(0,pi-1);
    int v = splineData->getGrid()->patchStart(1,pj-1);

    // Compute the frame for this point
    SMLVec3f Xdummy,Ndummy;
    //spline->interpolateMedialSurfacePoint(*splineData->getGrid(),u,v,Xdummy,mp.Xu(),mp.Xv(),Ndummy,mp.N());
    spline->interpolateMedialPoint02(*splineData->getGrid(),u,v,mp);
    mp.Xu().normalize();
    mp.Xv() = vnl_cross_3d(mp.Xu(),mp.N());
    XU = mp.Xu();XV = mp.Xv();
    }
  else
    {
    mp.Xu() = SMLVec3f(1,0,0);
    mp.Xv() = SMLVec3f(0,1,0);
    mp.N() = SMLVec3f(0,0,1);
    XU = mp.Xu();XV = mp.Xv();
    }
}

bool ControlPointEditMH::handleKeys(unsigned char key,int x,int y) {
  if (PointMovementMH::handleKeys(key,x,y))
    {
    return true;
    }
  if (key=='x' || key=='X')
    {
    frameIndex = (Frame)((frameIndex+1) % 3);

    if (havePoint())
      {
      T = T0 = SMLVec3f(0,0,0);
      motionStart = SMLVec2f(0,0);

      computeFrame();
      glutPostRedisplay();
      }

    return true;
    }
  return false;
}

