#include "TracerMainWindow.h"
#include <FL/gl.h>
#include <GL/glu.h>
#include <Fl/Fl.H>
#include <Fl/fl_draw.H>

#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>

#include "Trackball.h"

void 
TracerMainWindow
::ComputeDisplayList()
{
  
  // Generate a display list number if needed
  if(m_DisplayList < 0)
    m_DisplayList = glGenLists(1);

  // Get the triangle strip information.
  vtkCellArray *triStrips = m_Data->GetSourceMesh()->GetStrips();
    
  // Get the vertex information.
  vtkPoints *verts = m_Data->GetSourceMesh()->GetPoints();

  // Get the normal information.
  vtkDataArray *norms = m_Data->GetSourceMesh()->GetPointData()->GetNormals();

  // Build display list
  glNewList(m_DisplayList,GL_COMPILE_AND_EXECUTE);

  // Data for getting consecutive points
  vtkIdType ntris = 0;
  vtkIdType npts;
  vtkIdType *pts;
  for ( triStrips->InitTraversal(); triStrips->GetNextCell(npts,pts); ) 
    {
    ntris += npts-2;
    glBegin( GL_TRIANGLE_STRIP );
    for (vtkIdType j = 0; j < npts; j++) 
      {
      // Some ugly code to ensure VTK version compatibility
      double vx = verts->GetPoint(pts[j])[0];
      double vy = verts->GetPoint(pts[j])[1];
      double vz = verts->GetPoint(pts[j])[2];
      double nx = norms->GetTuple(pts[j])[0];
      double ny = norms->GetTuple(pts[j])[1];
      double nz = norms->GetTuple(pts[j])[2];
      double l = 1.0 / sqrt(nx*nx+ny*ny+nz*nz);

      // Specify normal.
      glNormal3d(l*nx, l*ny, l*nz);

      // Specify vertex.
      glVertex3d(vx, vy, vz);
      }
    glEnd();
    }
  glEndList();

  // Clear the dirty flag
  m_DisplayListDirty = false;
}

void 
TracerMainWindow
::InitializeGL()
{
  // Set up the clear color
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glEnable( GL_DEPTH_TEST );

  // Set up the materials
  GLfloat matAmb[4] = { 0.01, 0.01, 0.01, 1.00};
  GLfloat matDiff[4] = { 0.65, 0.65, 0.65, 1.00};
  GLfloat matSpec[4] = { 0.20, 0.20, 0.20, 1.00};
  GLfloat matShine = 10.0;
  
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, matAmb);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, matDiff);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, matSpec);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, matShine);
  glEnable(GL_COLOR_MATERIAL);

  // Setup Lighting
  GLfloat light0Pos[4] = { 0.0f, 0.0f, 1.0f, 0.0};
  GLfloat light0Amb[4] = { 0.1f, 0.1f, 0.1f, 1.0f};
  GLfloat light0Dif[4] = { 0.9f, 0.9f, 0.9f, 1.0f};

  glLightfv(GL_LIGHT0, GL_POSITION, light0Pos);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,  light0Dif);
  glLightfv(GL_LIGHT0, GL_AMBIENT,  light0Amb);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  glEnable(GL_NORMALIZE);

  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);

  // glShadeModel(GL_FLAT);
  // glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
}

void 
TracerMainWindow
::draw()
{
  if( !valid() )
    {
    // Set the view port
    glViewport(0, 0, w(), h());

    // Set up the trackball
    m_Trackball.Reset();

    // Rotate a little bit, so we see the three axes
    m_Trackball.StartRot(12,10,20,20);
    m_Trackball.TrackRot(10,8,20,20);
    m_Trackball.StopRot();

    // Set up the projection matrix
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    // float zoom = m_Trackball.GetZoom();
    // float x = 1.0 / zoom;
    glOrtho( -1 , 1, -1, 1, -100, 100 ); 

    // Set up the model view matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    }

  // Initialize GL if needed
  if(m_GLStateDirty)
    {
    InitializeGL();
    m_GLStateDirty = false;
    }
  
  // Clear the screen 
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glFlush();

  // If there is no mesh, don't do any of the following
  if(m_Data->GetSourceMesh() != NULL)
    {
    // Set up the model view matrix
    glPushMatrix();
    SetUpModelMatrix();

    // Draw the display list
    glColor3d(0.4,0.5,0.3);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    
    // Create the display list if needed
    if(m_DisplayListDirty) ComputeDisplayList();
    else glCallList(m_DisplayList);

    // Set the attributes for line drawing
    glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT | GL_COLOR_BUFFER_BIT);
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(2.5);

    // Display the existing curves
    for(unsigned int iCurve = 0; iCurve<m_Data->GetNumberOfCurves();iCurve++)
      {
      // Choose a color for the curve
      if(iCurve == m_Data->GetCurrentCurve()) glColor3d(1,0.5,0.5);
      else glColor3d(0.5,0,0);

      // Draw the curve
      glBegin(GL_LINE_STRIP);
      for(unsigned int iPoint = 0;
        iPoint < m_Data->GetNumberOfCurvePoints(iCurve);iPoint++)
        {
        TracerData::Vec x = m_Data->GetCurvePoint(iCurve,iPoint);
        glVertex3dv(x.data_block());
        }
      glEnd();
      }

    // Display the path to the current point
    if(m_CurrentPoint != -1 && m_EditMode == TRACER)
      {
      glColor3d(1,1,0);
      glPointSize(4.0);

      list<vtkIdType> lPoints = m_Data->GetPathToPoint(m_CurrentPoint);
      list<vtkIdType>::const_iterator it = lPoints.begin();

      if(lPoints.size() > 1) 
        {
        glBegin(GL_LINE_STRIP);
        TracerData::Vec x = m_Data->GetCurvePoint(
          m_Data->GetCurrentCurve(),
          m_Data->GetNumberOfCurvePoints(m_Data->GetCurrentCurve())-1);
        glVertex3d(x[0],x[1],x[2]);
        }
      else 
        {
        glBegin(GL_POINTS);
        }
      
      while(it != lPoints.end())
        {
        TracerData::Vec x = m_Data->GetPointCoordinate(*it);
        glVertex3d(x[0],x[1],x[2]);
        ++it;
        }
      
      glEnd();
      }

    // Restore the attributes
    glPopAttrib();

    // Pop the matrix
    glPopMatrix();
    }

  // Flush the pipeline
  glFlush();
}

void 
TracerMainWindow
::SetUpModelMatrix()
{
  glLoadIdentity();

  // Update the screen geometry
  glTranslatef( m_Trackball.GetPanX(), m_Trackball.GetPanY(), 0.0 );
  glMultMatrixf( m_Trackball.GetRot() );

  // Compute the scaling factor based on the zoom level
  float xScale = 
    2.0 / (m_Trackball.GetZoom() * m_Data->GetSourceMesh()->GetLength());
  glScaled(xScale,xScale,xScale);

  // Now, translate for the center of the mesh
  double *xCenter = m_Data->GetSourceMesh()->GetCenter();
  glTranslated(-xCenter[0],-xCenter[1],-xCenter[2]);
}

void 
TracerMainWindow
::FindPointUnderCursor()
{
  // Find the point currently under the cursoro
  make_current(); // update GL state

  // Set up the modelview matrix
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();  
  SetUpModelMatrix();

  // OpenGL matrices to query
  double mvmatrix[16];
  double projmatrix[16];
  GLint viewport[4];
 
  // Query the matrices
  glGetIntegerv( GL_VIEWPORT, viewport );
  glGetDoublev( GL_MODELVIEW_MATRIX, mvmatrix );
  glGetDoublev( GL_PROJECTION_MATRIX, projmatrix );

  // Restore the GL state
  glPopMatrix();

  // Get the event state
  int x = Fl::event_x();
  int y = viewport[3] - Fl::event_y() - 1;

  // Vectors for the starting and ending points
  TracerData::Vec xStart, xEnd;

  // Perform the unprojection
  gluUnProject( (GLdouble) x, (GLdouble) y, -1.0,
    mvmatrix, projmatrix, viewport,
    &(xStart[0]), &(xStart[1]), &(xStart[2]) );

  gluUnProject( (GLdouble) x, (GLdouble) y, 1.0,
    mvmatrix, projmatrix, viewport,
    &(xEnd[0]), &(xEnd[1]), &(xEnd[2]) );

  // Find the point corresponding to the intersection
  vtkIdType iPoint;
  if(m_Data->PickPoint(xStart,xEnd,iPoint))
    {
    // Point has been found
    m_CurrentPoint = iPoint;
    }
  else
    {
    m_CurrentPoint = -1;
    }
}

int
TracerMainWindow
::handle(int event)
{
  // Get the button number
  int button = Fl::event_button();
  
  // accomodate 2-button mice: Ctrl-right => middle
  if (Fl::event_state(FL_CTRL) && button == FL_RIGHT_MOUSE)
    button = FL_MIDDLE_MOUSE;

  // Perform navigation tasks
  if(m_EditMode == TRACKBALL) 
    {
    if(event == FL_PUSH)
      {
      if(button == FL_LEFT_MOUSE)
        {
        m_TrackballMode = ROTATE;
        m_Trackball.StartRot( Fl::event_x(), Fl::event_y(), w(), h());
        }
      else if(button == FL_MIDDLE_MOUSE)
        {
        m_TrackballMode = PAN;
        m_Trackball.StartPan( Fl::event_x(), Fl::event_y() );
        }
      else if(button == FL_RIGHT_MOUSE)
        {
        m_TrackballMode = ZOOM;
        m_Trackball.StartZoom( Fl::event_y() );
        }
      return 1;
      } 
    else if(event == FL_RELEASE) 
      {
      if(m_TrackballMode == ROTATE)
        {
        m_Trackball.StopRot();
        }
      else if(m_TrackballMode == ZOOM)
        {
        m_Trackball.StopZoom();
        }
      else if(m_TrackballMode == PAN)
        {
        m_Trackball.StopPan();
        }
      return 1;
      }
    else if(event == FL_DRAG)
      {
      if(m_TrackballMode == ROTATE)
        {
        m_Trackball.TrackRot(Fl::event_x(), Fl::event_y(), w(), h());
        redraw();
        }
      else if(m_TrackballMode == ZOOM)
        {
        m_Trackball.TrackZoom( Fl::event_y() );
        redraw();
        }
      else if(m_TrackballMode == PAN)
        {
        m_Trackball.TrackPan(Fl::event_x(), Fl::event_y(), w(), h(), 1, 1);
        redraw();
        }
      return 1;
      }
    }
  else if(m_EditMode == TRACER)
    {
    if(event == FL_ENTER)
      {
      Fl::belowmouse(this);
      fl_cursor(FL_CURSOR_CROSS,FL_BLACK,FL_WHITE);
      return 1;
      }
    if(event == FL_LEAVE)
      {
      fl_cursor(FL_CURSOR_DEFAULT);
      return 1;
      }
    if(event == FL_MOVE)
      {
      // Find the point under cursor
      FindPointUnderCursor();
      redraw();

      // Handled
      return 1;
      }
    else if(event == FL_RELEASE)
      {
      // Again, find the point under the cursor
      FindPointUnderCursor();

      // Add the point as to the list
      m_Data->AddNewPoint(m_CurrentPoint);
      }
    }

  return 0;
}
