#include "TracerMainWindow.h"
#include <FL/gl.h>
#include <GL/glu.h>
#include <Fl/Fl.H>
#include <Fl/fl_draw.H>
#include <Fl/Fl_Color_Chooser.h>

#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>

#include "Trackball.h"

TracerMainWindow::
TracerMainWindow(int x, int y, int w, int h, const char *label) 
: Fl_Gl_Window(x,y,w,h,label) 
{
  m_Data = NULL;
  m_DisplayListDirty = false;
  m_EdgeDisplayListDirty = false;
  m_DisplayList = -1;
  m_EdgeDisplayList = -1;
  m_EditMode = TRACKBALL;
  m_TrackballMode = NONE;
  m_GLStateDirty = true;
  m_CurrentPoint = -1;
  m_EdgeDisplayMode = EDGE_DISPLAY_NONE;
  m_SurfaceDisplayMode = SURFACE_DISPLAY_ALL;
} 

TracerMainWindow
::~TracerMainWindow()
{
  if(m_Data)
    m_Data->RemoveTracerDataListener(this);
}

void 
TracerMainWindow
::ComputeDisplayList()
{
  // Generate a display list number if needed
  if(m_DisplayList < 0)
    m_DisplayList = glGenLists(1);

  // Get the triangle strip information.
  vtkCellArray *triStrips = m_Data->GetDisplayMesh()->GetStrips();
    
  // Get the vertex information.
  vtkPoints *verts = m_Data->GetDisplayMesh()->GetPoints();

  // Get the normal information.
  vtkDataArray *norms = m_Data->GetDisplayMesh()->GetPointData()->GetNormals();

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
::SetGLColorHSV(double xHue, double xSaturation, double xValue)
{
  double r, g, b;
  Fl_Color_Chooser::hsv2rgb(xHue * 6.0, xSaturation, xValue, r, g, b);
  glColor3d(r,g,b);
}

void
TracerMainWindow
::SetGLEdgeColorFromDistance(double xDistance)
{
  // Use distance to source to paint the vertex in
  // a given color. For now, use the hue map
  double hue = fmod(xDistance * 0.01, 1.0);
  
  // Set the HSV color
  SetGLColorHSV(hue, 0.5, 1.0);
}

void
TracerMainWindow
::SetGLEdgeColorFromWeight(double xWeight)
{
  // Use distance to source to paint the vertex in
  // a given color. For now, use the hue map
  double hue = xWeight * 0.5;
  if(hue > 0.67) hue = 0.67;
  
  // Set the HSV color
  SetGLColorHSV(hue, 0.5, 1.0);
}

void
TracerMainWindow
::ComputeEdgeDisplayList()
{
  // Generate a display list number if needed
  if(m_EdgeDisplayList < 0)
    m_EdgeDisplayList = glGenLists(1);

  // Get the data needed to draw the edges
  const VTKMeshShortestDistance *dm = m_Data->GetDistanceMapper();

  // Build display list
  glNewList(m_EdgeDisplayList,GL_COMPILE_AND_EXECUTE);

  // Asjust the mode based on the current state
  EdgeDisplayMode xActualMode = m_EdgeDisplayMode;
  if(xActualMode == EDGE_DISPLAY_DISTANCE && !m_Data->IsPathSourceSet())
    xActualMode = EDGE_DISPLAY_PLAIN;

  // Only do something if the mode is set
  if(xActualMode != EDGE_DISPLAY_NONE)
    {

    // Draw the edges
    glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT | GL_COLOR_BUFFER_BIT);
    glDisable(GL_LIGHTING);

    // Set the default drawing color
    glColor3f(0.2f, 0.8f, 0.6f);
    
    // Start line drawing
    glBegin(GL_LINES);

    // Go through all the edges
    for(vtkIdType p = 0; p < dm->GetNumberOfVertices(); p++)
      {
      unsigned int nEdges = dm->GetVertexNumberOfEdges(p);
      for(unsigned int j = 0; j < nEdges; j++)
        {
        // Get the adjacent point
        vtkIdType q = dm->GetVertexEdge(p, j);
        
        // Only process each edge once
        if( p < q )
          {
          vnl_vector_fixed<double,3> x1, x2;
          m_Data->GetPointCoordinate(p, x1);
          m_Data->GetPointCoordinate(q, x2);

          // If in plain mode, don't set any color
          if(xActualMode == EDGE_DISPLAY_PLAIN)
            {
            glVertex3dv(x1.data_block());
            glVertex3dv(x2.data_block());
            }
          else if(m_EdgeDisplayMode == EDGE_DISPLAY_LENGTH)
            {
            // Get the weight of the edge
            float xWeight = dm->GetVertexEdgeWeight(p, j);
            
            // Compute and set the color based on length
            SetGLEdgeColorFromWeight(xWeight);

            // Draw the vertices
            glVertex3dv(x1.data_block());
            glVertex3dv(x2.data_block());
            }
          else if(m_EdgeDisplayMode == EDGE_DISPLAY_DISTANCE)
            {
            SetGLEdgeColorFromDistance(dm->GetVertexDistance(p));
            glVertex3dv(x1.data_block());
            
            SetGLEdgeColorFromDistance(dm->GetVertexDistance(q));
            glVertex3dv(x2.data_block());
            }
          }
        }
      }

    // Restore state
    glEnd();
    glPopAttrib();
    }

  // Close the display list
  glEndList();

  // The display list is 'clean'
  m_EdgeDisplayListDirty = false;
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

  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1.0,1.0);

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
  if(m_Data->GetDisplayMesh() != NULL)
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
    
    // Create the edge display list if needed
    if(m_EdgeDisplayListDirty) ComputeEdgeDisplayList();
    else glCallList(m_EdgeDisplayList);
    
    // Set the attributes for line drawing
    glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT | GL_COLOR_BUFFER_BIT);
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(2.5);

    // Get a list of curve id's
    typedef list<TracerCurves::IdType> IdList;
    typedef TracerCurves::MeshCurve MeshCurve;
    typedef const TracerCurves::MeshCurve ConstMeshCurve;
    
    IdList lCurves;
    m_Data->GetCurves()->GetCurveIdList(lCurves);
    
    // Display each of the curves
    for(IdList::iterator it = lCurves.begin(); it!=lCurves.end(); it++)
      {
      // Choose a color for the curve
      if(*it == m_Data->GetCurrentCurve()) 
        glColor3d(1,0.5,0.5);
      else 
        glColor3d(0.5,0,0);

      // Get the curve as a list of vertex IDs
      ConstMeshCurve &lPoints = m_Data->GetCurves()->GetCurveVertices(*it);

      // Draw the curve
      glBegin(GL_LINE_STRIP);
      
      ConstMeshCurve::const_iterator pit = lPoints.begin(); 
      while(pit != lPoints.end())
        {
        TracerData::Vec xPoint;
        m_Data->GetPointCoordinate(*pit, xPoint);
        glVertex3dv(xPoint.data_block());
        ++pit;
        }
      glEnd();
      }

    // Display the path to the current point
    if(m_EditMode == TRACER && m_Data->IsPathSourceSet() &&
      m_CurrentPoint != -1 && m_CurrentPoint != m_Data->GetPathSource())
      {
      // Get the path between the source point and the current point
      MeshCurve lPoints;
      m_Data->GetPathToPoint(m_CurrentPoint, lPoints);

      // Since the point list does not include the starting and 
      // ending points, we must include them ourselves
      lPoints.push_back(m_CurrentPoint);
      lPoints.push_front(m_Data->GetPathSource());

      // Need at least two points to draw
      MeshCurve::const_iterator it1 = lPoints.begin();
      MeshCurve::const_iterator it2 = it1; it2++;

      glBegin(GL_LINES);
      
      while(it2 != lPoints.end())
        {
        TracerData::Vec xPoint;

        // Get the distance measurement
        double xDistance = m_Data->GetMeshEdgeWeight(*it1, *it2);

        // Use xDistance to set the color
        double xColorBase = 1 * xDistance;
        glColor3d(1 - xColorBase, 1, xColorBase);

        // Plot the first points
        m_Data->GetPointCoordinate(*it1, xPoint);
        glVertex3dv(xPoint.data_block());

        m_Data->GetPointCoordinate(*it2, xPoint);
        glVertex3dv(xPoint.data_block());

        // Increase the iterators
        it1++; it2++;
        }

      glEnd();
      }

    // Restore the attributes
    glPopAttrib();

    // Draw control points (balls representing points)
    for(IdList::iterator it = lCurves.begin(); it!=lCurves.end(); it++)
      {
      // Get the control points in the curve
      const IdList &lControls = 
        m_Data->GetCurves()->GetCurveControls(*it);

      // Draw each of the control points
      IdList::const_iterator itControl = lControls.begin();
      while(itControl != lControls.end())
        {
        TracerData::Vec xPoint;

        // Display the point currently under the cursor as a little sphere
        glColor3d(0.2,0.6,0.8);
        
        glPushMatrix();
        m_Data->GetPointCoordinate(
          m_Data->GetCurves()->GetControlPointVertex(*itControl),xPoint);
        glTranslated(xPoint(0),xPoint(1),xPoint(2));
        
        GLUquadricObj *sphere = gluNewQuadric();
        gluSphere(sphere, 0.5, 10, 10);
        gluDeleteQuadric(sphere);

        glPopMatrix();

        ++itControl;
        }
      }
    
    if(m_CurrentPoint != -1 && m_EditMode == TRACER)
      {
      TracerData::Vec xPoint;
      
      // Display the point currently under the cursor as a little sphere
      glColor3d(0.8,0.2,0.2);
      glPushMatrix();
      m_Data->GetPointCoordinate(m_CurrentPoint,xPoint);
      
      glTranslated(xPoint(0),xPoint(1),xPoint(2));
      GLUquadricObj *sphere = gluNewQuadric();
      gluSphere(sphere, 0.5, 10, 10);
      gluDeleteQuadric(sphere);

      glPopMatrix();
      }

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
    2.0 / (m_Trackball.GetZoom() * m_Data->GetDisplayMesh()->GetLength());
  glScaled(xScale,xScale,xScale);

  // Now, translate for the center of the mesh
  double *xCenter = m_Data->GetDisplayMesh()->GetCenter();
  glTranslated(-xCenter[0],-xCenter[1],-xCenter[2]);
}

void 
TracerMainWindow
::FindPointUnderCursor()
{
  // Find the point currently under the cursor
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

  /** Tracer Mode */
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

      if(m_CurrentPoint != -1)
        {
        // Add the point as to the list
        m_Data->AddNewPoint(m_CurrentPoint);

        // The edge drawing is no longer clean if drawing accumulated distances
        if(m_EdgeDisplayMode == EDGE_DISPLAY_DISTANCE)
          m_EdgeDisplayListDirty = true;

        // Handled
        return 1;
        }
      }
    }

  return 0;
}
  
void 
TracerMainWindow
::OnMeshChange(TracerDataEvent *evt)
{
  m_DisplayListDirty = true;
  m_EdgeDisplayListDirty = true;
  m_CurrentPoint = -1;
  redraw();
}
  
void 
TracerMainWindow
::OnFocusCurveChange(TracerDataEvent *evt)
{
  m_CurrentPoint = -1;
  redraw();
}
  
void 
TracerMainWindow
::OnFocusPointChange(TracerDataEvent *evt)
{
  // Get the currently selected point from the data
  m_CurrentPoint = -1;

  // The edge drawing is no longer clean if drawing accumulated distances
  if(m_EdgeDisplayMode == EDGE_DISPLAY_DISTANCE)
    m_EdgeDisplayListDirty = true;

  // If surface display is in neighborhood mode, dirty it
  if(m_SurfaceDisplayMode == SURFACE_DISPLAY_NEIGHBORHOOD)
    m_DisplayListDirty = true;

  // Redraw the window
  redraw();
}

void 
TracerMainWindow
::OnFocusCurveDataChange(TracerDataEvent *evt)
{
  redraw();
}
  
void 
TracerMainWindow
::OnCurveListChange(TracerDataEvent *evt)
{
  redraw();
}
  
void 
TracerMainWindow
::OnEdgeWeightsUpdate(TracerDataEvent *evt)
{
  m_EdgeDisplayListDirty = true;
  redraw();
}
