#include "TracerMainWindow.h"
#include <FL/gl.h>
#include <FL/fl_ask.h>
#include <GL/glu.h>
#include <Fl/Fl.H>
#include <Fl/fl_draw.H>
#include <Fl/Fl_Color_Chooser.h>

#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>

#include "Trackball.h"

inline void glVertexNormal(double *x, double *n)
{
  glNormal3dv(n);
  glVertex3dv(x);
}

inline void glVertexNormalRGB(double *x, double *n, double *rgb)
{
  glColor3dv(rgb);
  glNormal3dv(n);
  glVertex3dv(x);
}

TracerMainWindow::
TracerMainWindow(int x, int y, int w, int h, const char *label) 
: Fl_Gl_Window(x,y,w,h,label) 
{
  m_Data = NULL;
  m_EdgeDisplayListDirty = false;
  m_FullMeshDisplayList = -1;
  m_NeighborhoodDisplayList = -1;
  m_EdgeDisplayList = -1;
  m_EditMode = TRACKBALL;
  m_TrackballMode = NONE;
  m_GLStateDirty = true;
  m_CurrentPoint = -1;
  m_CurrentMarkerCell = -1;
  m_EdgeDisplayMode = EDGE_DISPLAY_NONE;
  m_SurfaceDisplayMode = SURFACE_DISPLAY_ALL;
  m_CenterMesh = false;
  m_NeighborhoodSize = 10.0;
} 

TracerMainWindow
::~TracerMainWindow()
{
  if(m_Data)
    m_Data->RemoveTracerDataListener(this);
}

void
TracerMainWindow
::GLDrawStrippedPolyData(vtkPolyData *polys)
{
  // Get the triangle strip information.
  vtkCellArray *triStrips = polys->GetStrips();

  // Get the vertex information.
  vtkPoints *verts = polys->GetPoints();

  // Get the normal information.
  vtkDataArray *norms = polys->GetPointData()->GetNormals();

  // Data for getting consecutive points
  vtkIdType npts, *pts;
  for ( triStrips->InitTraversal(); triStrips->GetNextCell(npts,pts); ) 
    {
    glBegin( GL_TRIANGLE_STRIP );
    
    for (vtkIdType j = 0; j < npts; j++) 
      {
      // Some ugly code to ensure VTK version compatibility
      glNormal3dv(norms->GetTuple(pts[j]));
      glVertex3dv(verts->GetPoint(pts[j]));
      }
    
    glEnd();
    }
}

void 
TracerMainWindow
::ComputeFullMeshDisplayList()
{
  // Generate a display list number if needed
  if(m_FullMeshDisplayList < 0)
    m_FullMeshDisplayList = glGenLists(1);
  
  // Build display list
  glNewList(m_FullMeshDisplayList,GL_COMPILE_AND_EXECUTE);

  // Draw the base display list (non-markered)
  glColor3d(0.4,0.4,0.4);
  GLDrawStrippedPolyData(m_Data->GetDisplayMesh());

  // Draw each of the segmented triangle strips
  if(m_Data->GetNumberOfMarkers() > 0)
    {
    // Get an explicit pointer to TracerCurves
    const TracerCurves *curves = m_Data->GetCurves();

    // Get the marker ids
    TracerCurves::IdList lMarkers;
    curves->GetMarkerIdList(lMarkers);

    // Paint the corresponding meshes
    TracerCurves::IdList::iterator it = lMarkers.begin();
    while(it != lMarkers.end())
      {
      glColor3dv(curves->GetMarkerColor(*it).data_block());
      GLDrawStrippedPolyData(m_Data->GetMarkerMesh(*it));
      ++it;
      }
    }

  glEndList();

  m_FullMeshDisplayListDirty = false;
}

void 
TracerMainWindow
::ComputeNeighborhoodDisplayList()
{
  // Generate a display list number if needed
  if(m_NeighborhoodDisplayList < 0)
    m_NeighborhoodDisplayList = glGenLists(1);
  
  // Build display list
  glNewList(m_NeighborhoodDisplayList,GL_COMPILE_AND_EXECUTE);

  // Get the triangle strip information.
  vtkCellArray *polys = m_Data->GetInternalMesh()->GetPolys();

  // Get the vertex information.
  vtkPoints *verts = m_Data->GetInternalMesh()->GetPoints();

  // Get the normal information.
  vtkDataArray *norms = m_Data->GetInternalMesh()->GetPointData()->GetNormals();

  // Data for getting consecutive points
  vtkIdType ntris = 0;
  vtkIdType npts;
  vtkIdType *pts;
    
  glColor3d(0.4, 0.4, 0.4);
  glBegin(GL_TRIANGLES);
  
  vtkIdType iLastMarker = (vtkIdType) -2;
  for (unsigned int iCell = 0; iCell < polys->GetNumberOfCells(); iCell++)
    {
    // Get the triangle details
    m_Data->GetInternalMesh()->GetCellPoints(iCell, npts, pts);
    
    // Check if the triangle qualifies
    bool qual = true;
    for(vtkIdType iPoint = 0; qual && (iPoint < npts); iPoint++)
      if(m_Data->GetDistanceMapper()->GetVertexDistance(pts[iPoint]) 
        > m_NeighborhoodSize) qual = false;
    if(!qual) continue;

    // Check if we should color the triangle
    if(m_Data->IsMarkerSegmentationValid())
      {
      // Get the marker associated with the cell
      TracerCurves::IdType idMarker = m_Data->GetCellLabel(iCell);
      //cout << "Cell " << iCell << " maps to " << idMarker ;
      //cout << " distance " <<
      //  m_Data->GetVoronoiDiagram()->GetDiagram()->GetDistanceArray()[iCell] << endl;

      if(idMarker != iLastMarker)
        {
        if(idMarker != TracerData::NO_MARKER)
          glColor3dv(m_Data->GetCurves()->GetMarkerColor(idMarker).data_block());
        else
          glColor3d(0.4, 0.4, 0.4);

        iLastMarker = idMarker;
        }

      }

    // Display the qualified triangle
    for(vtkIdType iPoint = 0; iPoint < npts; iPoint++)
      {
      glNormal3dv(norms->GetTuple3(pts[iPoint]));
      glVertex3dv(verts->GetPoint(pts[iPoint]));
      }
    }

  glEnd();
  
  glEndList();

  // Clear the dirty flag
  m_NeighborhoodDisplayListDirty = false;
}

void 
TracerMainWindow
::GLDrawSphere(double *x, double r)
{
  glPushMatrix();
  glTranslated(x[0],x[1],x[2]);

  GLUquadricObj *sphere = gluNewQuadric();
  gluSphere(sphere, r, 10, 10);
  gluDeleteQuadric(sphere);

  glPopMatrix();
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
  SetGLColorHSV(hue, 0.5, 0.7);
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
  SetGLColorHSV(hue, 0.5, 0.7);
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

  // Special adjustments are needed when there is no source point
  if(!m_Data->IsPathSourceSet())
    {
    if(m_SurfaceDisplayMode == SURFACE_DISPLAY_NEIGHBORHOOD)
      xActualMode = EDGE_DISPLAY_NONE;
    else if(xActualMode == EDGE_DISPLAY_DISTANCE) 
      xActualMode = EDGE_DISPLAY_PLAIN;
    }

  // A similar flag for neighborhood display
  SurfaceDisplayMode xActualSurfaceMode = m_SurfaceDisplayMode;
  if(!m_Data->IsPathSourceSet())
    xActualSurfaceMode = SURFACE_DISPLAY_ALL;

  // Only do something if the mode is set
  if(xActualMode != EDGE_DISPLAY_NONE)
    {
    // Draw the edges
    glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT | GL_COLOR_BUFFER_BIT);

    // Set the default drawing color
    glColor3f(0.6f, 0.6f, 0.6f);
    
    // Start line drawing
    glBegin(GL_LINES);

    // Go through all the edges
    for(vtkIdType p = 0; p < dm->GetNumberOfVertices(); p++)
      {
      // If we are in neighborhood surface display mode and this point
      // is too far away, reject it
      if(xActualSurfaceMode == SURFACE_DISPLAY_NEIGHBORHOOD &&
        m_Data->GetDistanceMapper()->GetVertexDistance(p) > m_NeighborhoodSize)
        continue;

      // Get the number of edges
      unsigned int nEdges = dm->GetVertexNumberOfEdges(p);

      // Draw each edge
      for(unsigned int j = 0; j < nEdges; j++)
        {
        // Get the adjacent point
        vtkIdType q = dm->GetVertexEdge(p, j);
        
        // Only process each edge once
        if( p < q )
          {
          // If we are in neighborhood surface display mode and this point
          // is too far away, reject it
          if(xActualSurfaceMode == SURFACE_DISPLAY_NEIGHBORHOOD &&
            m_Data->GetDistanceMapper()->GetVertexDistance(q) > m_NeighborhoodSize)
            continue;

          vnl_vector_fixed<double,3> x1, x2, n1, n2;
          m_Data->GetPointCoordinateAndNormal(p, x1, n1);
          m_Data->GetPointCoordinateAndNormal(q, x2, n2);

          // If in plain mode, don't set any color
          if(xActualMode == EDGE_DISPLAY_PLAIN)
            {
            glVertexNormal(x1.data_block(), n1.data_block());
            glVertexNormal(x2.data_block(), n2.data_block());
            }
          else if(m_EdgeDisplayMode == EDGE_DISPLAY_LENGTH)
            {
            // Get the weight of the edge
            float xWeight = dm->GetVertexEdgeWeight(p, j);
            
            // Compute and set the color based on length
            SetGLEdgeColorFromWeight(xWeight);

            // Draw the vertices
            glVertexNormal(x1.data_block(), n1.data_block());
            glVertexNormal(x2.data_block(), n2.data_block());
            }
          else if(m_EdgeDisplayMode == EDGE_DISPLAY_DISTANCE)
            {
            SetGLEdgeColorFromDistance(dm->GetVertexDistance(p));
            glVertexNormal(x1.data_block(), n1.data_block());
            
            SetGLEdgeColorFromDistance(dm->GetVertexDistance(q));
            glVertexNormal(x2.data_block(), n2.data_block());
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
::DrawMesh()
{
  // Draw the display list
  glColor3d(0.4,0.5,0.3);
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);

  // If the mode is to display the entire surface, use the display mesh
  if(m_SurfaceDisplayMode == SURFACE_DISPLAY_ALL || !m_Data->IsPathSourceSet())
    {
    // Create the display list if needed
    if(m_FullMeshDisplayListDirty) ComputeFullMeshDisplayList();
    else glCallList(m_FullMeshDisplayList);
    }
  else 
    {
    // Create the display list if needed
    if(m_NeighborhoodDisplayListDirty) ComputeNeighborhoodDisplayList();
    else glCallList(m_NeighborhoodDisplayList);
    }

  // Create the edge display list if needed
  if(m_EdgeDisplayListDirty) ComputeEdgeDisplayList();
  else glCallList(m_EdgeDisplayList);
}

void
TracerMainWindow
::DrawCurves()
{
  // Set the attributes for line drawing
  glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT | GL_COLOR_BUFFER_BIT);
  glEnable(GL_BLEND);
  glEnable(GL_LINE_SMOOTH);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);    
  glLineWidth(3);

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
      Vec xPoint, xNormal;
      m_Data->GetPointCoordinateAndNormal(*pit, xPoint, xNormal);
      glVertexNormal(xPoint.data_block(), xNormal.data_block());
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


    glColor3d(0.8,0.6,0.2);
    glBegin(GL_LINE_STRIP);

    MeshCurve::const_iterator it = lPoints.begin();
    while(it != lPoints.end())
      {
      Vec xPoint, xNormal;
      m_Data->GetPointCoordinateAndNormal(*it, xPoint, xNormal);
      glVertexNormal(xPoint.data_block(), xNormal.data_block());
      it++;
      }

    glEnd();
    }

  // Restore the attributes
  glPopAttrib();

  // Determine the appropriate radius of the little balls
  double radius = 0.5;

  // Draw control points (balls representing points)
  for(IdList::iterator it = lCurves.begin(); it!=lCurves.end(); it++)
    {
    // Get the control points in the curve
    const IdList &lControls = 
      m_Data->GetCurves()->GetCurveControls(*it);

    // Draw each of the control points
    glColor3d(0.2,0.6,0.8);

    IdList::const_iterator itControl = lControls.begin();
    while(itControl != lControls.end())
      {
      // Don't draw the current control
      if(*itControl != m_Data->GetCurrentControlPoint())
        {
        // Get the coordinate of the control point
        Vec xPoint = m_Data->GetCurves()->GetControlPointPosition(*itControl);

        // Display the point currently under the cursor as a little sphere
        GLDrawSphere(xPoint.data_block(), radius);
        }
      ++itControl;
      }
    }

  // Draw the current path source
  if(m_Data->IsPathSourceSet())
    {
    // Get the coordinate of the current point
    Vec xPoint;
    m_Data->GetPointCoordinate(m_Data->GetPathSource(), xPoint);

    // Display the point currently under the cursor as a little sphere
    glColor3d(0.8,0.8,0.2);
    GLDrawSphere(xPoint.data_block(), radius);
    }

  // Draw the current point
  if(m_CurrentPoint != -1 && m_EditMode == TRACER)
    {
    // Get the coordinate of the current point
    Vec xPoint;
    m_Data->GetPointCoordinate(m_CurrentPoint,xPoint);

    // Display the point currently under the cursor as a little sphere
    glColor3d(0.8,0.2,0.2);
    GLDrawSphere(xPoint.data_block(), radius);
    }
}

void
TracerMainWindow
::GLDrawMarker(vtkIdType iCell, const Vec &color)
{
  // Get the corners of the marker
  Vec x, y, z, c, n, f[4], p;
  m_Data->GetCellGeometry(iCell, x, y, z, n, c);

  // Compute the direction of the flag
  Vec v(1.0,0.0,0.0);
  if(v == n)
    v = Vec(0.0,1.0,0.0);
  p = v - dot_product(v,n) * n;
  p.normalize();

  // Compute the corners of the flag
  f[0] = c + 10.0 * n;
  f[1] = c + 15.0 * n;
  f[2] = c + 15.0 * n + 8.0 * p;
  f[3] = c + 10.0 * n + 8.0 * p;
  
  // Save the state
  glPushAttrib(GL_LIGHTING_BIT);
  glDisable(GL_LIGHTING);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);    
  glLineWidth(1.0);

  // Draw the field of the flag
  glColor3dv(color.data_block());
  glBegin(GL_QUADS);
  glVertex3dv(f[0].data_block());
  glVertex3dv(f[1].data_block());
  glVertex3dv(f[2].data_block());
  glVertex3dv(f[3].data_block());
  glEnd();

  // Draw the stick and outline of the flag

  // Set the color of the flag based on selected state
  if(iCell == 
     m_Data->GetCurves()->GetMarkerFace(
       m_Data->GetCurrentMarker()))
    {
    glColor3d(1.0,1.0,1.0);
    }
  else
    {
    glColor3d(0.6,0.6,0.6);
    }

  // Draw the flag
  glBegin(GL_LINE_STRIP);
  glVertex3dv(c.data_block());
  glVertex3dv(f[0].data_block());
  glVertex3dv(f[1].data_block());
  glVertex3dv(f[2].data_block());
  glVertex3dv(f[3].data_block());
  glVertex3dv(f[0].data_block());
  glEnd();

  // Restore state
  glPopAttrib();
}

void 
TracerMainWindow
::DrawMarkers()
{
  // Get the curves object
  const TracerCurves *curves = m_Data->GetCurves();

  // Get the list of marker ID's
  TracerCurves::IdList lMarkers;
  curves->GetMarkerIdList(lMarkers);
  
  // Markers are displayed as little lolly-pops sticking out of the surface
  TracerCurves::IdList::iterator itMarker = lMarkers.begin();
  while(itMarker != lMarkers.end())
    {
    // Get the color of the marker
    Vec xColor = curves->GetMarkerColor(*itMarker);
    vtkIdType iFace = curves->GetMarkerFace(*itMarker);

    // Draw the marker
    GLDrawMarker(iFace, xColor);

    // On to the next
    ++itMarker;
    }

  // Draw a white marker at the current position of the mouse
  if(m_CurrentMarkerCell != -1 && m_EditMode == MARKER)
    {
    Vec xWhite(1.0, 1.0, 1.0);
    GLDrawMarker(m_CurrentMarkerCell, xWhite);
    }
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

    // Draw the mesh
    DrawMesh();

    // Draw markers on the surface
    DrawMarkers();

    // Draw the curves on the surface
    DrawCurves();
    
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

  // Determine which point should be the center
  Vec xCenter;
  if(m_CenterMesh && m_Data->IsPathSourceSet())
    m_Data->GetPointCoordinate(m_Data->GetPathSource(), xCenter);
  else
    xCenter.set(m_Data->GetDisplayMesh()->GetCenter());
    
  // Translate to the center
  glTranslated(-xCenter[0],-xCenter[1],-xCenter[2]);
}

void
TracerMainWindow
::ComputeClickRay(Vec &xStart, Vec &xEnd)
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

  // Perform the unprojection
  gluUnProject( (GLdouble) x, (GLdouble) y, -1.0,
    mvmatrix, projmatrix, viewport,
    &(xStart[0]), &(xStart[1]), &(xStart[2]) );

  gluUnProject( (GLdouble) x, (GLdouble) y, 1.0,
    mvmatrix, projmatrix, viewport,
    &(xEnd[0]), &(xEnd[1]), &(xEnd[2]) );
}

void
TracerMainWindow
::FindCellUnderCursor()
{
  // Vectors for the starting and ending points
  Vec xStart, xEnd;

  // Get the ray corrensponding to the click
  ComputeClickRay(xStart,xEnd);

  // Find the corresponding cell
  vtkIdType iCell;
  if(m_Data->GetDistanceMapper()->PickCell(xStart, xEnd, iCell))
    {
    // A cell has been located
    m_CurrentMarkerCell = iCell;
    }
  else
    {
    m_CurrentMarkerCell = -1;
    }
}

bool
TracerMainWindow
::CheckCell(vtkIdType iCell)
{
  if(m_Data->IsPathSourceSet() && 
    m_SurfaceDisplayMode == SURFACE_DISPLAY_NEIGHBORHOOD)
    {
    // We are in neighborhood display mode. Check that the distance
    // from all corners of the cell is close enough
    return
      m_Data->GetCellDistanceToPathSource(iCell) < m_NeighborhoodSize;
    }
  else
    return true;
}

void 
TracerMainWindow
::FindPointUnderCursor()
{
  // Vectors for the starting and ending points
  Vec xStart, xEnd;

  // Get the ray corrensponding to the click
  ComputeClickRay(xStart,xEnd);

  // Find the point corresponding to the intersection
  vtkIdType iPoint;
  if(m_Data->GetDistanceMapper()->PickPoint(xStart,xEnd,iPoint,this))
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

  /** Tracer Mode */
  else if(m_EditMode == MARKER)
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
      FindCellUnderCursor();
      redraw();

      // Handled
      return 1;
      }
    else if(event == FL_RELEASE)
      {
      // Again, find the point under the cursor
      FindCellUnderCursor();

      if(m_CurrentMarkerCell != -1)
        {
        // Make up a default marker name
        ostringstream sout;
        sout << "Marker " << m_CurrentMarkerCell;
        
        // Prompt for the name of the marker
        const char *sName = 
          fl_input("Please enter the name of the marker",sout.str().c_str());
        if(sName) 
          {
          // Pick a color for the marker
          double r, g, b;
          double h = rand() * 6.0 / RAND_MAX, s = 0.25 + rand() * 0.5 / RAND_MAX, v = 0.75;
          Fl_Color_Chooser::hsv2rgb(h,s,v,r,g,b);
          if(fl_color_chooser("Please select the color for the marker",r,g,b))
            {
            m_Data->AddMarker(m_CurrentMarkerCell, sName, r, g, b);
            }
          }

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
  m_NeighborhoodDisplayListDirty = true;
  m_FullMeshDisplayListDirty = true;
  m_EdgeDisplayListDirty = true;
  m_CurrentPoint = -1;
  m_CurrentMarkerCell = -1;
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
  // m_CurrentPoint = -1;

  // The edge drawing is no longer clean if drawing accumulated distances
  m_NeighborhoodDisplayListDirty = true;
  if(m_EdgeDisplayMode == EDGE_DISPLAY_DISTANCE || 
    m_SurfaceDisplayMode == SURFACE_DISPLAY_NEIGHBORHOOD)
    {
    m_EdgeDisplayListDirty = true;
    }
  
  // Special feature: if the current mode is to center on the path source, 
  // and the path source changes, remove the pan in the trackball
  if(m_CenterMesh && m_Data->IsPathSourceSet())
    {
    m_Trackball.ResetPan();
    }

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

void
TracerMainWindow
::OnMarkerListChange(TracerDataEvent *evt)
{
  m_CurrentMarkerCell = -1;
  redraw();
}

void
TracerMainWindow
::OnFocusMarkerChange(TracerDataEvent *evt)
{
  m_CurrentMarkerCell = -1;
  redraw();
}

void
TracerMainWindow
::OnSegmentationChange(TracerDataEvent *evt)
{
  m_FullMeshDisplayListDirty = true;
  m_NeighborhoodDisplayListDirty = true;
  redraw();
}
