#include "TracerMainWindow.h"
#include <FL/gl.h>

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
  glNewList(m_DisplayList,GL_COMPILE);

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

      // Specify normal.
      glNormal3d(nx, ny, nz);

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
::draw()
{
  if( !valid() )
    {
    // Set up the clear color
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glEnable( GL_DEPTH_TEST );

    // Set up the materials
    GLfloat light0Pos[4] = { 0.0, 0.0, 1.0, 0.0};
    GLfloat matAmb[4] = { 0.01, 0.01, 0.01, 1.00};
    GLfloat matDiff[4] = { 0.65, 0.65, 0.65, 1.00};
    GLfloat matSpec[4] = { 0.30, 0.30, 0.30, 1.00};
    GLfloat matShine = 10.0;
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, matAmb);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, matDiff);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, matSpec);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, matShine);
    glEnable(GL_COLOR_MATERIAL);

    // Setup Lighting
    glLightfv(GL_LIGHT0, GL_POSITION, light0Pos);
    glEnable(GL_LIGHT0);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glEnable(GL_LIGHTING);

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
    float zoom = m_Trackball.GetZoom();
    float x = 100 / zoom;
    glOrtho( -x , x, -x, x, -x, x ); 

    // Set up the model view matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    }

  // Create the display list if needed
  if(m_DisplayListDirty && m_Data->GetSourceMesh() != NULL)
    ComputeDisplayList();

  // Clear the screen 
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  // Set up the model view matrix
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();  
  glLoadIdentity();

  // Update the screen geometry
  glTranslatef( m_Trackball.GetPanX(), m_Trackball.GetPanY(), 0.0 );
  glMultMatrixf( m_Trackball.GetRot() );
  // glTranslatef( -m_CenterOfRotation[0], -m_CenterOfRotation[1], -m_CenterOfRotation[2] );

  // Draw the display list
  glCallList(m_DisplayList);
  
  glPopMatrix();
  glFlush();
}

int
  TracerMainWindow
::handle(int event)
{
}
