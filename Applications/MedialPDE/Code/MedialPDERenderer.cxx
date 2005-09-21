#include "MedialPDESolver.h"
#include "MedialPDERenderer.h"
#include "OptimizationTerms.h"

void glDrawQuadStripElements(unsigned short width,unsigned short height) 
{
  // Allocate the index array
  int size = width*2;
  unsigned short *index = new unsigned short[size];

  unsigned short iStart = 0,iStart2 = width;

  for (int j=0;j<height-1;j++)
    {
    int iIndex = 0;
    for (int i=0;i<width;i++)
      {
      index[iIndex++] = iStart++;
      index[iIndex++] = iStart2++;
      }
    glDrawElements(GL_QUAD_STRIP,size,GL_UNSIGNED_SHORT,index);
    }

  delete index;
}

void glDrawQuadElements(MedialAtomGrid *grid)
{
  // Begin drawing quads
  glBegin(GL_QUADS);
  
  // Iterate over the quads
  MedialQuadIterator *itQuad = grid->NewQuadIterator();
  while(!itQuad->IsAtEnd())
    {
    glArrayElement(itQuad->GetAtomIndex(0, 0));
    glArrayElement(itQuad->GetAtomIndex(0, 1));
    glArrayElement(itQuad->GetAtomIndex(1, 1));
    glArrayElement(itQuad->GetAtomIndex(1, 0));
    ++(*itQuad);
    }
  delete itQuad;

  // End quads
  glEnd();
}

void glDrawWireframeElements(unsigned short width, unsigned short height)
{
  // Horizontal elements as lines
  unsigned int size = height + width;
  unsigned short *index = new unsigned short[size];

  // Draw the vertical lines
  unsigned short iStart = 0;
  for(int j = 0; j < height; j++)
    {
    int iIndex = 0;
    for(int i = 0; i < width; i++)
      { index[iIndex++] = iStart++; }
    glDrawElements(GL_LINE_STRIP, width, GL_UNSIGNED_SHORT, index);
    }

  // Draw the horizontal lines
  for(int i = 0; i < width; i++)
    {
    iStart = i;
    int iIndex = 0;
    for(int j = 0; j < height; j++)
      { index[iIndex++] = iStart; iStart+=width; }
    glDrawElements(GL_LINE_STRIP, height, GL_UNSIGNED_SHORT, index);
    }
}

PDESplineRenderer
::PDESplineRenderer(MedialPDESolver *solver)
{
  this->solver = solver;
  // matMedial = new GLMaterial(GL_FRONT_AND_BACK, 
  //  GLColor(0.1), GLColor(0.4, 0.4, 0.8), GLColor(0.15), 64); 
  matMedial = new GLMaterial(GL_FRONT_AND_BACK, 
    GLColor(0.1), GLColor(0.4, 0.4, 0.8), GLColor(0.0, 0)); 

  matBoundary = new GLMaterial(GL_FRONT_AND_BACK, 
    GLColor(0.1), GLColor(0.4));
}

void PDESplineRenderer::DrawInternalPoints( size_t nCuts )
{
  // Generate the internal points using solutiondata
  SolutionData S(solver);
  S.UpdateInternalWeights( nCuts );

  // Pass the points as vertex pointers
  // glVertexPointer(3, GL_DOUBLE, sizeof(SMLVec3d), 
  //  S.xInternalPoints[0].data_block());

  // Start point rendering
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glDisable(GL_LIGHTING);
  glBegin(GL_POINTS);

  // Draw all the points
  MedialInternalPointIterator *it = 
    S.xAtomGrid->NewInternalPointIterator( nCuts );
  
  for( ; !it->IsAtEnd(); ++(*it))
    {
    double d = it->GetRelativeDistanceToMedialAxis();
    glColor3d(0.8, 0.2, 1 - d);
    glVertex3dv( S.xInternalPoints[it->GetIndex()].data_block() );
    }
    //glArrayElement(it->GetIndex());

  delete it;
  
  glEnd();

  glBegin(GL_LINES);
  MedialProfileIntervalIterator *itp = S.xAtomGrid->NewProfileIntervalIterator(nCuts);
  for(; !itp->IsAtEnd(); ++(*itp))
    {
    glColor3d(1,1,1);
    glVertex3dv(S.xInternalPoints[itp->GetInnerPointIndex()].data_block());
    glColor3d(0,1,0);
    glVertex3dv(S.xInternalPoints[itp->GetOuterPointIndex()].data_block());
    }
  delete itp;
  glEnd();



  
  glPopAttrib();
}

void
PDESplineRenderer
::build()
{
  // Set the display attributes
  glPushAttrib(GL_LIGHTING_BIT);
  glEnable(GL_LIGHTING);

  // Center the object
  glPushMatrix();
  // glTranslated(-0.5, -0.5, -0.0);

  // Enable vector arrays
  glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT); 
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);

  // Get the numtber of points
  unsigned int m = solver->GetNumberOfUPoints();
  unsigned int n = solver->GetNumberOfVPoints();

  // Initialize the medial colors
  matMedial->apply();

  // Supply the vertex list (medial surface)
  MedialAtom *mp = solver->GetAtomArray();
  glVertexPointer(3, GL_DOUBLE, sizeof(MedialAtom), mp->X.data_block());
  glNormalPointer(GL_DOUBLE, sizeof(MedialAtom), mp->X.data_block());

  // Build the quad array
  glColor3d(1, 0, 0);
  // glDrawWireframeElements(m, n);
  // glDrawQuadElements(solver->GetAtomGrid());
  
  // DrawInternalPoints( 5 );

  // Display the boundary
  matBoundary->apply();

  // Supply the first boundary array
  glVertexPointer(3, GL_DOUBLE, sizeof(MedialAtom), mp->xBnd[0].X.data_block());
  glNormalPointer(GL_DOUBLE, sizeof(MedialAtom), mp->xBnd[0].N.data_block());

  // Draw the wireframe
  // glDrawQuadElements(solver->GetAtomGrid());
  glDrawWireframeElements(m, n);

  // Supply the second boundary array
  glVertexPointer(3, GL_DOUBLE, sizeof(MedialAtom), mp->xBnd[1].X.data_block());
  glNormalPointer(GL_DOUBLE, sizeof(MedialAtom), mp->xBnd[1].N.data_block());

  // Draw the wireframe
  glDrawWireframeElements(m, n);

  // Junk
   
  glDisable(GL_LIGHTING);
  glBegin(GL_POINTS);
  for(size_t i = 0; i < solver->GetAtomGrid()->GetNumberOfAtoms(); i++)
    {
    glColor3d(1,0,0);
    glVertex3d(20 * mp[i].u - 10, 20 * mp[i].v - 10, mp[i].F);
    }
  glEnd();
  
  glBegin(GL_POINTS);
  for(size_t i = 0; i < solver->GetAtomGrid()->GetNumberOfAtoms(); i++)
    {
    if(mp[i].xBadness > 0)
        glColor3d(0,1,0);
    else 
        glColor3d(1,1,0);
    glVertex3d(20 * mp[i].u - 10, 20 * mp[i].v - 10, mp[i].xBadness);
    }
  glEnd();
  
  
  // Restore the client state
  glPopClientAttrib();

  // Restore the GL state
  glPopMatrix();
  glPopAttrib();
}
