/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date$
  Version:   $Revision$
  Copyright (c) 2003 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#ifndef __Trackball_h_
#define __Trackball_h_

#include <FL/gl.h>

/**
 * \class Trackball
 * \brief Virtual trackball for the 3D window
 *
 * \sa Window3D
 */
class Trackball  
{
private:
  GLboolean m_TrackingMotion;
  float m_Angle;
  float m_Axis[3];
  float m_LastPosition[3];
  GLfloat m_RotationMatrix[4][4];
  GLfloat m_Zoom, m_OldZoom;
  GLfloat m_PanX, m_PanY;
  GLfloat m_OldPanX, m_OldPanY;
  
  void PToV( int x, int y, int width, int height, float v[3] );

public:
  Trackball();
  Trackball( const Trackball& T );
  ~Trackball();

  void Reset();
  void ResetPan();
  void StartPan( int x, int y );
  void StopPan();
  void TrackPan( int x, int y, int w, int h, float ratew, float rateh );
  void StartZoom( int y );
  void StopZoom();
  void TrackZoom( int y );
  void StartRot( int x, int y, int w, int h );
  void StopRot();
  void TrackRot( int x, int y, int w, int h );
  inline GLfloat *GetRot() { return( (GLfloat*) m_RotationMatrix ); };
  inline GLfloat GetZoom() { return( m_Zoom ); };
  inline GLfloat GetPanX() { return( m_PanX ); };
  inline GLfloat GetPanY() { return( m_PanY ); };

  Trackball& operator= ( const Trackball& T );
};

#endif // __Trackball_h_

/*
 *$Log$
 *Revision 1.1  2005/09/21 14:07:10  pauly2
 *Moved and built BrainTracer
 *
 *Revision 1.2  2004/04/02 15:18:02  pauly2
 *Separate display lists for global and local meshes
 *
 *Revision 1.1  2004/03/05 23:10:42  pauly2
 *Added braintrace app
 *
 *Revision 1.4  2003/10/02 14:55:53  pauly
 *ENH: Development during the September code freeze
 *
 *Revision 1.1  2003/09/11 13:51:15  pauly
 *FIX: Enabled loading of images with different orientations
 *ENH: Implemented image save and load operations
 *
 *Revision 1.3  2003/08/27 14:03:24  pauly
 *FIX: Made sure that -Wall option in gcc generates 0 warnings.
 *FIX: Removed 'comment within comment' problem in the cvs log.
 *
 *Revision 1.2  2003/08/27 04:57:47  pauly
 *FIX: A large number of bugs has been fixed for 1.4 release
 *
 *Revision 1.1  2003/07/12 04:46:51  pauly
 *Initial checkin of the SNAP application into the InsightApplications tree.
 *
 *Revision 1.3  2003/07/12 01:34:18  pauly
 *More final changes before ITK checkin
 *
 *Revision 1.2  2003/07/11 23:25:33  pauly
 **** empty log message ***
 *
 *Revision 1.1  2003/03/07 19:29:48  pauly
 *Initial checkin
 *
 *Revision 1.1.1.1  2002/12/10 01:35:36  pauly
 *Started the project repository
 *
 *
 *Revision 1.2  2002/03/08 14:06:32  moon
 *Added Header and Log tags to all files
 **/
