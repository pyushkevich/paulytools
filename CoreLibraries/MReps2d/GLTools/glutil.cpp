#include "glutil.h"

#define min3(a,b,c) ((a < b && a < c) ? a : ((b < c) ? b : c))

#define max3(a,b,c) ((a > b && a > c) ? a : ((b > c) ? b : c))

void RGB_HSV(double r,double g,double b,
				 double &h,double &s,double &v) {
	double x = min3(r,g,b);
	v = max3(r,g,b);
	if(x==v) {
		h = 0;s = 0;
		return;
	}
	double f = (r==x) ? g-b : ((g==x) ? b-r : r-g);
	int i = (r==x)?3:((g==x)?5:1);
	h = (i-f/(v-x))/6;
	s = (v-x)/v;
	return;	
}

void HSV_RGB(double h,double s,double v,
				 double &r,double &g,double &b) {
	double m,n,f;
	int i;
	if(h==0 && s==0) {
		r = g = b = v;
		return;
	}
	h*=6;
	i = (int) h;
	f = h-i;
	
	if(!(i&1))
		f = 1-f;
	m = v * (1-s);
	n = v * (1-s*f);
	
	switch(i) {
	case 6 :
	case 0 :
		r = v;g = n;b = m;return;
	case 1 :
		r = n;g = v;b = m;return;
	case 2 :
		r = m;g = v;b = n;return;
	case 3 :
		r = m;g = n;b = v;return;
	case 4 :
		r = n;g = m;b = v;return;
	case 5 :
		r = v;g = m;b = n;return;
	}
}

const Gl_Color Gl_Color::white   = Gl_Color::rgb(1,1,1);
const Gl_Color Gl_Color::black   = Gl_Color::rgb(0,0,0);
const Gl_Color Gl_Color::red     = Gl_Color::rgb(1,0,0);
const Gl_Color Gl_Color::green   = Gl_Color::rgb(0,1,0);
const Gl_Color Gl_Color::blue    = Gl_Color::rgb(0,0,1);
const Gl_Color Gl_Color::yellow  = Gl_Color::rgb(1,1,0);
const Gl_Color Gl_Color::cyan    = Gl_Color::rgb(0,1,1);
const Gl_Color Gl_Color::magenta = Gl_Color::rgb(1,0,1);

void drawCircle(const Vector2D &x,double r) {
   static int dl = -1;
	double a;
   
   if(dl < 0) {
      dl = glGenLists(1);      
      glNewList(dl,GL_COMPILE);
      glBegin(GL_TRIANGLE_FAN);
      for(a=0;a<M_PI*2;a+=M_PI/5) {
         glVertex2d(cos(a),sin(a));
      }
      glEnd();
      glEndList();
   }
   
   
   glPushMatrix();
   glTranslated(x.x,x.y,0);
   glScaled(r,r,1);
   glCallList(dl);
   glPopMatrix();   
}

void drawOutlineCircle(const Vector2D &x,
							  double radius,
							  const Gl_Color &inside,
							  const Gl_Color &border)
{
	// Save the circle to a display list for future use
   static int dl = -1;
	double a;
   
   if(dl < 0) {
      dl = glGenLists(2);      
      glNewList(dl,GL_COMPILE);
      glBegin(GL_LINE_STRIP);
      for(a=0;a<=M_PI*2;a+=M_PI/5) {
         glVertex2d(cos(a),sin(a));
      }
      glEnd();
      glEndList();
      
      glNewList(dl+1,GL_COMPILE);
      glBegin(GL_TRIANGLE_FAN);
      for(a=0;a<M_PI*2;a+=M_PI/5) {
         glVertex2d(cos(a),sin(a));
      }
      glEnd();
      glEndList();
   };
   
   
   glPushMatrix();
   
	glTranslated(x.x,x.y,0);
   glScaled(radius,radius,1);

   
   border.apply();
	glCallList(dl);
   
   inside.apply();
	glCallList(dl+1);
   
   glPopMatrix();
   
}

void drawOutlineTriangle(const Vector2D &x,
								 double angle,
								 double radius,
								 const Gl_Color &inside,
								 const Gl_Color &border)
{
   static double s1 = 0.7;
   static double s2 = 0.4;
   static double x1 = -sqrt(3.0) * s1;
   static double x2 = -sqrt(3.0) * (s1-s2)/2;
   static double x3 = -sqrt(3.0) * (s1+s2)/2;
   
   glPushMatrix();
   glTranslated(x.x,x.y,0);
   glScaled(radius,radius,1);
   glRotated(-angle,0,0,1);
   
   glBegin(GL_TRIANGLES);
   
   border.apply();
   glVertex2d(0,0);
   glVertex2d(x1,-s1);
   glVertex2d(x1,s1);
   
   inside.apply();
   glVertex2d(x2,0);
   glVertex2d(x3,-s2);
   glVertex2d(x3,s2);
   
   glEnd();
   glPopMatrix();
}

void glDrawArrow(const Vector2D &T,const Vector2D &H,double maxArrowHeadSize) {
   static const double factor = sqrt(3.0) / 2.0;

   double len = H.distanceTo(T);
   double hypotenuse = maxArrowHeadSize < 0.4*len ? maxArrowHeadSize : 0.4*len;

   if(len == 0.0)
      return;

   Vector2D n = (H - T) / len;
   Vector2D n1 = n.getNormal();

   Vector2D C = T + n * (len - factor*hypotenuse);
   Vector2D A = C + n1 * (hypotenuse/2.0);
   Vector2D B = C - n1 * (hypotenuse/2.0);
   
   glBegin(GL_LINES);
   glVertex2d(T);
   glVertex2d(C);
   glEnd();

   glBegin(GL_POLYGON);
   glVertex2d(A);
   glVertex2d(B);
   glVertex2d(H);
   glEnd();

   glColor3d(0,0,0);
   glBegin(GL_LINE_STRIP);
   glVertex2d(A);
   glVertex2d(B);
   glVertex2d(H);
   glEnd();
}



