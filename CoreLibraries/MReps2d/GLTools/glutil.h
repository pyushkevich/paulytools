#ifndef _GL_UTIL_H_
#define _GL_UTIL_H_

#include <vector2d.h>
#include <Fl/gl.h>

void RGB_HSV(double r,double g,double b,
				 double &h,double &s,double &v);

void HSV_RGB(double h,double s,double v,
				 double &r,double &g,double &b);

// A color class - makes it easy to save and store GL colors.
class Gl_Color {
protected:
public:
	double r,g,b,a;

	// Some predefined colors
	const static Gl_Color white,black,red,green,blue,yellow,cyan,magenta;

	// Create a new color
	static Gl_Color rgb(double red,double green,double blue,double alpha=1.0) {
		Gl_Color c;
		c.r = red;
		c.g = green;
		c.b = blue;
		c.a = alpha;
		return c;
	}

	// Returns a color object that represents current color
	static Gl_Color current() {
		GLdouble clr[4];
		glGetDoublev(GL_CURRENT_COLOR,clr);
		return rgb(clr[0],clr[1],clr[2],clr[3]);
	}

	// Create a new HSV color
	static Gl_Color hsv(double hue,double saturation,double value,double alpha=1.0) {
		Gl_Color c;
		c.a = alpha;
		HSV_RGB(hue,saturation,value,c.r,c.g,c.b);		
		return c;
	}

	double h() {
		double H,S,V;
		RGB_HSV(r,g,b,H,S,V);
		return H;
	}

	double s() {
		double H,S,V;
		RGB_HSV(r,g,b,H,S,V);
		return S;
	}

	double v() {
		double H,S,V;
		RGB_HSV(r,g,b,H,S,V);
		return V;
	}

	// This method causes the color to be applied (using GL)
	void apply() const {
		glColor4d(r,g,b,a);
	}

	// This method returns the same color with an alpha change
	Gl_Color alpha(double newAlpha) const {
		Gl_Color c;
		c.r = r;c.g = g;c.b = b;c.a = newAlpha;
		return c;
	}
};

// Super helpful little cheat
inline void glVertex2d(const Vector2D &x) {
	glVertex2d(x.x,x.y);
}

// Draw a circle (using current color)
void drawCircle(const Vector2D &x,double r);

// Draw an outline circle using inside color and outline color
void drawOutlineCircle(const Vector2D &x,
							  double radius,
							  const Gl_Color &inside,
							  const Gl_Color &border = Gl_Color::black);

// Draw an outline equilateral triangle
void drawOutlineTriangle(const Vector2D &x,
								 double angle,
								 double radius,
								 const Gl_Color &inside,
								 const Gl_Color &border = Gl_Color::black);

// Draw an arrow from one point to another, (using current color)
void glDrawArrow(const Vector2D &T,const Vector2D &H,double maxArrowHeadSize);

#endif
