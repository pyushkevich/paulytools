#ifndef _IMPLICIT_H_
#define _IMPLICIT_H_

/**********************************************************************
 * DEFINITIONS FOR IMPLICIT SURFACE TRACKING									 *
 **********************************************************************/
#define TET	0  /* use tetrahedral decomposition */
#define NOTET	1  /* no tetrahedral decomposition  */

typedef struct {		   /* a three-dimensional point */
	double x, y, z;		   /* its coordinates */
} IMP_POINT;

typedef struct {		   /* test the function for a signed value */
	IMP_POINT p;			   /* location of test */
	double value;		   /* function value at p */
	int ok;			   /* if value is of correct sign */
} IMP_TEST;

typedef struct {		   /* surface vertex */
	IMP_POINT position, normal;	   /* position and surface normal */
} IMP_VERTEX;

typedef struct {	   /* list of vertices in polygonization */
	int count, max;		   /* # vertices, max # allowed */
	IMP_VERTEX *ptr;		   /* dynamically allocated */
} IMP_VERTICES;

extern "C" {
char *polygonize(	double (*f)(double,double,double),
						double size,
						int bounds,
						double  x,
						double  y,
						double  z,
						int (*triproc)(int,int,int,IMP_VERTICES),
						int mode);
};

#endif // _IMPLICIT_H_
