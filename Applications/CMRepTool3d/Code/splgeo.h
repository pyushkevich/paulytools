#ifndef _SPLINE_GEO_H_
#define _SPLINE_GEO_H_

#include "mspline.h"
#include "geodesic.h"

class SplineGeoMan : public GeodesicManifold {
private:
	MSpline *spline;
	SplineGridDefinition *grid;
public:
	SplineGeoMan(MSpline *spline, SplineGridDefinition *grid) : GeodesicManifold(2,3) {
		this->spline = spline;
		this->grid = grid;
	}

	void eval(double *x,double *y,double *J,double *H);
};

#endif