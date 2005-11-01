#include "MedialAtom.h"
#include <cmath>

MedialAtom::MedialAtom () {
	locus[0] = 0;
	locus[1] = 0;
	radius = 0;
	normal1[0] = 0;
	normal1[1] = 0;
	normal2[0] = 0;
	normal2[1] = 0;
}

MedialAtom::MedialAtom (const double x, const double y, const double r, const double xt, const double yt, const double rt) {
	set(x, y, r, xt, yt, rt);
}

void MedialAtom::set (const double x, const double y, const double r, const double xt, const double yt, const double rt) {
	locus[0] = x;
	locus[1] = y;
	radius = r;
	
	// compute the normals
	double l = xt * xt;
	l += yt * yt;
	double a = -rt / l;
	double b = l;
     b -=	rt * rt;
	if (b < 0) b = 0;
	b = sqrt(b);
	b /= l;
	double tmp1 = a * xt;
	double tmp2 = b * yt;
	normal1[0] = tmp1 - tmp2;
	normal2[0] = tmp1 + tmp2;
	tmp1 = a * yt; 
	tmp2 = b * xt;
	normal1[1] = tmp1 + tmp2;
	normal2[1] = tmp1 - tmp2;
}

void MedialAtom::phiToMedialAtoms (const double fx[], const double fy[], const double fxt[], const double fyt[], const double phi[], const int dim, MedialAtom m[]) {
	// the dimensions of the arrays fx, fy, fxt, fyt and m are (dim + 1)
	// the dimension of the array phi is (dim + 3)
	for (int i = 0; i < dim + 1; ++i) {
		double r = sqrt(phi[i + 1]);
		double rt = phi[i + 2] - phi[i];
		rt /= r;
		rt /= 4;
		rt *= dim;
		m[i].set(fx[i], fy[i], r, fxt[i], fyt[i], rt);
	}
	
	return;
}

ostream& operator<< (ostream& out, const MedialAtom& ma) {
	out << "locus = [" << ma.locus[0] << ", " << ma.locus[1] << "]" << endl;
	out << "radius = " << ma.radius << endl;
	out << "normal1 = [" << ma.normal1[0] << ", " << ma.normal1[1] << "]" << endl;
	out << "normal2 = [" << ma.normal2[0] << ", " << ma.normal2[1] << "]" << endl;
}

