
#ifndef _MedialAtom_H
#define _MedialAtom_H

#include <iostream>

using namespace std;

class MedialAtom {
	public:
	double locus[2];
	double radius;
	double normal1[2];
	double normal2[2];
	
	public:
	MedialAtom ();
	MedialAtom (const double x, const double y, const double r, const double xt, const double yt, const double rt);
	~MedialAtom () {}
	
	void set (const double x, const double y, const double r, const double xt, const double yt, const double rt);
	
	static void phiToMedialAtoms (const double fx[], const double fy[], const double fxt[], const double fyt[], const double phi[], const int dim, MedialAtom m[]);
	
	friend ostream& operator<< (ostream& out, const MedialAtom& ma);
};

#endif

