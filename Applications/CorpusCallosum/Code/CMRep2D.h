
#ifndef _CMRep2D_H
#define _CMRep2D_H

class MedialAtom;
class FunctionRep;

#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "matrix.h"

typedef float PixelType;
const int Dimension = 2;
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::LinearInterpolateImageFunction< ImageType, double > ImageInterpolator;

class CMRep2D {
	private:
	const int dim;
	double areaOfCMRep;
	double arcLenVar;
    	MedialAtom* medialAtoms;
	
	public:
	CMRep2D (const int _dim);
	~CMRep2D ();
	
	double  buildCMRep2D (const FunctionRep& fx, const FunctionRep& fy, const FunctionRep& frho, double *phi);
	void getBoundary(Vector &bx, Vector &by) ;
	double checkBoundaryFold () const;	
	double computeAreaOverlap (const ImageInterpolator *image, const double ratio, const int n, const int ne);
	double overlapAndCMRep (const ImageInterpolator *image, const double ratio, const int n, const int ne, const double dx); 
	double getAreaOfCMRep ();
	double getArcLenVar();
};

#endif

