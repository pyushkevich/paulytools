#ifndef __mathematicaio_h_
#define __mathematicaio_h_

#include <cstdio>
#include "LPWrapper.h"

class MathematicaIO {
public:
	static void readTo(FILE *f, char c);
	static char readToSymbol(FILE *f);
	static char readDouble(FILE *f,double &d);
	static int readMathematicaMatrix(FILE *f, Mat &A);

	static int readMathematicaVector(FILE *f, Vec &A);
};









#endif


