// #include <engine.h>
#include "matlab.h"

void Matlab::command(const char *command,ostream &out) {
	if(matEngine) {
		engEvalString(matEngine,command);
		out << buffer;
	}
	else {
		out << "MATLAB not available" << endl;
	}
}

Matlab::Matlab(const char *host) {	
	matEngine = engOpen(host);
	if(matEngine) {
		buffer = new char[BUFFER];
		engOutputBuffer(matEngine,buffer,BUFFER);
	}
}

Matlab::~Matlab() {	
	if(matEngine) {
		engClose(matEngine);
		delete buffer;
	}
}

bool Matlab::connected() {
	return (matEngine!=NULL);	
}

void Matlab::exportMatrix(const char *name,float *start,int rows,int cols,int stride) {
	mxArray *T = mxCreateDoubleMatrix(rows,cols,mxREAL);
	mxSetName(T,name);

	double *ptr = mxGetPr(T);
	for(int r=0;r<rows;r++) {
		for(int c=0;c<cols;c++) {
			*ptr = *start;
			ptr++;
			start+=stride;
		}
	}

	engPutArray(matEngine,T);
	mxDestroyArray(T);
}

