#ifndef _MATLAB_H_
#define _MATLAB_H_

#include <ostream>
#include <engine.h>

using namespace std;

class Matlab {
private:
	Engine *matEngine;
	enum {BUFFER = 0x10000};
	char *buffer;
public:
	Matlab(const char *host = "\0");
	~Matlab();

	bool connected();
	void command(const char *command,ostream &out);

	void exportMatrix(const char *name,float *start,int rows,int cols,int stride);
};


#endif