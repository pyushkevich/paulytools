#include "mathematicaio.h"

#include <map>
#include <vector>
#include <string>
using namespace std;

void MathematicaIO::readTo(FILE *f, char c) {
  while(fgetc(f)!=c);
}

char MathematicaIO::readToSymbol(FILE *f) {
  char c;
  do {
    c = fgetc(f);
  } while(c == ' ' || c == '\t' || c == '\n' || c == '\r');
  return c;
}

char MathematicaIO::readDouble(FILE *f,double &d) {
  char c;
  string str;
  do {
    c = fgetc(f);
    if(c == ' ' || c == '\t' || c == '\n' || c == '\r')
      continue;
    str+=c;
  } while(c != ',' && c != '}');

  d = atof(str.c_str());
  return c;
}

int MathematicaIO::readMathematicaMatrix(FILE *f, Mat &A) {
  vector<vector<double> > vec;

  // Find open brace
  char c = readToSymbol(f);
  if(c!='{')
    return -1;

  // Read each vector or closed brace
  while(!feof(f)) {
    c = readToSymbol(f);
    if(c == '}') {
      // Matrix is finished
      A.setSize(vec.size(),vec.front().size());
      for(int r=0;r<A.rows();r++) {
	for(int c=0;c<A.columns();c++) {
	  A(r,c) = vec[r][c];
	}
      }
      break;
    }
    else if(c == '{') {
      // Read numbers
      vector<double> v;
      do {
	double d;
	c = readDouble(f,d);
	v.push_back(d);
      }	while(c == ',');
      vec.push_back(v);
    }
    else if(c == ',') continue;
    else {
      return -1;
    }
  }

  return 0;
}

int MathematicaIO::readMathematicaVector(FILE *f, Vec &A) {
  
  // Find open brace
  char c = readToSymbol(f);
  if(c!='{')
    return -1;

  // Read numbers
  vector<double> v;
  do {
    double d;
    c = readDouble(f,d);
    v.push_back(d);
  } while(c == ',');

  // Save vector
  A.setSize(v.size());
  for(int i=0;i<A.rows();i++)
    A(i) = v[i];
  
  return 0;
}
