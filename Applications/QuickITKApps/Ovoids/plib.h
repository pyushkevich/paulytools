//---------------------------------------------------------------------------
#ifndef plibH
#define plibH

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>

// Comment this line if you don't want to read/write matrices
// with PStream
// #include <matrix/src/matrix.h>

class PPersistent;
class PStream;
class PHashTable;
class PList;

// FACTORY - add a function that will create instances of your Persistent
// classes based in magicNumber passed in.
const static int PMAGIC_PLIST = 0;
const static int PMAGIC_PHASHTABLE = 1;
const static int PSYSTEM_MAGIC_SIZE = 2;

typedef PPersistent* (*factoryFunction)(int magicNumber);
struct PFactoryFunction;

class PFactory {
protected:
	static PFactoryFunction *list, *p;
   static void initFactory();
public:
   static int addFunction(factoryFunction f,int mgcStart,int mgcEnd);

	// System factory - creates system classes
   static PPersistent *system(int magicNumber);

   // Create an instance for magic number
	static PPersistent* use(int magicNumber);
};

class PList {
protected:
	int iSize,dSize;
   int allocated;
   int lastElt;

   void **array;
	void reallocate();

public:
	PList(int initSize=10,int incrSize=10,int magic=PMAGIC_PLIST);
   ~PList();

   void set(int index,void *p);
   void *get(int index);
   int append(void *p);

	int getSize();
   int getAllocSize();
   int getIndex(void *p);

   void clear();
   void compact();
};

#endif
