//---------------------------------------------------------------------------
#include "plib.h"
#include <assert.h>

//---------------------------------------------------------------------------
void PList::reallocate() {
	allocated += dSize;

   void** temp = (void**)calloc(allocated,sizeof(void *));
   memcpy(temp,array,sizeof(void *)*getSize());
   free(array);
   array = temp;
	//memset(array,0,allocated*sizeof(PPersistent *));
}


PList::PList(int initSize,int incrSize,int magic) {
	iSize = initSize;
   dSize = incrSize;

   allocated = iSize;
   array = (void**)calloc(allocated,sizeof(void *));
   lastElt = -1;
}

PList::~PList() {
	free(array);
}

void PList::set(int index,void *p) {
	if(index >= allocated) {
   	allocated = index - (index % dSize);
      reallocate();
   }
 	if(index >= allocated) {
   	assert(0);
   }

   if(index > lastElt)
   	lastElt = index;
	array[index] = p;
}

void *PList::get(int index) {
	if(index <= lastElt && index >= 0)
   	return array[index];
   else
   	return NULL;
}

int PList::append(void *p) {
	int index = ++lastElt;
   set(index,p);
   return index;
}

int PList::getSize() {
	return (lastElt+1);
}

int PList::getAllocSize() {
	return allocated;
}

int PList::getIndex(void *p) {
	for(int i=0;i<getSize();i++) {
   	if(array[i]==p)
      	return i;
   }
   return -1;
}

void PList::clear() {
	free(array);
   allocated = iSize;
   array = (void**)calloc(allocated,sizeof(void *));
   lastElt = -1;
}

void PList::compact() {
	int cursor = 0;
	for(int i=0;i<getSize();i++) {
      if(array[i]) {
      	array[cursor] = array[i];
         if(i>cursor)
         	array[i] = NULL;
         cursor++;
      }
   }
   lastElt = cursor-1;
}




