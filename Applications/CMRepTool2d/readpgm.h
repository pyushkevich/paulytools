#ifndef READ_PGM_H
#define READ_PGM_H

#include "array2d.h"


#define PGM_OPEN_FAILED -1
#define PGM_READ_FAILED -2
#define PGM_INVALID_TYPE -3
#define PGM_OK 0

// Read PGM into a cimage
int readpgm(const char *,cimage&);

// Write a cimage to PGM
int writepgm(const char *,cimage&);

// Read into a byte array - you have to dealloc the array yourself
int readpgm(const char *fname,int *w,int *h,unsigned char **image);

// Write byte array to a PGM file
int writepgm(const char *fname,int w,int h,unsigned char *image);

void lineDDA(Array2D<unsigned char> &image,int x0, int y0, int x1, int y1, unsigned char value);

void floodFill(Array2D<unsigned char> &im,int x, int y, unsigned char value);

#endif
