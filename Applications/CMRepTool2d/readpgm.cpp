#include <stdio.h>
#include "readpgm.h"
#include <string.h>

int readString(FILE *f,char buffer[],int maxlen) {
	int bp = -1;
	
	while(!feof(f) && bp < maxlen-1) {
		char c = fgetc(f);
		
		// \handle r\n pairs
		if(c=='\r')
			c = fgetc(f);

		if(c==' ' || c=='\t' || c=='\r' || c=='\n') {
			if(bp >= 0) {
				buffer[++bp] = 0;
				return 0;
			}
		}
		else if(c=='#') {
			while(c != '\n') c = fgetc(f);
			if(bp >= 0) {
				buffer[++bp] = 0;
				return 0;
			}
		}
		else 
			buffer[++bp] = c;
	}

	return -1;
}

int readpgm(const char *fname,cimage &image) {
	// Open a file
	FILE *f = fopen(fname,"rb");
	if(f==NULL) return PGM_OPEN_FAILED;
	
	// Read header
	char buffer[256];
	if(0 > readString(f,buffer,256))
		return PGM_READ_FAILED;

	// Check for valid file type
	if(buffer[0]!='P' || (buffer[1]!='5' && buffer[1]!='4'))
		return PGM_INVALID_TYPE;
	
	// Read width and height and intensity max
	int iw,ih,imax;
	if(0 > readString(f,buffer,256))
		return PGM_READ_FAILED;
	iw = atoi(buffer);
	if(0 > readString(f,buffer,256))
		return PGM_READ_FAILED;
	ih = atoi(buffer);
	if(0 > readString(f,buffer,256))
		return PGM_READ_FAILED;
	imax = atoi(buffer);

	// Create an image
	image = cimage_new(CGREY,iw,ih,1);
	
	// Read the image
	unsigned char *ibuf = new unsigned char[iw*ih];
	int rsize = (int) fread(ibuf,1,iw*ih,f); 
	if(iw*ih > rsize)
		return PGM_READ_FAILED;

	// Copy to our image
	int bi = -1;
	for(int y=0;y<ih;y++)
		for(int x=0;x<iw;x++)		 
			cimage_set(image,x,y,ibuf[++bi]);
	
	// Delete the image
	delete ibuf;

	return 0;
}

int writepgm(const char *fname,cimage &image) {
	// Open a file
	FILE *f = fopen(fname,"wb");
	if(f==NULL) return PGM_OPEN_FAILED;

   // Scale the image to bytes
   int w = cimage_xdim(image);
   int h = cimage_ydim(image);
   
   cimage_max_min(image);
   double imax = cimage_imax(image);
   int newMax = 255;
   if(imax > 1.0 && imax < 256.0) {
      newMax = (int)imax;
      imax = 255.0;
   }

   unsigned char *array = new unsigned char[w*h];
   int offset = -1;
   for(int y=0;y<h;y++) {
      for(int x=0;x<w;x++) {
         array[++offset] = (unsigned char)(cimage_get(image,x,y) * 255.0 / imax);
      }
   }
   
   // Write header
   char buffer[256];
   sprintf(buffer,"P5 %d %d %d\n",w,h,newMax);
   fwrite(buffer,1,strlen(buffer),f);

   // Write the data
   fwrite(array,1,w*h,f);

   delete array;
   fclose(f);

   return 0;
}

int readpgm(const char *fname,int *w,int *h,unsigned char **image) {
	// Open a file
	FILE *f = fopen(fname,"rb");
	if(f==NULL) return PGM_OPEN_FAILED;
	
	// Read header
	char buffer[256];
	if(0 > readString(f,buffer,256))
		return PGM_READ_FAILED;

	// Check for valid file type
	if(buffer[0]!='P' || buffer[1]!='5')
		return PGM_INVALID_TYPE;
	
	// Read width and height and intensity max
	if(0 > readString(f,buffer,256))
		return PGM_READ_FAILED;
	*w = atoi(buffer);
	if(0 > readString(f,buffer,256))
		return PGM_READ_FAILED;
	*h = atoi(buffer);
	if(0 > readString(f,buffer,256))
		return PGM_READ_FAILED;
	int imax = atoi(buffer);

	// Create an image
	*image = new unsigned char[(*w)*(*h)];
	
	// Read the image
	int rsize = (int) fread(*image,1,(*w)*(*h),f); 
	if((*w)*(*h) > rsize)
		return PGM_READ_FAILED;

   // Close file
   fclose(f);

	return 0;
}

int writepgm(const char *fname,int w,int h,unsigned char *image) {
	// Open a file
	FILE *f = fopen(fname,"wb");
	if(f==NULL) return PGM_OPEN_FAILED;

   // Write header
   char buffer[256];
   sprintf(buffer,"P5 %d %d 255\n",w,h);
   fwrite(buffer,1,strlen(buffer),f);

   // Write the data
   fwrite(image,1,w*h,f);

   fclose(f);

   return 0;
}

inline int clipImageValue(int point,int size) {
	return (point < 0) ? 0 : (point >= size) ? size-1 : point;
}
 
void lineDDA(Array2D<unsigned char> &image,int x0, int y0, int x1, int y1, unsigned char value)
{
	x0 = clipImageValue(x0,image.width());
	x1 = clipImageValue(x1,image.width());
	y0 = clipImageValue(y0,image.height());
	y1 = clipImageValue(y1,image.height());

	int dy = y1 - y0;
	int dx = x1 - x0;
	float t = (float) 0.5;                      // offset for rounding

	image(x0,y0) = value;;
	if (abs(dx) > abs(dy)) {          // slope < 1
		float m = (float) dy / (float) dx;      // compute slope
		t += y0;
		dx = (dx < 0) ? -1 : 1;
		m *= dx;
		while (x0 != x1) {
			x0 += dx;                           // step to next x value
			t += m;                             // add slope to y value
			image(x0,(int)t) = value;
		}
	} else {                                    // slope >= 1
		float m = (float) dx / (float) dy;      // compute slope
		t += x0;
		dy = (dy < 0) ? -1 : 1;
		m *= dy;
		while (y0 != y1) {
			y0 += dy;                           // step to next y value
			t += m;                             // add slope to x value
			image((int)t,y0) = value;
		}
	}
}

void floodFill(Array2D<unsigned char> &im,int x, int y, unsigned char value) {
	if(x>=0 && y>=0 && x<im.width() && y<im.height() && im(x,y)!=value) {
		im(x,y) = value;
		floodFill(im,x+1,y,value);
		floodFill(im,x-1,y,value);
		floodFill(im,x,y-1,value);
		floodFill(im,x,y+1,value);
	}
}
	
