#ifndef _DANIELSON_H_
#define _DANIELSON_H_

int edt3ddan(char * input, int six, int siy, int siz, unsigned char metric,
			 short ** odx, short ** ody, short ** odz);

void edtSetVoxelSize(float vx,float vy,float vz);

#endif // _DANIELSON_H_

