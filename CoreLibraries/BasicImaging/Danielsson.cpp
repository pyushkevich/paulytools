#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>

// Voxel sizes
float edtVoxX2 = 1.0f;
float edtVoxY2 = 1.0f;
float edtVoxZ2 = 1.0f;

#define SQ(x) sqtab[x]

// #define DIS(x,y,z) (SQ(x)+SQ(y)+SQ(z))
#define DIS(x,y,z) (SQ(x)*edtVoxX2+SQ(y)*edtVoxY2+SQ(z)*edtVoxZ2)

#define MAX(x,y) (((x) > (y)) ? (x) : (y))


/* Extension of  Danielsson algorithm  to 3D */
/* CGIP 14, 227-248			     */
/* From Martin Styner			     */
/* Modified by Sarang Joshi to handle char input */

void edtSetVoxelSize(float vx,float vy,float vz) {
	edtVoxX2 = vx * vx;
	edtVoxY2 = vy * vy;
	edtVoxZ2 = vz * vz;
}



int edt3ddan(char * input, int six, int siy, int siz, unsigned char metric, short ** odx, short ** ody, short ** odz)
{
  int x=0,y=0,z=0,index=0,indexc=0, dd=0;
  int *sqtab = NULL;
  int *sq = NULL;
  int maxsiz = 0;
  short *dx=NULL, *dy=NULL, *dz=NULL;
  int stop=0, upy=0;
  int part=0;
  int dist=0, dst=0;

  maxsiz = MAX(six,siy);
  maxsiz = MAX(maxsiz,siz);
  maxsiz = 4 * maxsiz;
  if(!(sq=sqtab=(int *)malloc((4*maxsiz+1)*sizeof(int)))) {
    return(-1);
  }
  sqtab = &sqtab[2*maxsiz];
  for(index=(-2*maxsiz);index<=(2*maxsiz);index++) {
     sqtab[index] = index * index;
  }
  if(!(dx=(short *)malloc((six*siy*siz)*sizeof(short)))) {
    free(sq);
    return(-1);
  }
  if(!(dy=(short *)malloc((six*siy*siz)*sizeof(short)))) {
    free(dx);
    free(sq);
    return(-1);
  }
  if(!(dz=(short *)malloc((six*siy*siz)*sizeof(short)))) {
    free(dy);
    free(dx);
    free(sq);
    return(-1);
  }
  dd = six * siy;
  for(x=0;x<six;x++) {
    for(y=0;y<siy;y++) {
      for(z=0;z<siz;z++) {
        index = x + six * y + dd * z;
        if(input[index]) {
          dx[index] = dy[index] = dz[index] = 0;
        }
        else {
          dx[index] = dy[index] = dz[index] = maxsiz;
        }
      }
    }
  }
  for(z=1;z<siz;z++) {
    part = (int) (20.0 + 25.0 * ((float) z - 1.0) / ((float) siz - 1.0));
    upy = siy - 1;
    stop = six - 1;
    for(y=0;y<siy;y++) {
      index = z * dd + y * six;
      for(x=0;x<six;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index - dd;
        if((dst=DIS(dx[indexc],dy[indexc],dz[indexc]-1)) < dist) {
          dx[index] = dx[indexc];
          dy[index] = dy[indexc];
          dz[index] = dz[indexc] - 1;
          dist = dst;
        }
        if(metric != 6) {
          if(y > 0) {
            indexc = index - six - dd;
            if((dst=DIS(dx[indexc],dy[indexc]-1,dz[indexc]-1)) < dist) {
              dx[index] = dx[indexc];
              dy[index] = dy[indexc] - 1;
              dz[index] = dz[indexc] - 1;
              dist = dst;
            }
          }
          if(y < upy) {
            indexc = index + six - dd;
            if((dst=DIS(dx[indexc],dy[indexc]+1,dz[indexc]-1)) < dist) {
              dx[index] = dx[indexc];
              dy[index] = dy[indexc] + 1;
              dz[index] = dz[indexc] - 1;
              dist = dst;
            }
          }
          if(x > 0) {
            indexc = index - 1 - dd;
            if((dst=DIS(dx[indexc]-1,dy[indexc],dz[indexc]-1)) < dist) {
              dx[index] = dx[indexc] - 1;
              dy[index] = dy[indexc];
              dz[index] = dz[indexc] - 1;
              dist = dst;
            }
          }
          if(x < stop) {
            indexc = index + 1 - dd;
            if((dst=DIS(dx[indexc]+1,dy[indexc],dz[indexc]-1)) < dist) {
              dx[index] = dx[indexc] + 1;
              dy[index] = dy[indexc];
              dz[index] = dz[indexc] - 1;
              dist = dst;
            }
          }
          if(metric == 26) {
            if(y > 0) {
              if(x > 0) {
                indexc = index - dd - six - 1;
                if((dst=DIS(dx[indexc]-1,dy[indexc]-1,dz[indexc]-1)) < dist) {
                  dx[index] = dx[indexc] - 1;
                  dy[index] = dy[indexc] - 1;
                  dz[index] = dz[indexc] - 1;
                  dist = dst;
                }
              }
              if(x < stop) {
                indexc = index - dd - six + 1;    
                if((dst=DIS(dx[indexc]+1,dy[indexc]-1,dz[indexc]-1)) < dist) {
                  dx[index] = dx[indexc] + 1;
                  dy[index] = dy[indexc] - 1;
                  dz[index] = dz[indexc] - 1;
                  dist = dst;
                }
              }
            }
            if(y < upy) {
              if(x > 0) {
                indexc = index - dd + six - 1;
                if((dst=DIS(dx[indexc]-1,dy[indexc]+1,dz[indexc]-1)) < dist) {
                  dx[index] = dx[indexc] - 1;
                  dy[index] = dy[indexc] + 1;
                  dz[index] = dz[indexc] - 1;
                  dist = dst;
                }
              }
              if(x < stop) {
                indexc = index - dd + six + 1;
                if((dst=DIS(dx[indexc]+1,dy[indexc]+1,dz[indexc]-1)) < dist) {
                  dx[index] = dx[indexc] + 1;
                  dy[index] = dy[indexc] + 1;
                  dz[index] = dz[indexc] - 1;
                  dist = dst;
                }
              }
            }
          }
        }
      }
    }
    for(y=1;y<siy;y++) {
      index = z * dd + y * six;
      for(x=0;x<six;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index - six;
        if((dst=DIS(dx[indexc],dy[indexc]-1,dz[indexc])) < dist) {
          dx[index] = dx[indexc];
          dy[index] = dy[indexc] - 1;
          dz[index] = dz[indexc];
          dist = dst;
        }
        if(metric != 6) {
          if(x > 0) {
            indexc = index - six - 1;
            if((dst=DIS(dx[indexc]-1,dy[indexc]-1,dz[indexc])) < dist) {
              dx[index] = dx[indexc] - 1;
              dy[index] = dy[indexc] - 1;
              dz[index] = dz[indexc];
              dist = dst;
            }
          }
          if(x < stop) {
            indexc = index - six + 1;
            if((dst=DIS(dx[indexc]+1,dy[indexc]-1,dz[indexc])) < dist) {
              dx[index] = dx[indexc] + 1;
              dy[index] = dy[indexc] - 1;
              dz[index] = dz[indexc];
              dist = dst;
            }
          }
        }
      }
      index = z * dd + y * six + 1;
      for(x=1;x<six;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index - 1;
        if((dst=DIS(dx[indexc]-1,dy[indexc],dz[indexc])) < dist) {
          dx[index] = dx[indexc] - 1;
          dy[index] = dy[indexc];
          dz[index] = dz[indexc];
          dist = dst;
        }
      }
      index = z * dd + (y + 1) * six - 2;
      for(x=(six-2);x>=0;x--,index--) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index + 1;
        if((dst=DIS(dx[indexc]+1,dy[indexc],dz[indexc])) < dist) {
          dx[index] = dx[indexc] + 1;
          dy[index] = dy[indexc];
          dz[index] = dz[indexc];
          dist = dst;
        }
      }
    }
    for(y=(siy-2);y>=0;y--) {
      index = z * dd + y * six;
      for(x=0;x<six;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index + six;
        if((dst=DIS(dx[indexc],dy[indexc]+1,dz[indexc])) < dist) {
          dx[index] = dx[indexc];
          dy[index] = dy[indexc] + 1;
          dz[index] = dz[indexc];
          dist = dst;
        }
        if(metric != 6) {
          if(x > 0) {
            indexc = index + six - 1;
            if((dst=DIS(dx[indexc]-1,dy[indexc]+1,dz[indexc])) < dist) {
              dx[index] = dx[indexc] - 1;
              dy[index] = dy[indexc] + 1;
              dz[index] = dz[indexc];
              dist = dst;
            }
          }
          if(x < stop) {
            indexc = index + six + 1;
            if((dst=DIS(dx[indexc]+1,dy[indexc]+1,dz[indexc])) < dist) {
              dx[index] = dx[indexc] + 1;
              dy[index] = dy[indexc] + 1;
              dz[index] = dz[indexc];
              dist = dst;
            }
          }
        }
      }
      index = z * dd + y * six + 1;
      for(x=1;x<six;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index - 1;
        if((dst=DIS(dx[indexc]-1,dy[indexc],dz[indexc])) < dist) {
          dx[index] = dx[indexc] - 1;
          dy[index] = dy[indexc];
          dz[index] = dz[indexc];
          dist = dst;
        }
      }
      index = z * dd + (y + 1) * six - 2;
      for(x=(six-2);x>=0;x--,index--) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index + 1;
        if((dst=DIS(dx[indexc]+1,dy[indexc],dz[indexc])) < dist) {
          dx[index] = dx[indexc] + 1;
          dy[index] = dy[indexc];
          dz[index] = dz[indexc];
          dist = dst;
        }
      }
    }
  }
  for(z=siz-2;z>=0;z--) {
    part = (int) (45.0 + 25.0 * ((float) (siz - 2 - z) / ((float) siz - 1.0)));
    upy = siy - 1;
    stop = six - 1;
    for(y=0;y<siy;y++) {
      index = z * dd + y * six;
      for(x=0;x<six;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index + dd;
        if((dst=DIS(dx[indexc],dy[indexc],dz[indexc]+1)) < dist) {
          dx[index] = dx[indexc];
          dy[index] = dy[indexc];
          dz[index] = dz[indexc] + 1;
          dist = dst;
        }
        if(metric != 6) {
          if(y > 0) {
            indexc = index - six + dd;
            if((dst=DIS(dx[indexc],dy[indexc]-1,dz[indexc]+1)) < dist) {
              dx[index] = dx[indexc];
              dy[index] = dy[indexc] - 1;
              dz[index] = dz[indexc] + 1;
              dist = dst;
            }
          }
          if(y < upy) {
            indexc = index + six + dd;
            if((dst=DIS(dx[indexc],dy[indexc]+1,dz[indexc]+1)) < dist) {
              dx[index] = dx[indexc];
              dy[index] = dy[indexc] + 1;
              dz[index] = dz[indexc] + 1;
              dist = dst;
            }
          }
          if(x > 0) {
            indexc = index - 1 + dd;
            if((dst=DIS(dx[indexc]-1,dy[indexc],dz[indexc]+1)) < dist) {
              dx[index] = dx[indexc] - 1;
              dy[index] = dy[indexc];
              dz[index] = dz[indexc] + 1;
              dist = dst;
            }
          }
          if(x < stop) {
            indexc = index + 1 + dd;
            if((dst=DIS(dx[indexc]+1,dy[indexc],dz[indexc]+1)) < dist) {
              dx[index] = dx[indexc] + 1;
              dy[index] = dy[indexc];
              dz[index] = dz[indexc] + 1;
              dist = dst;
            }
          }
          if(metric == 26) {
            if(y > 0) {
              if(x > 0) {
                indexc = index + dd - six - 1;
                if((dst=DIS(dx[indexc]-1,dy[indexc]-1,dz[indexc]+1)) < dist) {
                  dx[index] = dx[indexc] - 1;
                  dy[index] = dy[indexc] - 1;
                  dz[index] = dz[indexc] + 1;
                  dist = dst;
                }
              }
              if(x < stop) {
                indexc = index + dd - six + 1;    
                if((dst=DIS(dx[indexc]+1,dy[indexc]-1,dz[indexc]+1)) < dist) {
                  dx[index] = dx[indexc] + 1;
                  dy[index] = dy[indexc] - 1;
                  dz[index] = dz[indexc] + 1;
                  dist = dst;
                }
              }
            }
            if(y < upy) {
              if(x > 0) {
                indexc = index + dd + six - 1;
                if((dst=DIS(dx[indexc]-1,dy[indexc]+1,dz[indexc]+1)) < dist) {
                  dx[index] = dx[indexc] - 1;
                  dy[index] = dy[indexc] + 1;
                  dz[index] = dz[indexc] + 1;
                  dist = dst;
                }
              }
              if(x < stop) {
                indexc = index + dd + six + 1;
                if((dst=DIS(dx[indexc]+1,dy[indexc]+1,dz[indexc]+1)) < dist) {
                  dx[index] = dx[indexc] + 1;
                  dy[index] = dy[indexc] + 1;
                  dz[index] = dz[indexc] + 1;
                  dist = dst;
                }
              }
            }
          }
        }
      }
    }
    for(y=1;y<siy;y++) {
      index = z * dd + y * six;
      for(x=0;x<six;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index - six;
        if((dst=DIS(dx[indexc],dy[indexc]-1,dz[indexc])) < dist) {
          dx[index] = dx[indexc];
          dy[index] = dy[indexc] - 1;
          dz[index] = dz[indexc];
          dist = dst;
        }
        if(metric != 6) {
          if(x > 0) {
            indexc = index - six - 1;
            if((dst=DIS(dx[indexc]-1,dy[indexc]-1,dz[indexc])) < dist) {
              dx[index] = dx[indexc] - 1;
              dy[index] = dy[indexc] - 1;
              dz[index] = dz[indexc];
              dist = dst;
            }
          }
          if(x < stop) {
            indexc = index - six + 1;
            if((dst=DIS(dx[indexc]+1,dy[indexc]-1,dz[indexc])) < dist) {
              dx[index] = dx[indexc] + 1;
              dy[index] = dy[indexc] - 1;
              dz[index] = dz[indexc];
              dist = dst;
            }
          }
        }
      }
      index = z * dd + y * six + 1;
      for(x=1;x<six;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index - 1;
        if((dst=DIS(dx[indexc]-1,dy[indexc],dz[indexc])) < dist) {
          dx[index] = dx[indexc] - 1;
          dy[index] = dy[indexc];
          dz[index] = dz[indexc];
          dist = dst;
        }
      }
      index = z * dd + (y + 1) * six - 2;
      for(x=(six-2);x>=0;x--,index--) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index + 1;
        if((dst=DIS(dx[indexc]+1,dy[indexc],dz[indexc])) < dist) {
          dx[index] = dx[indexc] + 1;
          dy[index] = dy[indexc];
          dz[index] = dz[indexc];
          dist = dst;
        }
      }
    }
    for(y=(siy-2);y>=0;y--) {
      index = z * dd + y * six;
      for(x=0;x<six;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index + six;
        if((dst=DIS(dx[indexc],dy[indexc]+1,dz[indexc])) < dist) {
          dx[index] = dx[indexc];
          dy[index] = dy[indexc] + 1;
          dz[index] = dz[indexc];
          dist = dst;
        }
        if(metric != 6) {
          if(x > 0) {
            indexc = index + six - 1;
            if((dst=DIS(dx[indexc]-1,dy[indexc]+1,dz[indexc])) < dist) {
              dx[index] = dx[indexc] - 1;
              dy[index] = dy[indexc] + 1;
              dz[index] = dz[indexc];
              dist = dst;
            }
          }
          if(x < stop) {
            indexc = index + six + 1;
            if((dst=DIS(dx[indexc]+1,dy[indexc]+1,dz[indexc])) < dist) {
              dx[index] = dx[indexc] + 1;
              dy[index] = dy[indexc] + 1;
              dz[index] = dz[indexc];
              dist = dst;
            }
          }
        }
      }
      index = z * dd + y * six + 1;
      for(x=1;x<six;x++,index++) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index - 1;
        if((dst=DIS(dx[indexc]-1,dy[indexc],dz[indexc])) < dist) {
          dx[index] = dx[indexc] - 1;
          dy[index] = dy[indexc];
          dz[index] = dz[indexc];
          dist = dst;
        }
      }
      index = z * dd + (y + 1) * six - 2;
      for(x=(six-2);x>=0;x--,index--) {
        dist = DIS(dx[index],dy[index],dz[index]);
        indexc = index + 1;
        if((dst=DIS(dx[indexc]+1,dy[indexc],dz[indexc])) < dist) {
          dx[index] = dx[indexc] + 1;
          dy[index] = dy[indexc];
          dz[index] = dz[indexc];
          dist = dst;
        }
      }
    }
  }
  for(z=0;z<siz;z++) {
    for(y=0;y<siy;y++) {
      index = z * dd + y * six;
      for(x=0;x<six;x++,index++) {
        if(((x+dx[index]) < 0) || ((x+dx[index]) > (six - 1)) ||
           ((y+dy[index]) < 0) || ((y+dy[index]) > (siy - 1)) ||
           ((z+dz[index]) < 0) || ((z+dz[index]) > (siz - 1))) {
          //printf("unclassified point [%d,%d,%d]\n",x,y,z);
          //return(-111);
        }
        if(!input[index]) {
          input[index] = input[(z+dz[index])*dd+(y+dy[index])*six+
                               (x+dx[index])];
        }
      }
    }
  }
  *odx = dx;
  *ody = dy;
  *odz = dz;
  free(sq);
  return(0);
}
