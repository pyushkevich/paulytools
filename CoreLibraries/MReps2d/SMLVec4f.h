#ifndef _SML_VEC_4f_
#define _SML_VEC_4f_

class MySMLVec4f {
    public:
	float x,y,z,w;
	float *data() { return (float *)this; }
	MySMLVec4f() {x=y=z=w=0;}
	MySMLVec4f(float x0,float y0,float z0,float w0) {x=x0;y=y0;z=z0;w=w0;}

	float Dot(const MySMLVec4f &v) {
	    return x*v.x + y*v.y + z*v.z + w*v.w;
	}

	void Set(float f1,float f2,float f3,float f4) {
	    x = f1;y = f2;z = f3;w = f4;
	}
};

class MySMLVec3f {
    public:
	float x,y,z;
	float *data() { return (float *)this; }
	MySMLVec3f() {x=y=z=0;}
	MySMLVec3f(float x0,float y0,float z0) {x=x0;y=y0;z=z0;}

	float Dot(const MySMLVec3f &v) {
	    return x*v.x + y*v.y + z*v.z;
	}

	void Set(float f1,float f2,float f3) {
	    x = f1;y = f2;z = f3;
	}

  MySMLVec3f operator+ (const MySMLVec3f &v) {
    return MySMLVec3f(x+v.x,y+v.y,z+v.z);
  }

  MySMLVec3f operator- (const MySMLVec3f &v) {
    return MySMLVec3f(x-v.x,y-v.y,z-v.z);
  }

  MySMLVec3f operator* (float k) {
    return MySMLVec3f(x*k,y*k,z*k);
  }
};

#endif
