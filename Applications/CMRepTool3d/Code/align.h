#ifndef _ALIGN_H_
#define _ALIGN_H_

#ifdef WIN32

#define ALIGN_PRE __declspec(align(16)) 
#define ALIGN_POST  

#else

#define _aligned_free(a) free(a);
#define ALIGN_PRE   
#define ALIGN_POST  __attribute__((aligned(16)))

inline void *_aligned_malloc(size_t n, size_t align) 
{
  void *out;
  if(posix_memalign(&out,align,n) == 0) 
    return out;
  else
    return NULL;
}

#endif // Platform

#endif // _ALIGN_H_
