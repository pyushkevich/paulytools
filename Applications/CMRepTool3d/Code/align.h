#ifndef _ALIGN_H_
#define _ALIGN_H_

#ifdef WIN32
  #define ALIGN32_PRE   ALIGN32
  #define ALIGN32_POST  
#else
  #define _aligned_free(a) free(a);
  #define ALIGN32_PRE   
  #define ALIGN32_POST  __attribute__((aligned(32)))

  inline void *_aligned_malloc(size_t n, size_t align) {
    void *out;
    if(posix_memalign(&out,align,n) == 0) 
      return out;
    else
      return NULL;
  }

#endif

#endif
