#ifndef ARRAY_2D_
#define ARRAY_2D_

#include <cassert>

/*****************************************************************
 A Fast 2D Array Implementation
  ****************************************************************/
template <class T> class Array2D {
private:
   T *array;
   T **rows;

   size_t w,h;
   
   void alloc() {
      if(array!=NULL) 
         dealloc();
      array = new T [h*w];
      rows = new T*[h];
      
      rows[0] = array;
      for(size_t y=1;y<h;y++) {
         rows[y]=rows[y-1] + w;
      }
   }

   void dealloc() {
      assert(array);
      delete [] array;
      delete [] rows;
      array = NULL;
      rows = NULL;
   }
   
   Array2D(const Array2D &copy) {
      *this = copy;
   };

public:
   Array2D() {
      array = NULL;
      rows = NULL;
      w = h = 0;
   }

   Array2D(int w,int h) {
      this->w = w;
      this->h = h;
      array = NULL;
      rows = NULL;
      alloc();
   }

   Array2D(int w,int h,T init_value) {
      this->w = w;
      this->h = h;
      array = NULL;
      rows = NULL;
      alloc();
      
      for(int x=0;x<w;x++) {
         for(int y=0;y<h;y++) {
            (*this)(x,y) = init_value;         
         }
      }
   }

   ~Array2D() {
      if(array)
         dealloc();
   }
   
   T& operator ()(size_t x,size_t y) { 
      assert(array && x<w && y<h);
      return rows[y][x];
   }
   
   T operator () (size_t x,size_t y) const { 
      assert(array && x<w && y<h);
      return rows[y][x];
   }

   // This method returns a row in the array as a pointer.  It's great for fast access
   T* row(size_t y) {
      return rows[y];
   }
   
   int width() const {
      return (int)w;
   }
   
   int height() const {
      return (int)h;
   }

   void resize(int w,int h) {
      this->w = w;
      this->h = h;
      alloc();
   }

   T* getData() {
      return array;
   }

private:

   Array2D& operator =(const Array2D &copy) {
      resize(copy.width(),copy.height());
      for(int j=0;j<h;j++)
         for(int i=0;i<w;i++)
            rows[j][i] = copy(i,j);
	  return *this;
   }
   
	
};


#endif

