#ifndef ARRAY_2D_
#define ARRAY_2D_

#include <cassert>

#include <xmmintrin.h>

#ifdef _DEBUG
#define DEBUG_ASSERT(expression) assert(expression)
#else
#define DEBUG_ASSERT(expression) ;
#endif



/*****************************************************************
A Fast 2D Array Implementation
****************************************************************/
template <class T> class Array2D {
private:
	T *array;
	T **rows;

	int w,h;

	void alloc() {
		if(array!=NULL) 
			dealloc();
		array = new T [h*w];
		rows = new T*[h];

		rows[0] = array;
		for(int y=1;y<h;y++) {
			rows[y]=rows[y-1] + w;
		}
	}

	void dealloc() {
		DEBUG_ASSERT(array);
		delete [] array;
		delete [] rows;
		array = NULL;
		rows = NULL;
	}

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

	Array2D(const Array2D &copy) {
		array = NULL;
		rows = NULL;
		w = h = 0;
		*this = copy;
	};

	~Array2D() {
		if(array)
			dealloc();
	}

	T& operator ()(int x,int y) { 
		DEBUG_ASSERT(array && x>=0 && y>=0 && x<w && y<h);
		return rows[y][x];
	}

	T operator () (int x,int y) const { 
		DEBUG_ASSERT(array && x>=0 && y>=0 && x<w && y<h);
		return rows[y][x];
	}

	// This method returns a row in the array as a pointer.  It's great for fast access
	T* row(int y) {
		return rows[y];
	}

	int width() const {
		return w;
	}

	int height() const {
		return h;
	}

	int size() const {
		return w*h;
	}

	void resize(int w,int h) {
		if(this->w != w || this->h != h) {
			this->w = w;
			this->h = h;
			alloc();
		}
	}

	void setAll(const T &value) {
		for(int j=0;j<h;j++)
			for(int i=0;i<w;i++)
				rows[j][i] = value;
	}

	T* getData() {
		return array;
	}

	Array2D& operator =(const Array2D &copy) {
		resize(copy.w,copy.h);
		for(int j=0;j<h;j++)
			for(int i=0;i<w;i++)
				rows[j][i] = copy.rows[j][i];
		return *this;
	}
};

/**
 * This is like a 2D array with an additional 1D array at the end (the auxillary array).
 * I use this for OpenGL vertex array functionality when there are extra vertices outside 
 * of the main grid
 */
template <class T> class Aux2D {  
private:
	T *array;
	T **rows;
	T *auxPtr;
	// int *rowIdx;

	int w,h;
	int auxLen;
	int len;

	void alloc() {
		if(array!=NULL) 
			dealloc();
		len = h*w+auxLen;
		
		// Allocate the array on an aligned boundary!
		array = (T*) _aligned_malloc(sizeof(T)*len,32);
		
		rows = new T*[h];

		rows[0] = array;
		for(int y=1;y<h;y++) {
			// rowIdx[y] = rowIdx
			rows[y]=rows[y-1] + w;
		}
		auxPtr = rows[h-1]+w;
	}

	void dealloc() {
		DEBUG_ASSERT(array);
		_aligned_free(array);
		delete [] rows;
		array = NULL;
		rows = NULL;
		auxPtr = NULL;
	}

public:
	Aux2D() {
		array = NULL;
		rows = NULL;
		auxPtr = NULL;
		len = w = h = auxLen = 0;
	}

	Aux2D(int w,int h,int auxLen) {
		this->w = w;
		this->h = h;
		this->auxLen = auxLen;
		array = NULL;
		rows = NULL;
		alloc();
	}

	Aux2D(const Aux2D &copy) {
		array = NULL;
		rows = NULL;
		auxPtr = NULL;
		len = w = h = auxLen = 0;
		*this = copy;
	};

	~Aux2D() {
		if(array)
			dealloc();
	}

	T& operator ()(int x,int y) { 
		DEBUG_ASSERT(array && x>=0 && y>=0 && x<w && y<h);
		return rows[y][x];
	}

	T& operator ()(int idx) { 
		DEBUG_ASSERT(array && idx < len);
		return array[idx];
	}

	const T &operator () (int x,int y) const { 
		DEBUG_ASSERT(array && x>=0 && y>=0 && x<w && y<h);
		return rows[y][x];
	}

	const T &operator ()(int idx) const { 
		DEBUG_ASSERT(array && idx < len);
		return array[idx];
	}

	// This method returns a row in the array as a pointer.  It's great for fast access
	T* row(int y) {
		return rows[y];
	}

	// This method returns the pointer to aux data
	T* aux() const {
		return auxPtr;
	}

	// This method returns the pointer to aux data
	T& aux(int idx) {
		return auxPtr[idx];
	}

	int width() const {
		return w;
	}

	int height() const {
		return h;
	}

	int auxSize() {
		return auxLen;
	}

	int size() {
		return len;
	}

	// Get the offset of an element
	int offset(int x,int y) {
		return (rows[y]-array) + x;
	}

	// Get the offset of an aux element
	int auxOffset(int z) {
		return (auxPtr - array) + z;
	}

	void resize(int w,int h,int auxLen) {
		if(this->w != w || this->h != h || this->auxLen != auxLen) {
			this->w = w;
			this->h = h;
			this->auxLen = auxLen;
			alloc();
		}
	}

	void setAll(const T &value) {
		for(int j=0;j<len;j++)
			array[j] = value;
	}

	T* getData() {
		return array;
	}

	Aux2D& operator =(const Aux2D &copy) {
		resize(copy.w,copy.h,copy.auxLen);
		for(int j=0;j<len;j++)
			array[j] = copy.array[j];
		return *this;
	}
};



#endif