#ifndef _DSL_UNDO_H_
#define _DSL_UNDO_H_

#include <deque>

template <class T>
class UndoBuffer {
private:
   deque <T> undoStack;
   deque <T> redoStack;
   int undoSize,redoSize;
   
public:
   UndoBuffer(int undoSize=50)
   {
      this->undoSize = undoSize;
   }
   
   void push(const T &t) {
      int uSize = undoStack.size();

      // Keep the stack constant size
      if(uSize >= undoSize) 
         undoStack.pop_front();
      
      // Clear the redo stack
      redoStack.clear();
      
      // Push the current model onto the undo stack
      undoStack.push_back(t);
   }
   
   T undo(const T &t){
      if(canUndo()) {
         // Current model goes onto the redo stack
         redoStack.push_back(t);
         
         // Get the top of the undo stack
         T rtn = undoStack.back();
         
         // Pop the guy off the undo stack
         undoStack.pop_back();
         
         // Return the guy
         return rtn;
      }
      else 
         return t;
   }

   T redo(const T &t){
      if(canRedo()) {
         // If redo has elements on it that means undo stack can't be full
         assert(undoStack.size() <= undoSize);
            
         // Current guy goes on the undo stack
         undoStack.push_back(t);
         
         // Get the top of redo stack
         T rtn = redoStack.back();
         
         // Pop the guy off the undo stack
         redoStack.pop_back();
         
         // Return the guy
         return rtn;
      }
      else 
         return t;
   }   

	// Returns the pointer to the next object to be undone or NULL
	T* getTail() {
		if(canUndo())
			return &undoStack.back();
		return NULL;
	}
   
   bool canUndo(){
      return (undoStack.size() > 0);
   }

   bool canRedo(){
      return (redoStack.size() > 0);
   }
};

/*
void pushModel();
void undo();
void redo();
void purgeUndo();
*/

#endif

