#include "dsltool.h"
#include "ui.h"

#include <deque>








/**********************************************************************
* UNDO Related
**********************************************************************
int undoSize = 50;
int undoOldest=-1,undoNewest=-1,undoCurrent = -1;
vector <DSLObject2D *> undoRing;

bool canUndo() {
	return (undoOldest != undoCurrent);
}

bool canRedo() {
	return (undoNewest != undoCurrent);
}

void updateUndoMenus() {
   if(canUndo()) {
		miUndo->activate();
      bUndo->activate();
   }
   else {
		miUndo->deactivate();
      bUndo->deactivate();
   }
	
   if(canRedo()) {
		miRedo->activate();
      bRedo->activate();
   }
   else {
      miRedo->deactivate();
      bRedo->deactivate();
   }
}

void pushModel() {
	// If this is being done for the first time
	if(undoCurrent < 0) {
		undoOldest = 0;
		undoNewest = 0;
		undoCurrent = 0;
		
		for(int i=0;i<undoSize;i++)
			undoRing.push_back(NULL);
	}
	else {
		undoCurrent = (undoCurrent+1) % undoSize;
		undoNewest = undoCurrent;
		if(undoOldest == undoCurrent) 
			undoOldest = (undoOldest+1) % undoSize;
		
      undoRing[undoCurrent] = new DSLObject2D(model);
      removeAllDataFromModel(*undoRing[undoCurrent]);

		miUndo->activate();
      bUndo->activate();
	}
	
	updateUndoMenus();
}

void undo() {
	if(canUndo()) {
		// Make copy of current model
      DSLObject2D *copy = new DSLObject2D(model);
      removeAllDataFromModel(*copy);
      
      // Set model to stored model
      model = *undoRing[undoCurrent];
      addAllDataToModel(model);

      // Store model being undone onto undo list
      undoRing[undoCurrent] = copy;

		// 
      undoCurrent = (undoCurrent-1) % undoSize;
	}
	updateUndoMenus();	
}


void redo() {
	if(canRedo()) {
		undoCurrent = (undoCurrent+1) % undoSize;
      model = undoRing[undoCurrent];
      addAllDataToModel(model);

	}
	updateUndoMenus();	
}

void purgeUndo() {
   undoRing.clear();
}
*/