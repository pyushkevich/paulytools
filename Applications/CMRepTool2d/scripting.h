#ifndef _UI_SCRIPTING_
#define _UI_SCRIPTING_

#include <vector>

using namespace std;

class StringStack {
   vector <char *> stk;
public:
   void push(const char *s) {
      char *cpy = new char[strlen(s)+1];
      strcpy(cpy,s);
      stk.push_back(cpy);
   }

   void pop() {
      if(stk.size()>0) {
         delete stk[stk.size()-1];
         stk.pop_back();
      }
   }

   void clear() {
      while(size())
         pop();
   }

   char *operator[](int i) {
      return stk[i];
   }

   int size() {
      return (int) stk.size();
   }
};

class UserInterfaceScript {

private:
   bool record;
   bool paused;
   // vector <String> macroScript;
   StringStack macroScript;
   StringStack completeScript;

   void *sliderTag;
   void (*commandProcessor)(const char *,bool);
   void (*commandProcessorChunks)(int argc,char **,bool);

   // Run a command (used internally, or can be done from a command line
   void command(char *cmd,bool fromScript);

public:
   UserInterfaceScript(void (*commandProcessor)(const char *,bool fromScript));
   UserInterfaceScript(void (*commandProcessorChunks)(int argc,char **,bool fromScript));

   // Record macro script
   void startRecording();
   void pauseRecording();
   void stopRecording();

   // Run script
   void run(int argc = 0,char *argv[] = NULL);
	
   // Load/Save script
   bool load(const char *fileName);
   bool save(const char *fileName);
   bool saveHistory(const char *fileName);

   // Called by the user interface to execute command and add it to a script
   void uiCommand(const char *fmt,...);

   // Called by a widget like a slider that changes its value multiple times, so that only the last call is 
   // stored in the script.  It uses the tag variable to keep track of the last caller.  So if the tag matches between
   // two calls to this routine (and no other command is executed in between) only the second call is used.
   void uiSliderCommand(void *tag,const char *fmt,...);

   void runCommand(const char *cmd);

   bool isRecording();
   bool isPaused();
};

#endif // _UI_SCRIPTING_

