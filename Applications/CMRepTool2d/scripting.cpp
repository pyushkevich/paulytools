#include "scripting.h"
#include <stdio.h>
#include <stdarg.h>

void UserInterfaceScript::command(char *cmd,bool fromScript) {
   // If simple command interface then just call and be done
   if(commandProcessor) {
      commandProcessor(cmd,fromScript);
      return;
   }

   char *src = cmd;
   //static char **chunks = new char *[40];
   static char *chunks[40];
   int chunk = -1;   
   
   bool inquotes = false;
   bool lastblank = true;
   

   while(*src != 0) {
      if(inquotes) {
         if(*src=='\"') {
            inquotes = false;
            lastblank = true;
            *src = 0;
         }
      }
      else if(*src == ' ' || *src == '\t' ||  *src == ',' || *src == '\n'  || *src == '\r' ) {
         if(!lastblank) {
            *src = 0;
            lastblank = true;
         }
      }
      else if(*src == '\"') {
         if(!lastblank) {
            *src = 0;
            lastblank = true;
         }
         inquotes = true;
         chunks[++chunk] = src+1;
      }
      else {
         if(lastblank) {
            chunks[++chunk] = src;
            lastblank=false;
         }
      }
      src++;
   }

   commandProcessorChunks(chunk+1,chunks,fromScript);
}

void UserInterfaceScript::startRecording() {
   // If already recording, return
   if(record)
      return;

   // If paused then just resume recording.  If not paused then clear the script first
   if(!paused) 
      macroScript.clear();

   record = true;
   paused = false;

   sliderTag = this;
}

void UserInterfaceScript::pauseRecording() {
   // If not recording, or already paused, return
   if(!record || paused)
      return;

   // Pause
   paused = true;
   record = false;
}

void UserInterfaceScript::stopRecording() {
   // To stop we must be either recording or paused
   if(record || paused) {
      paused = false;
      record = false;
   }
}

string commLinSub(string command,string *args) {
	int idx = 0;
	string result;
	
	try {
		while (1) {
			int dix = (int) command.find("$",idx);
			if(dix < 0) {
				result.append(command.substr(idx,-1));
				return result;
			}
			else {
				result.append(command.substr(idx,dix-idx));
				char c = command.at(dix+1);
				if(c=='$')
					result.append("$");
				else
					result.append(args[c-'1']);
				idx = dix + 2;
			}
		}
	}
	catch(...) {
		return command;
	}
}

void UserInterfaceScript::run(int argc,char *argv[]) {
	int i;

	// Need a copy of the args
	string args[10];
	for(i=0;i<argc;i++)
		args[i] = argv[i];

   for(i=0;i<macroScript.size();i++) {      
		// Do command line parameter substitution
		string sub = commLinSub(macroScript[i],args);
		completeScript.push(sub.c_str());
		
		char *tmp = new char[sub.length()+256];
		strcpy(tmp,sub.c_str());
      command(tmp,true);      
		delete tmp;
   }
}

bool UserInterfaceScript::save(const char *fileName) {
   FILE *f = fopen(fileName,"wt");

   if(f==NULL)
      return false;

   for(int i=0;i<macroScript.size();i++) 
      fprintf(f,"%s\n",macroScript[i]);

   fclose(f);
   return true;
}

bool UserInterfaceScript::saveHistory(const char *fileName) {
   FILE *f = fopen(fileName,"wt");

   if(f==NULL)
      return false;

   for(int i=0;i<completeScript.size();i++) 
      fprintf(f,"%s\n",completeScript[i]);

   fclose(f);
   return true;
}

bool UserInterfaceScript::load(const char *fileName) {
   char command[256];
   FILE *f = fopen(fileName,"rt");

   if(f==NULL)
      return false;

   record = paused = false;
   macroScript.clear();

   while(!feof(f)) {
      fgets(command,255,f);
      macroScript.push(command);      
   }
   
   fclose(f);
   return true;
}

bool UserInterfaceScript::isRecording() {
   return record;
}

bool UserInterfaceScript::isPaused() {
   return paused;
}

UserInterfaceScript::UserInterfaceScript(void (*commandProcessor)(const char *,bool fromScript)) {
   this->commandProcessor = commandProcessor;
   this->commandProcessorChunks = NULL;

   record = false;
   paused = false;
   sliderTag = this;
}

UserInterfaceScript::UserInterfaceScript(void (*commandProcessorChunks)(int argc,char **,bool fromScript)) {
   this->commandProcessorChunks = commandProcessorChunks;
   this->commandProcessor = NULL;

   record = false;
   paused = false;
   sliderTag = this;
}

void UserInterfaceScript::uiCommand(const char *fmt,...) {
   va_list marker;   
   char cmd[256];
   
   // Make a command out of passed in ...
   va_start(marker, fmt);     
   vsprintf(cmd,fmt,marker);
   va_end(marker);

   // Record the command in history and macro
   if(record)
      macroScript.push(cmd);
   completeScript.push(cmd);

   // Process the command
   command(cmd,false);

   sliderTag = this;
}

void UserInterfaceScript::uiSliderCommand(void *tag,const char *fmt,...) {
   va_list marker;   
   char cmd[256];

   // Make a command out of passed in ...
   va_start(marker, fmt);     
   vsprintf(cmd,fmt,marker);
   va_end(marker);

   
   // See if the command is repetative and if it is, pull the last one from the script
   if(tag==sliderTag && tag!=this) {
      if(record)
         macroScript.pop();
      completeScript.pop();
   }

   // Record the command in history and macro
   if(record)
      macroScript.push(cmd);
   completeScript.push(cmd);

   // Process the command
   command(cmd,false);

   sliderTag = tag;
}

void UserInterfaceScript::runCommand(const char *incmd) {
   char cmd[256];
   strcpy(cmd,incmd);
   command(cmd,true);
}

   
