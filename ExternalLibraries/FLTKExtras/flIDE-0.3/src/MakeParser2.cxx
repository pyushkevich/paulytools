#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
#include <fcntl.h>
#include <ctype.h>
#include <string.h>
#include <regex.h>
#include <sys/types.h> 
#include <sys/wait.h>
#include <pthread.h>
#include <FL/Fl.H> 
#include <FL/fl_ask.H>
#include "MakeParser.H"

#include "xpopen.h" /* an extended popen implementation to easily kill the created child(s) */

/* ID's of the different categories of the recognized/displayed messages */
// plain text (no special meaning/display type)
#define FLMAKE_MSG_PLAIN	0
// internal flmake status messages
#define FLMAKE_MSG_ISTATUS	1
#define FLMAKE_MSG_IERROR	2
// make recognized message categories
#define FLMAKE_MSG_WARNING	3
#define FLMAKE_MSG_ERROR	4
#define FLMAKE_MSG_INFO		5
#define FLMAKE_MSG_RUNTARGET0	6
#define FLMAKE_MSG_RUNTARGET1	7


pthread_t makeThread=0;
pthread_mutex_t fastmutex = PTHREAD_MUTEX_INITIALIZER;

bool make_busy;
bool make_request;

char cmd[1024];


/** insert string newtxt into the string srcdest replacing the substring that
   starts at nschars of size ntlen chars, returns srcdest, srcdest must be long 
   enough to  contain the whole output string */
char *strinstr(char *srcdest,char *newtxt, size_t posinsrc,size_t nschars,size_t ntlen);

void MakeParser::PrevErrorMsg(void)
{
	mErrorTab->show();
	PrevMessage(mMsgHeaders[FLMAKE_MSG_ERROR]);
}

void MakeParser::NextErrorMsg(void)
{
	mErrorTab->show();
	NextMessage(mMsgHeaders[FLMAKE_MSG_ERROR]);
}

void MakeParser::PrevWarningMsg(void)
{
	mErrorTab->show();
	PrevMessage(mMsgHeaders[FLMAKE_MSG_WARNING]);
}

void MakeParser::NextWarningMsg(void)
{
	mErrorTab->show();
	NextMessage(mMsgHeaders[FLMAKE_MSG_WARNING]);
}

void MakeParser::NextMessage(const char* identTok)
{
	int sline=mErrorBrowser->value();
	int nlines=mErrorBrowser->size();
	const char* text;
	do {
	 sline++;
	 text=mErrorBrowser->text(sline);
	 if(text!=NULL && strstr(text,identTok)){
	     mErrorBrowser->value(sline);
	     if(mEditOnMsgBrowse) {
	        SelectErrorLine();
	     }
	     text=NULL;
	  }
	} while( sline<=nlines && text!=NULL);
} 

void MakeParser::PrevMessage(const char* identTok)
{
	int sline=mErrorBrowser->value();
	const char* text;

	do {
	  sline--;
	  text=mErrorBrowser->text(sline);
	  if(text!=NULL && strstr(text,identTok)){
	     mErrorBrowser->value(sline);
	     if(mEditOnMsgBrowse) {
	        SelectErrorLine();
	     }
	     text=NULL;
	  }
	} while( sline>0 && text!=NULL);
}

char* MakeParser::LookForArg(const char* arg,const char* makeoutline)
{
	char* found1=strstr(makeoutline,arg);
	char* ret=0;
	
	if (found1) {
		char* ptr=found1+strlen(arg);
		char* found2=strstr(ptr," ");
		int len=0;	
		if (found2) {
			len=found2-ptr;
		}else if((found2=strstr(ptr,"\t"))){
			len=found2-ptr;
		}else if((found2=strstr(ptr,"\n"))){
			len=found2-ptr;
		}else {
			len=strlen(ptr);
		}
		len+=strlen(mCurrDir);
		ret=new char[len+2];
		strcpy(ret,mCurrDir);
		strcat(ret,"/");
		strncpy(ret+strlen(mCurrDir),ptr,len-strlen(mCurrDir));
		ret[len]='\0';
	}
	return ret;
}

int MakeParser::FindTokInBrowser(Fl_Browser *mBrowser,const char *identTok) 
{
	int sline=0;
	int nlines=mBrowser->size();
	const char* text;
	if(nlines==0) 
	  return false;
	  
	do {
	 sline++;
	 text=mBrowser->text(sline);
	 if(text!=NULL && strstr(text,identTok)){
	     text=NULL;
	  }
	} while( sline<=nlines && text!=NULL);
        return (text==NULL && sline<=nlines);
}

void MakeParser::LookForExecutable(const char* makeoutline)
{
	char buf[1024]; 
	char *msgStart;
	char* found=LookForArg(" -o ",makeoutline);
	if (found) {
	    if(strstr(found,".o")==NULL && 
	       strstr(found,".so")==NULL) {
		pthread_mutex_lock(&fastmutex);		   
		    strcpy(buf,((mRunTargetsBrowser->size())&(1))?
		                 mMsgHeaders[FLMAKE_MSG_RUNTARGET0]:
				 mMsgHeaders[FLMAKE_MSG_RUNTARGET1]);
		    msgStart=buf+strlen(buf)-1;		 
		    strcat(buf,found);
		    strcat(buf,"\t");
                    if(mCleanRunTargets ||
		       !FindTokInBrowser(mRunTargetsBrowser,msgStart)) {	
		       mRunTargetsBrowser->add(buf);
		    }
		pthread_mutex_unlock(&fastmutex);
	    }
	}
}

void* make_thread_routine(void* ptr)
{
	((MakeParser*)(ptr))->MakeThread();
}

/* Kill the popened process, before killing the thread.
   This is required since killing a thread doesn't kills the forked child
   processes belonging to him. */
void kill_thread_routine(void *arg)
{
	XPopenDescr *makeout=(XPopenDescr *)arg;
	xpkillclose_XPD(makeout);

	make_request = 0;
	make_busy=0;
}

/* pthread safe fgets ... 
  This is necesary to be sure that the cancel request is procesed in the systems
  which don't have a thread safe gets func (i.e. linux+libc5) */
char *tsfgets(char *s, int size, FILE *stream)
{
	char *resp;
	pthread_testcancel();
	resp=fgets(s,size,stream);
	pthread_testcancel();
	return resp;
}

void MakeParser::MakeThread() 
{
	char buf[1024];
	XPopenDescr *makeout;

	pthread_mutex_lock(&fastmutex);
	  mErrorBrowser->clear();
	  if(mCleanRunTargets)
	      mRunTargetsBrowser->clear();
	pthread_mutex_unlock(&fastmutex);

	/* We use here our own extended popen,pclose functions (xpkillclose here) since it's
	  the easiest way of knowing the process id of the child to be killed,
	  (the standard system don't give back the create child's id making overly
	  complicated to kill the child process) */
	makeout = xpopen_XPD(cmd,"r");

	if(makeout==NULL) {
	   pthread_mutex_lock(&fastmutex);
	     strcpy(buf,mMsgHeaders[FLMAKE_MSG_IERROR]);
	     strcat(buf,"ERROR: Couldn't execute make command");
	     mErrorBrowser->add(buf);
	   pthread_mutex_unlock(&fastmutex);
	   make_busy=0;
	}

	 /* set cleanup routine to be executed if the thread gets killed before
	 finishing the make. (i.e. using the stop make button) */
	 pthread_cleanup_push(kill_thread_routine, (void *)makeout);
	
	 while (tsfgets(buf,1024,(*makeout).fd)) {
	       char filename[1024];
	       int linenumber;
	       int msgtype;
	       const char* msg;
	       char cpy[1024];
	       char expline[1024];
	       filename[0]='\0';

	       char* addstr=0;
	
	       msg=SplitLineRegExp(buf,filename,expline,&linenumber,&msgtype,1);

	       if (linenumber) {
	           switch(msgtype) {
	              case FLMAKE_MSG_WARNING:
	        	    strcpy(cpy,mMsgHeaders[FLMAKE_MSG_WARNING]);
	        	    strcat(cpy,expline);
	        	    addstr = cpy;
	        	    break;
	              case FLMAKE_MSG_ERROR:
	        	    strcpy(cpy,mMsgHeaders[FLMAKE_MSG_ERROR]);
	        	    strcat(cpy,expline);
	        	    addstr = cpy;
	        	    break;
	              case FLMAKE_MSG_INFO:
	        	    strcpy(cpy,mMsgHeaders[FLMAKE_MSG_INFO]);
	        	    strcat(cpy,expline);
	        	    addstr = cpy;
	        	    break;
	              default:
	        	    addstr = expline;
	        	    break;
	           }
	       } else {
	           addstr = expline;
	       }
	       if (addstr) {
	           pthread_mutex_lock(&fastmutex);
	             mErrorBrowser->add(addstr);
	           pthread_mutex_unlock(&fastmutex);
	       }
	       LookForExecutable(buf);
	       pthread_testcancel();
	}
	pthread_mutex_lock(&fastmutex);
	  strcpy(buf,mMsgHeaders[FLMAKE_MSG_ISTATUS]);
	  strcat(buf,"DONE");
	  mErrorBrowser->add(buf);
	pthread_mutex_unlock(&fastmutex);
	pthread_cleanup_pop(0);
	make_busy=0;
	mCleanRunTargets=false;
	xpclose_XPD(makeout);
}

void MakeParser::RunMake(char* cmdin) 
{
	strncpy(cmd,cmdin,1024-6);
	strncat(cmd," 2>&1",1024-strlen(cmd));	
	make_request = 1;
	if (make_busy) {
	       // recursive call: cancel the thread and issue a request.
	       pthread_cancel(makeThread);
	       pthread_join(makeThread,0);
	       make_request = 1;
	       make_busy=0;
	       return;
	}

	while (make_request) {
	       make_request=0;
	       make_busy = 1;

	       pthread_create(&makeThread,NULL,make_thread_routine,(void *)this);

	       while (make_busy) {
	               pthread_mutex_lock(&fastmutex);
	               Fl::check();
	               pthread_mutex_unlock(&fastmutex);
	       }
	}
}

void MakeParser::Make(void)
{
	char makecmd[1024];
	mErrorTab->show();
	strcpy(makecmd,"make ");
	strcat(makecmd,mMakeArgs);
	RunMake(makecmd);
}

void MakeParser::StopMake(void) 
{
	char msgtmp[1024];
	if (make_busy) {
	   /* cancel the thread and kill the child make  proces if there's one running
	     this is actually done in the kill_thread_routine func that process the thread cancel
	     request */
	   mErrorTab->show();
	   pthread_cancel(makeThread);
	   pthread_join(makeThread,NULL);
	   make_request = 0;
	   make_busy=0;
	   strcpy(msgtmp,mMsgHeaders[FLMAKE_MSG_ISTATUS]);
	   strcat(msgtmp,"STOPPED");
	   mErrorBrowser->add(msgtmp);
	   mCleanRunTargets=false;
	}
}

void MakeParser::RunWithConsole(void)
{
//TODO: add a shell, xterm variable to add to the command line..
}

void MakeParser::RunInDebugger(void)
{
//TODO: add a debug load/start to add to the command line..
}

int MakeParser::SplitRunTargetsLine(const char *text,char *exename,char *params)
{
  	int i,j;
	if(text!=0) {
	   i=0;
	   if (text[i]=='@') {
	       i++;
	       while (text[i] && !(text[i-1]=='@' && text[i]=='.')) {
	     	      i++;
	       }
	       if (text[i]=='.')
	     	     i++;
	   }
	   j=i;
	   while (text[i] && !(text[i]=='\t')) {
	     	      i++;
	   }
	   strncpy(exename,&text[j],i-j);
	   exename[i-j]='\0';
	   if (text[i]=='\t')
	     	     i++;	
	   strcpy(params,&text[i]);
	   return true;
	} else {
 	   return false;
	}

}

void MakeParser::ChangeRunTargetsParams(void)
{
	const char* text;
	int line=mRunTargetsBrowser->value();
	
	mRunTab->show();

        if(line==0) 
	  line=mRunTargetsBrowser->size();

        mRunTargetsBrowser->select(line);
	 
	if (mRunTargetsBrowser->size()>0 && line<=mRunTargetsBrowser->size()) {
	   char buf[1024];
	   char buf1[1024];
	   const char *buft;
           text=mRunTargetsBrowser->text(line);
           if(SplitRunTargetsLine(text,buf,buf1)) {
	      buft=fl_input("Command line parameters for \"%s\":",buf1,buf);
	      if(buft!=NULL) {
	          char buf3[1024];
		  strcpy(buf3,(line&(1))?
		                 mMsgHeaders[FLMAKE_MSG_RUNTARGET1]:
				 mMsgHeaders[FLMAKE_MSG_RUNTARGET0]);
	          strcat(buf3,buf);
	          strcat(buf3,"\t");
	          strcat(buf3,buft);
    		  mRunTargetsBrowser->text(line,buf3);
	     } else if (buft==NULL) {
	              return;
	     }
           }
	}
}

void MakeParser::Run(void)
{
	const char *new_argv[4];
	const char* text;
	int line=mRunTargetsBrowser->value();
	
	mRunTab->show();

        if(line==0) 
	  line=mRunTargetsBrowser->size();
	 
	if (mRunTargetsBrowser->size()>0 && line<=mRunTargetsBrowser->size()) {
	   char buf[1024];
	   char buf1[1024];
	   buf1[0]='\0';
           text=mRunTargetsBrowser->text(line);
           if(SplitRunTargetsLine(text,buf,buf1)) {
              char buf3[1024];
	      if(buf[0]!='/') {
	        strcpy(buf3,buf);
	        new_argv[0]=buf;
                new_argv[1]=(buf1[0]!='\0')?buf1:NULL;
	        new_argv[2]=NULL;
	      } else {
	        int i;
	        for(i=strlen(buf)-1; i>=0 && buf[i]!='/'; i--);
		buf[i]='\0';
	        sprintf(buf3,"cd %s ; %s %s",buf,&buf[i+1],buf1);
		strcpy(buf,buf3);
		strcpy(buf3,"/bin/sh");
		new_argv[0]="sh";
	        new_argv[1]="-c";
	        new_argv[2]=buf;
	        new_argv[3]=NULL;
	      }	
 	      if (fork()==0) {
	         (void) execv(buf3,(char *const *)new_argv);
	      }
 	      return;
	   } 
	}
	fl_alert("Flmake Wasn't capable of finding a runnable file:\n\n"
	           " -Either the make process didn't create one or,\n"
		   " -Flmake wasn't able to identify one.");
}

void MakeParser::MakeClean(void)
{
	mErrorTab->show();
	mCleanRunTargets=true;
	RunMake("make clean");
}

int MakeParser::CompileRegEx(regex_t *preg,const char* regex)
{
	int err,errsize;
	char ebuff[255];
	err=regcomp(preg,regex,REG_EXTENDED);
	if(err==0)
	  return 1;
	errsize=regerror(err,preg,ebuff,255);
	fprintf(stderr,"RegExp Error: %s\n",ebuff);
	return 0;
}

void MakeParser::InitRegex(void)
{
 mNPreg=0;

/* the first two are the make enter and leaving directoryes*/
 CompileRegEx(&mPreg[mNPreg++],"^make[[]([0-9]+)[]]: Entering directory.+`(.+)'");
 CompileRegEx(&mPreg[mNPreg++],"^make[[]([0-9]+)[]]: Leaving directory.+`(.+)'");

// informative messages
 CompileRegEx(&mPreg[mNPreg++],"^In file included from (.+):([0-9]+)[,:]");
 CompileRegEx(&mPreg[mNPreg++],"^[ \t]+from (.+):([0-9]+)[,:]");

 mINPreg=mNPreg;

/* cc, CC, flex, bison warnings*/
 CompileRegEx(&mPreg[mNPreg++],"^.*\"(.+)\".*line ([0-9]+).*warning.*");  
/* gcc,g++,g77, gmake, javac warnings*/
 CompileRegEx(&mPreg[mNPreg++],"^(.+):([0-9]+):.+warning.*");
 CompileRegEx(&mPreg[mNPreg++],"^(.+):([0-9]+):.+Note.*");

 mWNPreg=mNPreg;
 
/* perl */
 CompileRegEx(&mPreg[mNPreg++],"^.+at (.+) line ([0-9]+).*");  

/* gcc,g++,g77, gmake, javac errors*/
 CompileRegEx(&mPreg[mNPreg++],"^.+[ \t]+(.+):([0-9]+)[:,].*"); 
 CompileRegEx(&mPreg[mNPreg++],"^(.+):([0-9]+)[:,].*"); 

/* cc, CC, flex, bison*/
 CompileRegEx(&mPreg[mNPreg++],"^.*\"(.+)\".*line.*([0-9]+).*");  
 
/*make ignored error */
 CompileRegEx(&mPreg[mNPreg++],
 	"^(.+):[[]?([0-9]+)[]]?:[*][*][*].*\\(ignored\\).*"); 
 
/*make error */
 CompileRegEx(&mPreg[mNPreg++],"^(.+):[[]?([0-9]+)[]]?: [*][*][*] .*"); 
}

int MakeParser::TestRegSet(const char *buff,int nmatch, regmatch_t pmatch[])
{
	int i, resp;

	for(i=0, resp=-1; i<mNPreg && resp==-1; i++ ){
	  resp=(regexec(&mPreg[i],buff,nmatch,pmatch,0)==0)?i:-1;
	}
	return resp;
}


const char* MakeParser::SplitLineRegExp(char* intext,
	char* filename,char *text, int* linenumber,int *msgtype,int concatdir)
{
 int i=0,j,nmatch=4,tlen;
 regmatch_t pmatch[12];
 char *txt;
 int regextype,bufflen,cdirlen;

 strcpy(text,intext);

 if (text[i]=='@') {
 	 i++;
 	 while (text[i] && !(text[i-1]=='@' && text[i]=='.')) {
 		 i++;
 	 }
 	 if (text[i]=='.')
 		 i++;
 }
 j=i;
 txt=&text[i];
 bufflen=strlen(txt);
 cdirlen=strlen(mCurrDir);
 regextype=TestRegSet(txt,nmatch,pmatch);
 *linenumber=0;
 filename[0]='\0';
 if(regextype!=-1) {
      if(regextype<2 ) {
       /*entering or leaving directory */
 	  *msgtype=FLMAKE_MSG_PLAIN;
 	  for(i=0; i<nmatch && (pmatch[i].rm_so!=-1 && pmatch[i].rm_so<=bufflen
 			    && pmatch[i].rm_eo<=bufflen
 			    && pmatch[i].rm_eo>=pmatch[i].rm_so ) ; i++);
 	  if(i>=3) {
 	    *linenumber=0;
 	    i=2;
 	    strncpy(mCurrDir,txt+(pmatch[i].rm_so),pmatch[i].rm_eo-pmatch[i].rm_so);
 	    mCurrDir[pmatch[i].rm_eo-pmatch[i].rm_so]='/';
 	    mCurrDir[pmatch[i].rm_eo-pmatch[i].rm_so+1]='\0';
 	  }
      } else {
 	 /*it's an info, warning or error*/
 	  if(regextype<mINPreg)  {
 	     *msgtype=FLMAKE_MSG_INFO;
 	  } else if(regextype<mWNPreg)  {
             *msgtype=FLMAKE_MSG_WARNING;
 	  } else {
             *msgtype=FLMAKE_MSG_ERROR;
 	  }

 	  for(i=0; i<nmatch && (pmatch[i].rm_so!=-1 && pmatch[i].rm_so<=bufflen
 			    && pmatch[i].rm_eo<=bufflen
 			    && pmatch[i].rm_eo>=pmatch[i].rm_so ) ; i++);
 	  if(i>=3 ) {
 	    *linenumber=0;
 	    i=2;
 	    strncpy(filename,txt+(pmatch[i].rm_so),pmatch[i].rm_eo-pmatch[i].rm_so);
 	    filename[pmatch[i].rm_eo-pmatch[i].rm_so]='\0';
 	    *linenumber=atoi(filename);
 	    if(*linenumber!=0) {
 	      i=1;
 	      if(concatdir && ((txt+(pmatch[i].rm_so))[0])!='/' ) {
 		strcpy(filename,mCurrDir);
 		strncpy(filename+cdirlen,txt+(pmatch[i].rm_so),pmatch[i].rm_eo-pmatch[i].rm_so);
 		tlen=cdirlen+pmatch[i].rm_eo-pmatch[i].rm_so;
 		filename[tlen++]='\0';
 		strinstr(txt,filename,pmatch[i].rm_so,pmatch[i].rm_eo-pmatch[i].rm_so,tlen);
 	      } else {
 		strncpy(filename,txt+(pmatch[i].rm_so),pmatch[i].rm_eo-pmatch[i].rm_so);
 		filename[pmatch[i].rm_eo-pmatch[i].rm_so]='\0';
 	      }
 	    } else {
 	      i=1;
 	      strncpy(filename,txt+(pmatch[i].rm_so),pmatch[i].rm_eo-pmatch[i].rm_so);
 	      filename[pmatch[i].rm_eo-pmatch[i].rm_so]='\0';
 	      *linenumber=atoi(filename);
 	      if(*linenumber!=0) {
 		i=2;
 		if(concatdir && ((txt+(pmatch[i].rm_so))[0])!='/' ) {
 		  strcpy(filename,mCurrDir);
 		  strncpy(filename+cdirlen,txt+(pmatch[i].rm_so),pmatch[i].rm_eo-pmatch[i].rm_so);
 		  tlen=cdirlen+pmatch[i].rm_eo-pmatch[i].rm_so;
 		  filename[tlen++]='\0';
 		  strinstr(txt,filename,pmatch[i].rm_so,pmatch[i].rm_eo-pmatch[i].rm_so,tlen);
 		} else {
 		  strncpy(filename,txt+(pmatch[i].rm_so),pmatch[i].rm_eo-pmatch[i].rm_so);
 		  filename[pmatch[i].rm_eo-pmatch[i].rm_so]='\0';
 		}
 	      }
 	    }
 	  }
      }
 } else {
   *msgtype=FLMAKE_MSG_PLAIN;
 }

 return &text[j];
}

void MakeParser::SelectErrorLine(void) 
{ 
	int msgtype;
	int line=mErrorBrowser->value();
	if (line>=mErrorBrowser->size()) return;
	const char* text=mErrorBrowser->text(line);
	if (text==0) return;
	char filename[1024];
	char resptxt[1024];
	int linenumber=0;
	SplitLineRegExp((char*)text,filename,resptxt,&linenumber,&msgtype,0);
	char cmd[1024];
	if (linenumber) {
	    if(msgtype) {
		snprintf(cmd,1024,mInvEditor,linenumber,filename);
		system(cmd);
	     }
	}else{
	     char* filenameptr=LookForArg(" -c ",text);
	     if (filenameptr) {
		snprintf(cmd,1024,mInvEditor,0,filenameptr);
		system(cmd);
		delete [] filenameptr;
	     }
	}
}

char *strinstr(char *srcdest,char *newtxt, size_t posinsrc,size_t nschars,size_t ntlen) 
{
	char *tmpstr;
	tmpstr=strdup(&srcdest[posinsrc+nschars]);
	strncpy(&srcdest[posinsrc],newtxt,ntlen);
	strcat(&srcdest[posinsrc+ntlen-1],tmpstr);
	free(tmpstr);
	return srcdest;
}

void MakeParser::InitMessageHeaders(void) 
{
	char *ts;
	/* Initializes the Message color/types for the different message groups that
	 can be displayed by flmake.
	 TODO: this defaults should be changed with a dialog and saved in the
	 prefs file. */
	mNumMsgHeaders=0;
	mMsgHeaders[FLMAKE_MSG_PLAIN]=strdup("@C0@.");
	mNumMsgHeaders++;
	mMsgHeaders[FLMAKE_MSG_ISTATUS]=strdup("@B3@.");
	mNumMsgHeaders++;
	mMsgHeaders[FLMAKE_MSG_IERROR]=strdup("@B1@C7@.");
	mNumMsgHeaders++;
	mMsgHeaders[FLMAKE_MSG_WARNING]=strdup("@C4@.");
	mNumMsgHeaders++;
	mMsgHeaders[FLMAKE_MSG_ERROR]=strdup("@C1@.");
	mNumMsgHeaders++;
	mMsgHeaders[FLMAKE_MSG_INFO]=strdup("@C11@.");
	mNumMsgHeaders++;
	mMsgHeaders[FLMAKE_MSG_RUNTARGET0]=strdup("@C0@.");
	mNumMsgHeaders++;
	mMsgHeaders[FLMAKE_MSG_RUNTARGET1]=strdup("@C0@B53@.");
	mNumMsgHeaders++;

	/* set the error and warning arrows to the same color as the corresponding
	 message text color */
	if((ts=strstr(mMsgHeaders[FLMAKE_MSG_ERROR],"@C"))!=NULL) {
	  mNextError->labelcolor(atoi(ts+2));
	  mPrevError->labelcolor(atoi(ts+2));
	}
	if((ts=strstr(mMsgHeaders[FLMAKE_MSG_WARNING],"@C"))!=NULL) {
	  mNextWarning->labelcolor(atoi(ts+2));
	  mPrevWarning->labelcolor(atoi(ts+2));
	}
	mMainGroup->labelcolor(FL_BLACK);
}

void MakeParser::Init(int argc, char **argv) 
{
	strcpy(mMakeArgs,"");
	for(int i=1;i <argc; i++) {
		strcat(mMakeArgs,argv[i]);
		strcat(mMakeArgs," ");
	}
	strcpy(mCurrDir,"");
	strcpy(mInvEditor,"nc -noask -line %d %s");

	/* this is true if the user want flmake to open the file when using the
	  error/warning browse arrows or if it must click to open the file.
	  by deafult is must click to open==false
	  TODO: put this an options dialog and add to the preferences file.
	*/
	mEditOnMsgBrowse=false;

	/* clear the run targets browser control toggle */
	mCleanRunTargets=true;

	InitMessageHeaders();
	InitRegex();

	fl_message_font(mRunTargetsBrowser->textfont(),
	                mRunTargetsBrowser->textsize());
}


