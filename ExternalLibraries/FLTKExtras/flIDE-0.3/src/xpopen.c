/**
 Extended popen, pclose commands to easly get the forked child pid for fine
 tuned control of the pipe, based on the standard popen, pclose implementation 
 of the gnu glibc2, it consist of the following set of functions:
 
 commands using standard types:
 xpopen, xpclose, xpkillclose
 
 commands using the XPopenDescr pipe file descriptor definded in the .h file:
 xpopen_XPD, xpclose_XPD, xpkillclose_XPD  
*/
#include <errno.h>
#include <stddef.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <fcntl.h>
#include "xpopen.h"

#define SH_PATH "/bin/sh"	/* Shell to run.  */
#define SH_NAME "sh"		/* Name to give it.  */
/*--------------------------------------------------------------------------*/
/* Basic Extended xpopen, pclose pkillclose basic functions   */

/**  Extended popen command.
 Same functionality as the standard popen command but returns also the created
 child process id in *pid..
*/
FILE *xpopen(const char *command,const char *mode,pid_t *pid)
{
 int pipedes[2];
 FILE *stream;
 pid_t cpid;
 
 /* verify parameters */
 if(command==NULL || mode == NULL || pid == NULL || (*mode!='r' && *mode!='w')) {
    return NULL;
 }   
 /* create the pipe */
 if(pipe(pipedes) <0) 
   return NULL;

 /* fork off the child */
 cpid=vfork();
 *pid=cpid;
 if(cpid==(pid_t) -1) {
   /*the fork failed */
   (void)close(pipedes[0]);
   (void)close(pipedes[1]);
   return NULL;
 } else if (cpid == (pid_t) 0) {
    /* We are in the child side.
       Make the write side of the pipe be the stdin or the read side the
       stdout*/
    const char *new_argv[4];
    if ((*mode == 'w' ? dup2(pipedes[STDIN_FILENO], STDIN_FILENO) :
           dup2 (pipedes[STDOUT_FILENO], STDOUT_FILENO)) < 0) {
        _exit (127);
    }    
    /* close the pipe descriptors */
    (void)close(pipedes[STDIN_FILENO]);
    (void)close(pipedes[STDOUT_FILENO]);
    /*exec the shell */
    new_argv[0]=SH_NAME;
    new_argv[1]="-c";
    new_argv[2]=command;
    new_argv[3]=NULL;
    (void) execv (SH_PATH, (char *const *) new_argv);
    /*die if it failed*/  
    _exit (127);
     
 }
 /* here we are on the parent side */
 
 /* close the irrelevan side of the pipe and open the relevan one as a new stream
    mark our side of the pipe to close on exec, so the new children won't be
    able to see it */
 if(*mode=='r') {
   (void)close(pipedes[STDOUT_FILENO]);
   (void) fcntl (pipedes[STDIN_FILENO], F_SETFD, FD_CLOEXEC);
   stream = fdopen (pipedes[STDIN_FILENO], mode);
 } else {
   (void) close (pipedes[STDIN_FILENO]);
   (void) fcntl (pipedes[STDOUT_FILENO], F_SETFD, FD_CLOEXEC);
   stream = fdopen (pipedes[STDOUT_FILENO], mode);
 }
 if(stream==NULL) {
  /*the stream couldn't be opened*/
/*  int save = errno;*/
  (void)kill(cpid,SIGKILL);
  (void) close (pipedes[*mode == 'r' ? STDOUT_FILENO : STDIN_FILENO]);
#ifndef NO_WAITPID
  (void) waitpid (cpid, (int *) NULL, 0);
#else
  {
   pid_t dead;
   do
    dead = wait ((int *) NULL);
   while (dead > 0 && dead != cpid);
  }
 #endif
/*   __set_errno (save); */
   return NULL;  
 }
 return stream;
}

/**  Extended pclose command with child process id.
 Same functionality as the standard pclose command, but also requires the
 xpopened child id, complements the xpopen command.
*/
int xpclose (FILE *stream, pid_t pid)
{
  pid_t dead;
  int status=pid;
  if(stream == NULL || pid <0 ) {
    return -1;
  }
#ifndef NO_WAITPID
  dead= waitpid (pid, (int *) NULL, 0);
#else
   do
    dead = wait ((int *) NULL);
   while (dead > 0 && dead != pid);
#endif  
 if(dead!=pid)
   status=-1;
  if(fclose(stream))
    return -1;
   
 return status;
}

/** Augmented xpclose command to kill the child process, always before closing
 the pipe.  
 Same functionality as the xpclose command, but also requires the
 xpopened child id so to be able to kill it if it's still running..
*/
int xpkillclose (FILE *stream, pid_t pid)
{
  if(stream == NULL || pid <0 ) {
    return -1;
  }
  (void)kill(pid,SIGKILL);
  return xpclose(stream,pid);
}

/*--------------------------------------------------------------------------*/
/* Useful eXtended xpopen, xpclose xpkillclose functions to use a new popened 
 file  descriptor struct to ease & standarize use in threaded programs  */

/**  Extended xpopen command.
 Same functionality as the xpopen command but returns a XPopenDescr pipe file
 descriptor
*/
XPopenDescr *xpopen_XPD(const char *command,const char *mode)
{
 XPopenDescr *presp;
 presp=(XPopenDescr*)malloc(sizeof(XPopenDescr));
 if(((*presp).fd=xpopen(command,mode,&((*presp).cpid)))==NULL)  {
   free(presp);
   return NULL;
 }      
 (*presp).mode=*mode;
 return presp;    
}

/**  Extended xpclose command.
 Same functionality as the xpclose command but returns a XPopenDescr pipe file
 descriptor
*/
int xpclose_XPD (XPopenDescr *pdescr)
{
 int status;
 if(pdescr==NULL)
   return -1;
 if((status=xpclose((*pdescr).fd, (*pdescr).cpid))!=-1) {
    free(pdescr);
    return status;
 }
 return -1;
}

/**  Extended xpkillclose command.
 Same functionality as the xpkillclose command but returns a XPopenDescr pipe 
 file descriptor
*/
int xpkillclose_XPD (XPopenDescr *pdescr)
{
 int status;
 if(pdescr==NULL)
   return -1;
 if((status=xpkillclose((*pdescr).fd, (*pdescr).cpid))!=-1) {
    free(pdescr);
    return status;
 } 
 return -1;
}
