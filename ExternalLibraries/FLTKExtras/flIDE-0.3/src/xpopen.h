/**
 Extended popen, pclose commands to easly get the forked child pid for fine
 tuned control of the pipe, based on the standard popen, pclose implementation 
 of the gnu glibc2, it consist of the following set of functions:
 
 commands using standard types:
 xpopen, xpclose, xpkillclose
 
 commands using the XPopenDescr pipe file descriptor definded in this file:
 xpopen_XPD, xpclose_XPD, xpkillclose_XPD  
*/
#ifndef __XPOPEN_H_
#define __XPOPEN_H_
#include <errno.h>
#include <stddef.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <fcntl.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Basic eXtended popen, pclose pkillclose basic functions */
FILE *xpopen(const char *command,const char *mode,pid_t *pid);
int xpclose (FILE *stream, pid_t pid);
int xpkillclose (FILE *stream, pid_t pid);

/*--------------------------------------------------------------------------*/
/* Useful eXtended xpopen, xpclose xpkillclose functions to use a new popened 
 file  descriptor struct to ease & standarize use in threaded programs  */

/* Popened file descriptor type/struct */
typedef struct __XPopenDescrStruct {
  FILE *fd;    /* pipe file descriptor */
  pid_t cpid;  /* child process id */
  char mode;   /* io mode (r/w) */
} XPopenDescr;

/* eXtended xpopen, xpclose xpkillclos functions to return a XPopenDescr
 (that's why the suffix (XPD))*/
XPopenDescr *xpopen_XPD(const char *command,const char *mode);
int xpclose_XPD (XPopenDescr *pdescr);
int xpkillclose_XPD (XPopenDescr *pdescr);

#ifdef __cplusplus
}
#endif

#endif /* __XPOPEN_H_ */
